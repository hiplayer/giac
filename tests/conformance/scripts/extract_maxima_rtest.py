#!/usr/bin/env python3
"""
Extract limit / integrate / diff / solve / taylor cases from Maxima rtest .mac
files into giac-rs conformance fixture JSON (fixture-v1).

Spec: .doc/external-test-resources.md §5
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Iterator

EXTRACTOR_VERSION = "1"

DOMAIN_DEFAULT_VERIFY: dict[str, str] = {
    "limit": "sympy_limit",
    "integrate": "integrate_derivative",
    "diff": "diff_inverse",
    "solve": "solve_residual",
    "taylor": "sympy_limit",  # series compare — extend sympy_verify later
    "series": "sympy_limit",
    "partfrac": "partfrac_expand",
    "desolve": "desolve_odesol",
    "poly": "literal",
    "trig": "literal",
    "other": "parse_only",
}

# Maxima builtins we translate to giac calls in `line`.
CALL_RE = re.compile(
    r"\b(limit|integrate|diff|solve|taylor|partfrac|desolve)\s*\(",
    re.IGNORECASE,
)

CK_INT_RE = re.compile(r"CK-INT-(\d+)", re.IGNORECASE)


@dataclass
class ExtractedCase:
    upstream_input: str
    expected_hint: str | None
    source_line: int
    comment: str | None = None
    extra_tags: list[str] = field(default_factory=list)


def maxima_to_giac(expr: str) -> str:
    """Best-effort Maxima → giac surface syntax."""
    s = expr.strip()
    if s.endswith(";"):
        s = s[:-1].strip()
    if s.endswith("$"):
        s = s[:-1].strip()

    # log → ln (giac)
    s = re.sub(r"\blog\b", "ln", s)

    replacements = [
        (r"\b%pi\b", "pi"),
        (r"\b%e\b", "exp(1)"),
        (r"\b\+?inf\b", "+infinity"),
        (r"\binfinity\b", "+infinity"),
        (r"\bminf\b", "+infinity"),
        (r"\b\+infinity\b", "+infinity"),
        (r"\b-inf\b", "-infinity"),
        (r"\bminus_infinity\b", "-infinity"),
    ]
    for pat, repl in replacements:
        s = re.sub(pat, repl, s, flags=re.IGNORECASE)

    # taylor → series for giac (configurable later)
    s = re.sub(r"\btaylor\b", "series", s, flags=re.IGNORECASE)

    return s


def call_kind(expr: str) -> str | None:
    m = CALL_RE.search(expr)
    if not m:
        return None
    name = m.group(1).lower()
    return "taylor" if name == "taylor" else name


def detect_domain(expr: str, default: str) -> str:
    kind = call_kind(expr)
    if kind is None:
        return default
    if kind == "taylor":
        return "taylor"
    return kind if kind in DOMAIN_DEFAULT_VERIFY else default


def build_line(expr: str, domain: str) -> str | None:
    g = maxima_to_giac(expr)
    if not CALL_RE.search(g):
        return None
    return g


def tags_from_context(comment: str | None, expr: str) -> list[str]:
    tags: list[str] = []
    blob = " ".join(filter(None, [comment, expr]))
    if re.search(r"gruntz", blob, re.I):
        tags.append("gruntz")
    if re.search(r"\+infinity|-infinity|\binf\b", blob, re.I):
        tags.append("infinity")
    m = CK_INT_RE.search(blob)
    if m:
        tags.append(f"CK-INT-{m.group(1)}")
    return tags


def strip_maxima_comment(line: str) -> tuple[str, str | None]:
    """Return (code, trailing comment)."""
    if "/*" in line:
        return line.split("/*", 1)[0].strip(), line
    if ";" in line and not line.strip().startswith("/*"):
        parts = line.split(";", 1)
        if len(parts) == 2 and parts[1].strip().startswith("/*"):
            return parts[0].strip() + ";", parts[1].strip()
    return line.strip(), None


def iter_cases(path: Path, verbose: bool) -> Iterator[ExtractedCase]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    pending_comment: str | None = None
    pending_expr: tuple[str, int] | None = None

    for i, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("/*") and line.endswith("*/"):
            if "CK-INT" in line or "gruntz" in line.lower():
                pending_comment = line
            continue

        if line.startswith("/*"):
            pending_comment = line
            continue

        code, inline_comment = strip_maxima_comment(line)
        comment = inline_comment or pending_comment
        pending_comment = None

        if not code or code.startswith("--"):
            continue

        # Skip setup / teardown
        if re.match(r"^(block|kill|reset|load|batch)\b", code, re.I):
            continue

        if CALL_RE.search(code):
            if pending_expr is not None:
                if verbose:
                    print(
                        f"{path}:{pending_expr[1]}: warning: no expected for {pending_expr[0]!r}",
                        file=sys.stderr,
                    )
            pending_expr = (code, i)
            continue

        # Expected line (no CAS call)
        if pending_expr is not None:
            expected = maxima_to_giac(code)
            if expected.lower() in ("true", "false", "error"):
                pending_expr = None
                continue
            yield ExtractedCase(
                upstream_input=pending_expr[0],
                expected_hint=expected,
                source_line=pending_expr[1],
                comment=comment,
                extra_tags=tags_from_context(comment, pending_expr[0]),
            )
            pending_expr = None

    if pending_expr is not None and verbose:
        print(
            f"{path}:{pending_expr[1]}: warning: trailing expr without expected",
            file=sys.stderr,
        )


def merge_entries(
    cases: list[ExtractedCase],
    *,
    path: Path,
    domain: str,
    prefix: str,
    start_id: int,
    tag: str | None,
) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    n = start_id
    for case in cases:
        case_domain = detect_domain(case.upstream_input, domain)
        if case_domain != domain:
            continue
        line = build_line(case.upstream_input, domain)
        if line is None:
            continue
        tags = list(case.extra_tags)
        if tag and tag not in tags:
            tags.append(tag)
        entry: dict[str, Any] = {
            "id": f"{prefix}-{n:03d}",
            "line": line,
            "enabled": False,
            "api": case_domain,
            "verify": DOMAIN_DEFAULT_VERIFY.get(case_domain, "parse_only"),
            "source_ref": f"{path.name}:{case.source_line}",
            "upstream_input": case.upstream_input.rstrip(";").strip(),
            "maxima_input": case.upstream_input.rstrip(";").strip(),
            "source_file": path.name,
            "source_line": case.source_line,
        }
        if case.expected_hint:
            entry["expected_hint"] = case.expected_hint
            entry["maxima_expected"] = case.expected_hint
        if tags:
            entry["tags"] = tags
        if case.comment:
            entry["notes"] = case.comment[:200]
        entries.append(entry)
        n += 1
    return entries


def build_fixture(
    entries: list[dict[str, Any]],
    *,
    path: Path,
    domain: str,
    verify: str | None,
    source_url: str | None,
) -> dict[str, Any]:
    return {
        "version": 1,
        "source": {
            "name": "maxima",
            "file": str(path),
            "url": source_url or "",
            "license": "GPL-2.0-or-later",
            "extracted_at": date.today().isoformat(),
            "extractor": "extract_maxima_rtest.py",
            "extractor_version": EXTRACTOR_VERSION,
        },
        "domain": domain,
        "verify": verify or DOMAIN_DEFAULT_VERIFY.get(domain, "parse_only"),
        "entries": entries,
    }


def load_manifest(manifest_path: Path) -> list[dict[str, Any]]:
    try:
        import yaml  # type: ignore
    except ImportError as e:
        raise SystemExit(
            "manifest mode requires PyYAML: pip install pyyaml"
        ) from e
    data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise SystemExit(f"manifest must be a YAML list: {manifest_path}")
    return data


def run_extract(
    input_path: Path,
    *,
    domain: str,
    prefix: str,
    start_id: int,
    tag: str | None,
    verify: str | None,
    source_url: str | None,
    verbose: bool,
) -> dict[str, Any]:
    cases = list(iter_cases(input_path, verbose))
    entries = merge_entries(
        cases,
        path=input_path,
        domain=domain,
        prefix=prefix,
        start_id=start_id,
        tag=tag,
    )
    return build_fixture(
        entries,
        path=input_path,
        domain=domain,
        verify=verify,
        source_url=source_url,
    )


def default_prefix(domain: str) -> str:
    return {
        "limit": "LIM-M",
        "integrate": "INT-M",
        "diff": "DIF-M",
        "solve": "SOL-M",
        "taylor": "SER-M",
        "series": "SER-M",
    }.get(domain, "MAX-M")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract Maxima rtest .mac → giac-rs fixture-v1 JSON"
    )
    parser.add_argument("--input", "-i", action="append", type=Path, help="Maxima .mac file")
    parser.add_argument("--manifest", "-m", type=Path, help="YAML batch manifest")
    parser.add_argument("--domain", "-d", default="limit", choices=list(DOMAIN_DEFAULT_VERIFY))
    parser.add_argument("--output", "-o", type=Path, help="Output JSON path")
    parser.add_argument("--prefix", "-p", help="Entry id prefix (default from domain)")
    parser.add_argument("--verify", choices=list(DOMAIN_DEFAULT_VERIFY.values()))
    parser.add_argument("--tag", help="Tag added to every entry")
    parser.add_argument("--start-id", type=int, default=1)
    parser.add_argument("--source-url", help="Upstream URL stored in fixture source")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    if not args.input and not args.manifest:
        parser.error("one of --input or --manifest is required")

    jobs: list[dict[str, Any]] = []
    if args.manifest:
        for item in load_manifest(args.manifest):
            jobs.append(item)
    if args.input:
        for p in args.input:
            jobs.append(
                {
                    "input": str(p),
                    "domain": args.domain,
                    "output": str(args.output) if args.output and len(args.input) == 1 else None,
                    "prefix": args.prefix,
                    "tag": args.tag,
                    "verify": args.verify,
                    "source_url": args.source_url,
                    "start_id": args.start_id,
                }
            )

    if len(jobs) > 1 and args.output and not args.manifest:
        parser.error("--output with multiple --input requires --manifest")

    total_entries = 0
    for job in jobs:
        input_path = Path(job["input"])
        if not input_path.is_file():
            print(f"error: not found: {input_path}", file=sys.stderr)
            return 1
        domain = job.get("domain", args.domain)
        prefix = job.get("prefix") or args.prefix or default_prefix(domain)
        fixture = run_extract(
            input_path,
            domain=domain,
            prefix=prefix,
            start_id=int(job.get("start_id", args.start_id)),
            tag=job.get("tag") or args.tag,
            verify=job.get("verify") or args.verify,
            source_url=job.get("source_url") or args.source_url,
            verbose=args.verbose,
        )
        total_entries += len(fixture["entries"])
        out = Path(job["output"]) if job.get("output") else args.output
        if args.dry_run:
            print(json.dumps(fixture, indent=2, ensure_ascii=False))
        elif out:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(
                json.dumps(fixture, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            print(f"wrote {len(fixture['entries'])} entries → {out}", file=sys.stderr)
        else:
            parser.error("--output required unless --dry-run")

    if total_entries == 0:
        print("warning: no entries extracted", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
