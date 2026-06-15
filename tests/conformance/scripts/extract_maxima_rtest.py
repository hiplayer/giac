#!/usr/bin/env python3
"""Extract integrate/limit/solve cases from Maxima rtest .mac files.

Reads Maxima regression tests (`input;` + `expected$` pairs) and optional
standalone API lines, translates Maxima dialect to giac/Xcas-style input,
and writes a JSON fixture draft for giac-conformance.

Usage:
  extract_maxima_rtest.py rtest_limit_wester.mac -o /tmp/limit.json
  extract_maxima_rtest.py ../maxima-5.49.0/tests/rtest_limit_wester.mac \\
      ../maxima-5.49.0/tests/wester_problems/test_equations.mac \\
      --apis integrate,limit,solve --include-standalone

Regenerate committed fixture (from giac-rs/tests/conformance):

  python3 scripts/extract_maxima_rtest.py \\
    ../../../maxima-5.49.0/tests/rtest_limit_wester.mac \\
    -o fixtures/phase4_maxima_rtest.json --stats

Output entries default to ``enabled: false``; review and flip after probing
with giac-cli / ``run_line`` + ``sympy_verify.py``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Literal

ApiName = Literal["integrate", "limit", "solve", "diff"]

SUPPORTED_APIS = ("integrate", "limit", "solve", "diff")

# Maxima-only constructs we cannot translate or verify yet.
SKIP_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bassume\b"), "assume"),
    (re.compile(r"\bforget\b"), "forget"),
    (re.compile(r"\bblock\s*\("), "block"),
    (re.compile(r"\bev\s*\("), "ev"),
    (re.compile(r"\bmatchdeclare\b"), "matchdeclare"),
    (re.compile(r"\berrcatch\b"), "errcatch"),
    (re.compile(r"\bunit_step\b"), "unit_step"),
    (re.compile(r"\bgamma\s*\("), "gamma"),
    (re.compile(r"\berf[i]?\s*\("), "erf"),
    (re.compile(r"\b%s\b"), "string_format"),
    (re.compile(r"\bnounify\b"), "nounify"),
    (re.compile(r"\bclosedform\b"), "closedform"),
    (re.compile(r"\bantidiff\b"), "antidiff"),
    (re.compile(r"\b'integrate\b"), "quoted_integrate"),
    (re.compile(r"\b'limit\b"), "quoted_limit"),
    (re.compile(r"\bintegrate\s*\(\s*%"), "integrate_previous"),
    (re.compile(r"\b[a-zA-Z_][\w]*\s*\(\s*x\s*\)\s*:="), "function_definition"),
    (re.compile(r"\bdepends\b"), "depends"),
    (re.compile(r"\bremfunction\b"), "remfunction"),
    (re.compile(r"\bload\s*\("), "load"),
    (re.compile(r"\bdeclare\b"), "declare"),
    (re.compile(r"\bremove\b"), "remove"),
    (re.compile(r"\bmultiplicities\b"), "multiplicities"),
    (re.compile(r",\s*(?:minus|plus)\s*\)"), "limit_direction"),
    (re.compile(r"\b(?:und|ind|pinf|ninf)\b"), "limit_und"),
    (re.compile(r"\b%th\b"), "previous_output"),
    (re.compile(r"\b%piargs\b"), "piargs"),
    (re.compile(r"\b%iargs\b"), "iargs"),
    (re.compile(r"\bii\b"), "maxima_imag_unit"),
    (re.compile(r"\btrue\b|\bfalse\b"), "boolean_literal"),
    (re.compile(r"\bn!"), "factorial"),
    (re.compile(r"\bgcd\s*\("), "gcd_call"),
    (re.compile(r"\bdiffcheck\b"), "diffcheck"),
    (re.compile(r"\beqn\d+\b"), "wester_eqn_ref"),
]

def strip_mac_comment(line: str) -> str:
    if "/*" in line:
        return line.split("/*", 1)[0].strip()
    return line.strip()

API_TOP_RE = re.compile(
    r"^(?P<api>integrate|limit|solve|diff)\s*\(", re.IGNORECASE
)


@dataclass
class ExtractedCase:
    api: ApiName
    maxima_input: str
    line: str
    maxima_expected: str | None
    source_file: str
    source_line: int
    verify: str
    enabled: bool = False
    skip_reason: str | None = None
    maxima_postfix: list[str] = field(default_factory=list)
    id: str = ""

    def to_json(self) -> dict:
        out = {
            "id": self.id,
            "line": self.line,
            "api": self.api,
            "verify": self.verify,
            "enabled": self.enabled,
            "maxima_input": self.maxima_input,
            "maxima_expected": self.maxima_expected,
            "source_file": self.source_file,
            "source_line": self.source_line,
            "skip_reason": self.skip_reason,
        }
        if self.maxima_postfix:
            out["maxima_postfix"] = self.maxima_postfix
        return out


@dataclass
class ExtractStats:
    pairs_seen: int = 0
    standalone_seen: int = 0
    extracted: int = 0
    skipped_api: int = 0
    skipped_pattern: int = 0
    skipped_translate: int = 0
    by_api: dict[str, int] = field(default_factory=dict)


MAXIMA_POSTFIX_FLAGS = frozenset(
    {
        "radcan",
        "factor",
        "ratsimp",
        "ratexpand",
        "expand",
        "fullratsimp",
        "trigsimp",
        "logexpand",
        "logcontract",
    }
)


def strip_maxima_postfix(expr: str) -> tuple[str, list[str]]:
    """Drop trailing `, flag` Maxima modifiers (e.g. integrate(...), radcan)."""
    depth = 0
    for i, ch in enumerate(expr):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                rest = expr[i + 1 :].strip()
                if rest.startswith(","):
                    flags = [f.strip() for f in rest[1:].split(",") if f.strip()]
                    return expr[: i + 1].strip(), flags
                return expr.strip(), []
    return expr.strip(), []


def has_unsupported_postfix(flags: list[str]) -> bool:
    for flag in flags:
        low = flag.lower()
        if ":" in flag:
            return True
        if low in MAXIMA_POSTFIX_FLAGS:
            return True
    return False


def strip_mac_comment(line: str) -> str:
    if "/*" in line:
        return line.split("/*", 1)[0].strip()
    return line.strip()


def skip_reason(expr: str) -> str | None:
    for pat, name in SKIP_PATTERNS:
        if pat.search(expr):
            return name
    return None


def maxima_to_giac(expr: str) -> str:
    """Best-effort Maxima → giac/Xcas expression translation."""
    s = expr.strip()
    s = re.sub(r"\s+", " ", s)
    # trailing option clauses: integrate(...), domain:real
    s = re.sub(r",\s*domain\s*:\s*\w+", "", s, flags=re.IGNORECASE)
    s = s.replace("**", "^")
    s = re.sub(r"\blog\b", "ln", s)
    s = re.sub(r"\binf\b", "+infinity", s)
    s = re.sub(r"\bminf\b", "-infinity", s)
    s = re.sub(r"%e\^", "exp(1)^", s)
    s = re.sub(r"%e\b", "exp(1)", s)
    s = re.sub(r"%i\b", "i", s)
    s = re.sub(r"%pi\b", "pi", s)
    s = re.sub(r"\bpi\b", "pi", s)
    # Maxima 'diff is quoted derivative — rare in our targets
    s = s.replace("'diff", "diff")
    s = s.replace("'integrate", "integrate")
    s = s.replace("'limit", "limit")
    return s


def detect_api(expr: str) -> ApiName | None:
    m = API_TOP_RE.match(expr.strip())
    if not m:
        return None
    api = m.group("api").lower()
    if api == "int":
        return "integrate"
    return api  # type: ignore[return-value]


def normalize_giac_line(maxima_expr: str) -> str:
    """Produce a giac conformance line (no trailing semicolon — matches fixtures)."""
    return maxima_to_giac(maxima_expr).rstrip(";")


def suggest_verify(api: ApiName, giac_line: str) -> str:
    if api == "integrate":
        if re.search(r"integrate\([^,]+,x,[^,]+,[^)]+\)", giac_line):
            return "integrate_definite"
        return "integrate_derivative"
    if api == "limit":
        return "sympy_limit"
    if api == "solve":
        return "sympy_solve"
    if api == "diff":
        return "sympy_diff"
    return "sympy_equiv"


def parse_rtest_pairs(path: Path) -> Iterator[tuple[str, str | None, int]]:
    """Yield (input_expr, expected_or_none, line_no) from rtest .mac files."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        raw = strip_mac_comment(lines[i])
        if not raw or raw in ("0", "0$", "true$", "false$", "done$", "done;"):
            i += 1
            continue
        if raw.endswith(";"):
            inp = raw[:-1].strip()
            expected: str | None = None
            if i + 1 < len(lines):
                nxt = strip_mac_comment(lines[i + 1])
                if nxt.endswith("$"):
                    expected = nxt[:-1].strip()
                    i += 2
                    yield inp, expected, i - 1
                    continue
            yield inp, None, i + 1
            i += 1
            continue
        if raw.endswith("$") and not raw.endswith(";"):
            # orphan expected line
            i += 1
            continue
        i += 1


def extract_from_file(
    path: Path,
    apis: set[ApiName],
    include_standalone: bool,
    id_prefix: str,
    stats: ExtractStats,
) -> list[ExtractedCase]:
    cases: list[ExtractedCase] = []
    seq = 0

    for inp, expected, line_no in parse_rtest_pairs(path):
        if expected is not None:
            stats.pairs_seen += 1
        else:
            if include_standalone:
                stats.standalone_seen += 1
            else:
                continue

        api = detect_api(inp)
        if api is None or api not in apis:
            stats.skipped_api += 1
            continue

        core_inp, postfix_flags = strip_maxima_postfix(inp)
        if postfix_flags and has_unsupported_postfix(postfix_flags):
            stats.skipped_pattern += 1
            continue

        reason = skip_reason(core_inp)
        if reason:
            stats.skipped_pattern += 1
            continue
        if expected is not None and skip_reason(expected or ""):
            stats.skipped_pattern += 1
            continue

        try:
            line = normalize_giac_line(core_inp)
        except Exception:
            stats.skipped_translate += 1
            continue

        if skip_reason(line):
            stats.skipped_pattern += 1
            continue

        seq += 1
        case_id = f"{id_prefix}-{api[:3].upper()}-{seq:03d}"
        case = ExtractedCase(
            id=case_id,
            api=api,
            maxima_input=inp,
            line=line.rstrip(";"),
            maxima_expected=expected,
            source_file=path.name,
            source_line=line_no,
            verify=suggest_verify(api, line),
            enabled=False,
            skip_reason=None,
            maxima_postfix=postfix_flags,
        )
        cases.append(case)
        stats.extracted += 1
        stats.by_api[api] = stats.by_api.get(api, 0) + 1

    return cases


def default_id_prefix(path: Path) -> str:
    stem = path.stem.upper().replace("RTEST_", "").replace("TEST_", "")
    stem = re.sub(r"[^A-Z0-9]+", "_", stem).strip("_")
    return f"MX_{stem[:20]}"


def build_fixture(
    paths: Iterable[Path],
    apis: set[ApiName],
    include_standalone: bool,
) -> tuple[dict, ExtractStats]:
    stats = ExtractStats()
    all_cases: list[ExtractedCase] = []

    for path in paths:
        prefix = default_id_prefix(path)
        all_cases.extend(
            extract_from_file(path, apis, include_standalone, prefix, stats)
        )

    # stable unique ids across files
    seen: set[str] = set()
    for case in all_cases:
        base = case.id
        n = 1
        while case.id in seen:
            n += 1
            case.id = f"{base}_{n}"
        seen.add(case.id)

    fixture = {
        "version": 1,
        "generator": "extract_maxima_rtest.py",
        "note": (
            "Draft fixture: entries are disabled by default. "
            "Review giac line translation, then enable and wire into conformance."
        ),
        "entries": [c.to_json() for c in all_cases],
    }
    return fixture, stats


def parse_apis(s: str) -> set[ApiName]:
    parts = [p.strip().lower() for p in s.split(",") if p.strip()]
    out: set[ApiName] = set()
    for p in parts:
        if p == "int":
            p = "integrate"
        if p not in SUPPORTED_APIS:
            raise argparse.ArgumentTypeError(f"unknown api: {p}")
        out.add(p)  # type: ignore[arg-type]
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mac_files",
        nargs="+",
        type=Path,
        help="Maxima .mac rtest files to scan",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write JSON fixture (default: stdout)",
    )
    parser.add_argument(
        "--apis",
        type=parse_apis,
        default=parse_apis("integrate,limit,solve"),
        help="Comma-separated APIs to extract (default: integrate,limit,solve)",
    )
    parser.add_argument(
        "--include-standalone",
        action="store_true",
        help="Include API lines without a following expected$ line",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print extraction stats to stderr",
    )
    args = parser.parse_args(argv)

    paths = [p.resolve() for p in args.mac_files]
    for p in paths:
        if not p.is_file():
            print(f"error: not a file: {p}", file=sys.stderr)
            return 1

    fixture, stats = build_fixture(paths, args.apis, args.include_standalone)
    text = json.dumps(fixture, indent=2, ensure_ascii=False) + "\n"

    if args.output:
        args.output.write_text(text, encoding="utf-8")
        print(f"wrote {len(fixture['entries'])} entries → {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(text)

    if args.stats:
        print(
            f"pairs={stats.pairs_seen} standalone={stats.standalone_seen} "
            f"extracted={stats.extracted} skip_api={stats.skipped_api} "
            f"skip_pattern={stats.skipped_pattern} by_api={stats.by_api}",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
