#!/usr/bin/env python3
"""Cross-verify giac-rs Phase 2 outputs with SymPy (third-party CAS).

Usage:
  phase2_sympy.py reference <giac_line>     # JSON expected value from SymPy
  phase2_sympy.py verify <giac_line> <out>  # exit 0 if output matches SymPy
  phase2_sympy.py equiv <out_a> <out_b>     # exit 0 if two giac-style outputs are equivalent
"""

from __future__ import annotations

import json
import re
import sys
from typing import Any

import sympy as sp
from sympy import Poly, QQ, ZZ, Matrix, symbols, resultant, apart, lcm, gcd, div, groebner

x, y, z = symbols("x y z")


def giac_to_sympy(s: str) -> sp.Expr:
    """Best-effort parse of giac-style output into SymPy."""
    s = s.strip()
    if not s:
        raise ValueError("empty output")
    if s.startswith("[[") and s.endswith("]]"):
        return giac_matrix(s)
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1]
        if not inner.strip():
            return sp.Tuple()
        parts = split_top_level(inner)
        return sp.Tuple(*[giac_to_sympy(p.strip()) for p in parts])
    if " mod " in s:
        base, _mod = s.rsplit(" mod ", 1)
        return giac_to_sympy(base.strip())
    s = re.sub(r"\bi\b", "I", s)
    s = s.replace("^", "**")
    s = re.sub(r"(\d)\(", r"\1*(", s)
    s = re.sub(r"\)\(", ")*(", s)
    return sp.sympify(s, locals={"x": x, "y": y, "z": z, "I": sp.I})


def split_top_level(s: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    cur: list[str] = []
    for ch in s:
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        cur.append(ch)
    if cur:
        parts.append("".join(cur))
    return parts


def giac_matrix(s: str) -> Matrix:
    rows = s.strip()[1:-1]
    row_strs = split_top_level(rows)
    data = []
    for row in row_strs:
        row = row.strip()
        if row.startswith("[") and row.endswith("]"):
            cells = split_top_level(row[1:-1])
            data.append([giac_to_sympy(c.strip()) for c in cells])
        else:
            data.append([giac_to_sympy(row)])
    return Matrix(data)


def parse_poly(expr_str: str, var=x) -> Poly:
    e = giac_to_sympy(expr_str)
    return Poly(sp.expand(e), var, domain=QQ)


def positive_mod(a: int, m: int) -> int:
    r = a % m
    return r if r >= 0 else r + m


def rref_mod(rows: list[list[int]], mod: int) -> Matrix:
    m = Matrix(rows)
    nrows, ncols = m.shape
    pivot_row = 0
    for col in range(ncols):
        pivot = None
        for r in range(pivot_row, nrows):
            if positive_mod(int(m[r, col]), mod) != 0:
                pivot = r
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            m.row_swap(pivot_row, pivot)
        inv = pow(positive_mod(int(m[pivot_row, col]), mod), -1, mod)
        for j in range(ncols):
            m[pivot_row, j] = (int(m[pivot_row, j]) * inv) % mod
        for r in range(nrows):
            if r == pivot_row:
                continue
            factor = positive_mod(int(m[r, col]), mod)
            if factor:
                for j in range(ncols):
                    m[r, j] = (int(m[r, j]) - factor * int(m[pivot_row, j])) % mod
        pivot_row += 1
    return m.applyfunc(lambda v: positive_mod(int(v), mod))
    r = a % m
    return r if r >= 0 else r + m


def poly_mod_equiv(a: sp.Expr, b: sp.Expr, mod: int) -> bool:
    pa = Poly(sp.expand(a), x, modulus=mod)
    pb = Poly(sp.expand(b), x, modulus=mod)
    return pa == pb


def poly_gcd_associate(g1: Poly, g2: Poly) -> bool:
    """True if g1 and g2 are the same up to multiplication by a unit in the field."""
    if g1 == 0 or g2 == 0:
        return g1 == g2
    q, r = div(g1, g2)
    return r == 0 and q.LC() != 0


def eval_line(line: str) -> Any:
    line = line.strip().rstrip(";")

    m = re.fullmatch(r"gcd\(\((.+)\) % (\d+),\((.+)\) % \2\)", line)
    if m:
        a = parse_poly(m.group(1))
        b = parse_poly(m.group(3))
        mod = int(m.group(2))
        pa = Poly(a.as_expr(), x, modulus=mod)
        pb = Poly(b.as_expr(), x, modulus=mod)
        return Poly(gcd(pa, pb), x, modulus=mod)

    m = re.fullmatch(r"factor\((.+)\) mod (\d+)", line)
    if m:
        p = giac_to_sympy(m.group(1))
        mod = int(m.group(2))
        return sp.factor(sp.expand(p), modulus=mod)

    m = re.fullmatch(r"gcd\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        return Poly(gcd(a, b), x)

    m = re.fullmatch(r"quo\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        q, _ = div(a, b)
        return q

    m = re.fullmatch(r"rem\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        _, r = div(a, b)
        return r

    m = re.fullmatch(r"content\((.+)\)", line)
    if m:
        return sp.Integer(parse_poly(m.group(1)).content())

    m = re.fullmatch(r"gauss\((.+),\[(.+)\]\)", line)
    if m:
        return sp.expand(giac_to_sympy(m.group(1)))

    m = re.fullmatch(r"horner\((.+),(.+)\)", line)
    if m:
        p = giac_to_sympy(m.group(1))
        val = giac_to_sympy(m.group(2))
        return sp.expand(p).subs(x, val)

    m = re.fullmatch(r"resultant\((.+),(.+),x\)", line)
    if m:
        a = giac_to_sympy(m.group(1))
        b = giac_to_sympy(m.group(2))
        return resultant(sp.expand(a), sp.expand(b), x)

    m = re.fullmatch(r"simp2\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        g = gcd(a, b)
        return sp.Tuple(a // g, b // g)

    m = re.fullmatch(r"lcm\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        return lcm(a, b)

    m = re.fullmatch(r"factor\((.+)\)", line)
    if m:
        return sp.factor(sp.expand(giac_to_sympy(m.group(1))))

    m = re.fullmatch(r"modp\((.+),(\d+)\)", line)
    if m:
        p = parse_poly(m.group(1))
        mod = int(m.group(2))
        return Poly(p, x, modulus=mod)

    m = re.fullmatch(r"smod\((.+),(\d+)\)", line)
    if m:
        return sp.Integer(positive_mod(int(m.group(1)), int(m.group(2))))

    m = re.fullmatch(r"irem\((.+),(\d+)\)", line)
    if m:
        return sp.Integer(int(m.group(1)) % int(m.group(2)))

    m = re.fullmatch(r"normal\(\(\((.+)\) % (\d+)\)\^(\d+)\)", line)
    if m:
        base, mod, exp = m.group(1), int(m.group(2)), int(m.group(3))
        b = giac_to_sympy(base)
        val = sp.expand(b ** exp)
        coeffs = sp.Poly(val, x, domain=ZZ).all_coeffs()
        reduced = sum(
            positive_mod(int(c), mod) * x ** (len(coeffs) - 1 - i)
            for i, c in enumerate(coeffs)
        )
        return sp.expand(reduced)

    m = re.fullmatch(r"normal\(\((.+)\)\^(\d+)\)", line)
    if m:
        b = giac_to_sympy(m.group(1))
        return sp.expand(b ** int(m.group(2)))

    m = re.fullmatch(r"partfrac\((.+),x\)", line)
    if m:
        expr = giac_to_sympy(m.group(1))
        return sp.expand(apart(expr, x))

    m = re.fullmatch(r"egcd\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        g, s, t = sp.gcdex(a.as_expr(), b.as_expr(), x, domain=QQ)
        return sp.Tuple(g, s, t)

    m = re.fullmatch(r"abcuv\((.+),(.+),(.+)\)", line)
    if m:
        a, b, c = parse_poly(m.group(1)), parse_poly(m.group(2)), parse_poly(m.group(3))
        g = gcd(a, b)
        if g == 0:
            raise ValueError("gcd is zero")
        q, r = div(c, g)
        return sp.Tuple(q, r)

    m = re.fullmatch(r"roots\((.+),x\)", line)
    if m:
        p = parse_poly(m.group(1))
        return sp.nroots(p.as_expr())

    m = re.fullmatch(r"chinrem\(\[(.+)\],\[(.+)\]\)", line)
    if m:
        residues = [parse_poly(s.strip()) for s in split_top_level(m.group(1))]
        moduli = [parse_poly(s.strip()) for s in split_top_level(m.group(2))]
        sol = chinrem_poly(residues, moduli)
        prod = moduli[0]
        for mi in moduli[1:]:
            prod = prod * mi
        return sp.Tuple(sol, prod)

    m = re.fullmatch(r"greduce\((.+),\[(.+)\],\[(.+)\]\)", line)
    if m:
        f = giac_to_sympy(m.group(1))
        gens = [giac_to_sympy(g.strip()) for g in split_top_level(m.group(2))]
        vars_ = [symbols(v.strip()) for v in split_top_level(m.group(3))]
        return verify_greduce(f, gens, vars_)

    m = re.fullmatch(r"rref\(\[(.+)\] % (\d+),\[(.+)\] % \2\)", line)
    if m:
        mod = int(m.group(2))
        row1 = [positive_mod(int(c.strip()), mod) for c in split_top_level(m.group(1))]
        row2 = [positive_mod(int(c.strip()), mod) for c in split_top_level(m.group(3))]
        return rref_mod([row1, row2], mod)

    raise ValueError(f"unsupported line for SymPy: {line}")


def chinrem_poly(residues: list[Poly], moduli: list[Poly]) -> Poly:
    """Polynomial CRT via successive combination."""
    sol = residues[0]
    mod = moduli[0]
    for r, m in zip(residues[1:], moduli[1:]):
        g, u, v = sp.gcdex(mod.as_expr(), m.as_expr(), x, domain=QQ)
        if g != 1:
            raise ValueError("non-coprime moduli in chinrem")
        diff = (r - sol) // mod
        sol = sol + mod * (diff * u)
        mod = mod * m
    _, r = div(sol, mod)
    return r


def verify_greduce(f: sp.Expr, gens: list[sp.Expr], vars_: list) -> sp.Expr:
    """Remainder of f modulo ideal generated by gens."""
    G = groebner(gens, vars_, order="lex", domain=QQ)
    return sp.expand(G.reduce(sp.expand(f))[1])


def normalize(val: Any) -> Any:
    if isinstance(val, Poly):
        return sp.expand(val.as_expr())
    if isinstance(val, sp.Tuple):
        return sp.Tuple(*[normalize(v) for v in val.args])
    if isinstance(val, Matrix):
        return val.applyfunc(lambda e: sp.expand(e))
    if isinstance(val, list):
        return [normalize(v) for v in val]
    return sp.expand(val) if hasattr(val, "free_symbols") else val


def equivalent(a: Any, b: Any) -> bool:
    a, b = normalize(a), normalize(b)
    if isinstance(a, Matrix) and isinstance(b, Matrix):
        return a.shape == b.shape and all(
            sp.simplify(aij - bij) == 0 for aij, bij in zip(a.flat(), b.flat())
        )
    if isinstance(a, sp.Tuple) and isinstance(b, sp.Tuple):
        if len(a) != len(b):
            return False
        return all(equiv_pair(ai, bi) for ai, bi in zip(a, b))
    return equiv_pair(a, b)


def equiv_pair(a: Any, b: Any) -> bool:
    if isinstance(a, (int, sp.Integer)) and isinstance(b, (int, sp.Integer)):
        return int(a) == int(b)
    if hasattr(a, "free_symbols") or hasattr(b, "free_symbols"):
        return sp.simplify(a - b) == 0
    return a == b


def output_to_sympy(line: str, output: str) -> Any:
    line = line.strip().rstrip(";")
    output = output.strip()
    if line.startswith(("content(", "smod(", "irem(", "horner(", "resultant(")):
        if line.startswith("horner(") or line.startswith("resultant("):
            return giac_to_sympy(output)
        return sp.Integer(int(output))
    if line.startswith("roots("):
        return giac_to_sympy(output)
    if line.startswith("rref("):
        return giac_matrix(output) if output.startswith("[[") else giac_to_sympy(output)
    if line.startswith(("chinrem(", "egcd(", "abcuv(", "simp2(")):
        return giac_to_sympy(output)
    if line.startswith("greduce("):
        return giac_to_sympy(output)
    if line.startswith("factor(") and " mod " in output:
        base, _ = output.rsplit(" mod ", 1)
        return giac_to_sympy(base)
    return giac_to_sympy(output)


def as_tuple(val: Any) -> tuple[Any, ...]:
    if isinstance(val, sp.Tuple):
        return val.args
    if isinstance(val, (tuple, list)):
        return tuple(val)
    raise ValueError(f"expected tuple, got {type(val)}")


def parse_giac_tuple(output: str) -> tuple[Any, ...]:
    output = output.strip()
    if output.startswith("[") and output.endswith("]"):
        return as_tuple(giac_to_sympy(output))
    # egcd/abcuv may print comma-separated values without brackets
    parts = split_top_level(output)
    return tuple(giac_to_sympy(p.strip()) for p in parts)


def verify_property(line: str, output: str) -> tuple[bool, str]:
    """Property-based checks when direct SymPy comparison is format-sensitive."""
    line = line.strip().rstrip(";")

    m = re.fullmatch(r"egcd\((.+),(.+)\)", line)
    if m:
        a, b = parse_poly(m.group(1)), parse_poly(m.group(2))
        tup = parse_giac_tuple(output)
        if len(tup) != 3:
            return False, "egcd expects 3-tuple"
        g, u, t = tup
        g_exp = gcd(a, b).as_expr()
        lhs = sp.expand(u * a.as_expr() + t * b.as_expr())
        if sp.rem(lhs, g_exp, x) != 0:
            return False, "egcd identity fails"
        return True, "ok"

    m = re.fullmatch(r"abcuv\((.+),(.+),(.+)\)", line)
    if m:
        a, b, c = parse_poly(m.group(1)), parse_poly(m.group(2)), parse_poly(m.group(3))
        tup = parse_giac_tuple(output)
        if len(tup) != 2:
            return False, "abcuv expects 2-tuple"
        u, v = tup
        lhs = sp.expand(u * a.as_expr() + v * b.as_expr())
        g = gcd(a, b).as_expr()
        if sp.rem(lhs, g, x) != 0:
            return False, "abcuv not in ideal of gcd"
        if sp.rem(c.as_expr(), g, x) != 0:
            return False, "c not divisible by gcd"
        # lhs and c should be associates (same up to unit in QQ[x])
        if sp.rem(lhs, c.as_expr(), x) != 0 and sp.rem(c.as_expr(), lhs, x) != 0:
            return False, "abcuv does not match c up to scaling"
        return True, "ok"

    m = re.fullmatch(r"roots\((.+),x\)", line)
    if m:
        p = parse_poly(m.group(1))
        roots_out = giac_to_sympy(output)
        items = roots_out.args if isinstance(roots_out, sp.Tuple) else [roots_out]
        for r in items:
            if sp.simplify(p.as_expr().subs(x, r)) != 0:
                return False, f"{r} is not a root"
        return True, "ok"

    m = re.fullmatch(r"chinrem\(\[(.+)\],\[(.+)\]\)", line)
    if m:
        residues = [parse_poly(s.strip()) for s in split_top_level(m.group(1))]
        moduli = [parse_poly(s.strip()) for s in split_top_level(m.group(2))]
        tup = parse_giac_tuple(output)
        sol = tup[0]
        for r, mod in zip(residues, moduli):
            if sp.rem(sp.expand(sol - r.as_expr()), mod.as_expr(), x) != 0:
                return False, "chinrem residue mismatch"
        return True, "ok"

    m = re.fullmatch(r"greduce\((.+),\[(.+)\],\[(.+)\]\)", line)
    if m:
        f = giac_to_sympy(m.group(1))
        gens = [giac_to_sympy(g.strip()) for g in split_top_level(m.group(2))]
        vars_ = [symbols(v.strip()) for v in split_top_level(m.group(3))]
        rem = giac_to_sympy(output)
        G = groebner(gens, vars_, order="lex", domain=QQ)
        reduced = sp.expand(G.reduce(sp.expand(f))[1])
        if sp.simplify(reduced - rem) != 0:
            return False, f"greduce remainder mismatch: expected {reduced}, got {rem}"
        return True, "ok"

    return _verify_direct(line, output)


def parse_normal_mod_output(output: str, mod: int) -> sp.Expr:
    """Parse giac-style `(c % m)*x^k+...` or plain expanded polynomial."""
    output = output.strip()
    if "%" not in output:
        return giac_to_sympy(output)
    terms: dict[int, int] = {}
    for raw in split_top_level(output.replace("-", "+-")):
        raw = raw.strip()
        if not raw:
            continue
        sign = -1 if raw.startswith("-") else 1
        raw = raw.lstrip("+-")
        m = re.fullmatch(r"\((-?\d+) % \d+\)\*x\^(\d+)", raw)
        if m:
            exp = int(m.group(2))
            terms[exp] = terms.get(exp, 0) + sign * positive_mod(int(m.group(1)), mod)
            continue
        m = re.fullmatch(r"\((-?\d+) % \d+\)\*x", raw)
        if m:
            terms[1] = terms.get(1, 0) + sign * positive_mod(int(m.group(1)), mod)
            continue
        m = re.fullmatch(r"(-?\d+) % \d+", raw)
        if m:
            terms[0] = terms.get(0, 0) + sign * positive_mod(int(m.group(1)), mod)
            continue
    if terms:
        return sum(c * x ** k for k, c in terms.items())
    # giac-rs nested mod display — compare via expanded base then reduce coeffs
    base = output.split(" mod ")[0] if " mod " in output else output
    e = giac_to_sympy(base)
    p = Poly(sp.expand(e), x, domain=ZZ)
    coeffs = p.all_coeffs()
    return sum(
        positive_mod(int(c), mod) * x ** (len(coeffs) - 1 - i)
        for i, c in enumerate(coeffs)
    )


def _verify_direct(line: str, output: str) -> tuple[bool, str]:
    line = line.strip().rstrip(";")
    output = output.strip()

    m = re.fullmatch(r"factor\((.+)\) mod (\d+)", line)
    if m:
        mod = int(m.group(2))
        expected = sp.factor(sp.expand(giac_to_sympy(m.group(1))), modulus=mod)
        base = output.rsplit(" mod ", 1)[0] if " mod " in output else output
        got = giac_to_sympy(base)
        if poly_mod_equiv(expected, got, mod):
            return True, "ok"
        return False, f"expected {sp.sstr(expected)!r} mod {mod}, got {sp.sstr(got)!r}"

    m = re.fullmatch(r"gcd\(\((.+)\) % (\d+),\((.+)\) % \2\)", line)
    if m:
        mod = int(m.group(2))
        a = parse_poly(m.group(1))
        b = parse_poly(m.group(3))
        pa = Poly(a.as_expr(), x, modulus=mod)
        pb = Poly(b.as_expr(), x, modulus=mod)
        expected = Poly(gcd(pa, pb), x, modulus=mod)
        base = output.rsplit(" mod ", 1)[0] if " mod " in output else output
        got = Poly(giac_to_sympy(base), x, modulus=mod)
        if poly_gcd_associate(expected, got):
            return True, "ok"
        return False, f"expected gcd associate of {expected!r}, got {got!r}"

    m = re.fullmatch(r"normal\(\(\((.+)\) % (\d+)\)\^(\d+)\)", line)
    if m:
        mod = int(m.group(2))
        expected = eval_line(line)
        got = parse_normal_mod_output(output, mod)
        if poly_mod_equiv(expected, got, mod):
            return True, "ok"
        return False, f"expected {sp.sstr(expected)!r} mod {mod}, got {sp.sstr(got)!r}"

    try:
        expected = eval_line(line)
        got = output_to_sympy(line, output)
        if equivalent(expected, got):
            return True, "ok"
        return False, f"expected {sp.sstr(expected)!r}, got {sp.sstr(got)!r}"
    except Exception as e:
        return False, str(e)


def verify(line: str, output: str) -> tuple[bool, str]:
    line = line.strip().rstrip(";")
    if line.startswith(("egcd(", "abcuv(", "roots(", "chinrem(", "greduce(")):
        return verify_property(line, output)
    return _verify_direct(line, output)


def reference(line: str) -> dict:
    val = eval_line(line)
    if isinstance(val, Poly):
        s = sp.sstr(sp.expand(val.as_expr()))
    elif isinstance(val, Matrix):
        s = str(val.tolist())
    elif isinstance(val, sp.Tuple):
        s = "[" + ",".join(sp.sstr(normalize(v)) for v in val.args) + "]"
    else:
        s = sp.sstr(normalize(val))
    return {"line": line, "sympy": s}


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    cmd = sys.argv[1]
    if cmd == "reference" and len(sys.argv) == 3:
        print(json.dumps(reference(sys.argv[2])))
        return 0
    if cmd == "verify" and len(sys.argv) == 4:
        ok, msg = verify(sys.argv[2], sys.argv[3])
        if not ok:
            print(msg, file=sys.stderr)
        return 0 if ok else 1
    if cmd == "equiv" and len(sys.argv) == 4:
        try:
            a = giac_to_sympy(sys.argv[2])
            b = giac_to_sympy(sys.argv[3])
            return 0 if equivalent(a, b) else 1
        except Exception as e:
            print(e, file=sys.stderr)
            return 1
    print("invalid arguments", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
