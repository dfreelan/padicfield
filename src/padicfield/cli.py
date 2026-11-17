"""Command-line interface for padicfield.

Usage examples::

    padicfield val 360 5
    padicfield expand 42 5 --prec 8
    padicfield teich 3 7 --prec 6
    padicfield sqrt 2 7 --prec 8
    padicfield hensel 1 0 -2 7 --root 3 --prec 6
    padicfield legendre 10 13
    padicfield hilbert 2 5 7
    padicfield info 5
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from typing import Sequence


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="padicfield",
        description="p-adic arithmetic and local-field numerics",
    )
    parser.add_argument("--version", action="version", version="padicfield 0.1.0")
    sub = parser.add_subparsers(dest="command", required=True)

    # val
    p_val = sub.add_parser("val", help="compute v_p(n) for an integer or rational")
    p_val.add_argument("n", type=str, help="integer or rational (a/b)")
    p_val.add_argument("p", type=int, help="prime")

    # norm
    p_norm = sub.add_parser("norm", help="compute |n|_p")
    p_norm.add_argument("n", type=str)
    p_norm.add_argument("p", type=int)

    # expand
    p_exp = sub.add_parser("expand", help="base-p expansion of an integer in Z_p")
    p_exp.add_argument("n", type=int)
    p_exp.add_argument("p", type=int)
    p_exp.add_argument("--prec", type=int, default=10)

    # teich
    p_teich = sub.add_parser("teich", help="Teichmüller representative ω(a) mod p^prec")
    p_teich.add_argument("a", type=int)
    p_teich.add_argument("p", type=int)
    p_teich.add_argument("--prec", type=int, default=8)

    # sqrt
    p_sqrt = sub.add_parser("sqrt", help="p-adic square root of a rational")
    p_sqrt.add_argument("n", type=str)
    p_sqrt.add_argument("p", type=int)
    p_sqrt.add_argument("--prec", type=int, default=10)

    # hensel
    p_hensel = sub.add_parser(
        "hensel",
        help="lift a simple root of a polynomial (coeffs c0 c1 ... cd)",
    )
    p_hensel.add_argument("coeffs", nargs="+", type=int, help="coefficients c0 c1 ...")
    p_hensel.add_argument("--p", type=int, required=True)
    p_hensel.add_argument("--root", type=int, required=True, help="root modulo p")
    p_hensel.add_argument("--prec", type=int, default=8)

    # legendre
    p_leg = sub.add_parser("legendre", help="Legendre/Kronecker symbol (a/p)")
    p_leg.add_argument("a", type=int)
    p_leg.add_argument("p", type=int)

    # hilbert
    p_hilb = sub.add_parser("hilbert", help="Hilbert symbol (a,b)_p")
    p_hilb.add_argument("a", type=str)
    p_hilb.add_argument("b", type=str)
    p_hilb.add_argument("p", type=int)

    # info
    p_info = sub.add_parser("info", help="print structural info about Q_p")
    p_info.add_argument("p", type=int)
    p_info.add_argument("--prec", type=int, default=8)

    # units
    p_units = sub.add_parser("units", help="describe (Z/p^kZ)^×")
    p_units.add_argument("p", type=int)
    p_units.add_argument("k", type=int)

    # add / mul demo
    p_arith = sub.add_parser("arith", help="add or multiply two integers in Z_p")
    p_arith.add_argument("op", choices=["add", "mul", "div"])
    p_arith.add_argument("a", type=int)
    p_arith.add_argument("b", type=int)
    p_arith.add_argument("p", type=int)
    p_arith.add_argument("--prec", type=int, default=10)

    return parser


def _parse_rational(s: str) -> Fraction:
    return Fraction(s)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``padicfield`` console script."""
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    try:
        return _dispatch(args)
    except Exception as exc:  # noqa: BLE001 — CLI surfaces all errors cleanly
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace) -> int:
    cmd = args.command

    if cmd == "val":
        from padicfield.valuation import valuation

        x = _parse_rational(args.n)
        print(valuation(x, args.p))
        return 0

    if cmd == "norm":
        from padicfield.valuation import padic_norm

        x = _parse_rational(args.n)
        print(padic_norm(x, args.p))
        return 0

    if cmd == "expand":
        from padicfield.padic_int import PadicInt

        x = PadicInt(args.n, args.p, args.prec)
        digs = x.digits()
        print(f"digits (LSB first): {digs}")
        print(x)
        return 0

    if cmd == "teich":
        from padicfield.teichmuller import teichmuller_representative

        w = teichmuller_representative(args.a, args.p, args.prec)
        print(f"ω({args.a}) ≡ {w.residue} (mod {args.p}^{args.prec})")
        print(w)
        return 0

    if cmd == "sqrt":
        from padicfield.sqrt import sqrt_rational

        r = sqrt_rational(_parse_rational(args.n), args.p, args.prec)
        print(r)
        return 0

    if cmd == "hensel":
        from padicfield.hensel import hensel_lift_simple

        result = hensel_lift_simple(args.coeffs, args.root, args.p, args.prec)
        print(f"lifted root ≡ {result.root.residue} (mod {args.p}^{args.prec})")
        return 0

    if cmd == "legendre":
        from padicfield.residue import quadratic_residue_symbol

        print(quadratic_residue_symbol(args.a, args.p))
        return 0

    if cmd == "hilbert":
        from padicfield.residue import hilbert_symbol

        print(hilbert_symbol(_parse_rational(args.a), _parse_rational(args.b), args.p))
        return 0

    if cmd == "info":
        from padicfield.local_field import LocalField, describe_unit_group
        from padicfield.residue import square_classes_qp

        K = LocalField(args.p, args.prec)
        print(K)
        print(f"unit group mod p^2: {describe_unit_group(args.p, 2)}")
        print(f"square classes Q_p^×/(Q_p^×)^2: {square_classes_qp(args.p)}")
        print(f"exp converges for v_p(x) > {1 / (args.p - 1)}")
        return 0

    if cmd == "units":
        from padicfield.local_field import describe_unit_group
        from padicfield.zeta import count_units_mod_pk

        print(describe_unit_group(args.p, args.k))
        print(f"order = {count_units_mod_pk(args.p, args.k)}")
        return 0

    if cmd == "arith":
        from padicfield.padic_int import PadicInt

        a = PadicInt(args.a, args.p, args.prec)
        b = PadicInt(args.b, args.p, args.prec)
        if args.op == "add":
            print(a + b)
        elif args.op == "mul":
            print(a * b)
        else:
            print(a / b)
        return 0

    print(f"unknown command {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
