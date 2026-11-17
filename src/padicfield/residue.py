"""Quadratic residue symbols modulo ``p^k`` and the Hilbert symbol."""

from __future__ import annotations

from fractions import Fraction

from padicfield.exceptions import ResidueSymbolError
from padicfield.padic_number import PadicNumber
from padicfield.sqrt import is_quadratic_residue_mod_p
from padicfield.utils import jacobi_symbol, kronecker_symbol, legendre_symbol, require_prime
from padicfield.valuation import remove_factor, valuation_int


def quadratic_residue_symbol(a: int, p: int) -> int:
    """Legendre ``(a/p)`` for odd ``p``; Kronecker for ``p=2``."""
    require_prime(p)
    return kronecker_symbol(a, 2) if p == 2 else legendre_symbol(a, p)


def _unit_square_symbol(u: int, p: int, prec: int) -> int:
    if u % p == 0:
        raise ResidueSymbolError("expected a unit")
    if p == 2:
        if prec == 1:
            return 1
        if prec == 2:
            return 1 if u % 4 == 1 else -1
        return 1 if u % 8 == 1 else -1
    return 1 if is_quadratic_residue_mod_p(u % p, p) else -1


def residue_symbol_mod_pk(a: int, p: int, k: int) -> int:
    """Square symbol of ``a`` modulo ``p^k`` in ``{-1, 0, 1}``."""
    require_prime(p)
    if k < 1:
        raise ResidueSymbolError("k must be >= 1")
    a %= p**k
    if a == 0:
        return 0
    v = valuation_int(a, p)
    if v >= k:
        return 0
    if v % 2 == 1:
        return -1
    return _unit_square_symbol(a // (p**v), p, k - v)


def is_square_mod_pk(a: int, p: int, k: int) -> bool:
    require_prime(p)
    a %= p**k
    if a == 0:
        return True
    v = valuation_int(a, p)
    if v % 2 == 1:
        return False
    if v >= k:
        return True
    return _unit_square_symbol((a // (p**v)) % (p ** (k - v)), p, k - v) == 1


def _eps(u: int) -> int:
    return 1 if ((u - 1) // 2) % 2 == 0 else -1


def _omega(u: int) -> int:
    return 1 if ((u * u - 1) // 8) % 2 == 0 else -1


def _val_unit_fraction(x: Fraction, p: int) -> tuple[int, int]:
    vn, un = remove_factor(x.numerator, p)
    vd, ud = remove_factor(x.denominator, p)
    return vn - vd, un * ud


def _hilbert_at_odd(va: int, ua: int, vb: int, ub: int, p: int) -> int:
    sign = 1
    if (va * vb) % 2 == 1 and _eps(p) == -1:
        sign = -1
    la, lb = legendre_symbol(ua % p, p), legendre_symbol(ub % p, p)
    if vb % 2 == 1:
        sign *= la
    if va % 2 == 1:
        sign *= lb
    return sign


def _hilbert_at_2(va: int, ua: int, vb: int, ub: int) -> int:
    if ua % 2 == 0 or ub % 2 == 0:
        raise ResidueSymbolError("unit parts must be odd at p=2")
    e = 0
    if _eps(ua) == -1 and _eps(ub) == -1:
        e += 1
    if va % 2 == 1 and _omega(ub) == -1:
        e += 1
    if vb % 2 == 1 and _omega(ua) == -1:
        e += 1
    return -1 if e % 2 else 1


def hilbert_symbol(a: int | Fraction, b: int | Fraction, p: int) -> int:
    """Hilbert symbol ``(a, b)_p ∈ {±1}``."""
    require_prime(p)
    a, b = Fraction(a), Fraction(b)
    if a == 0 or b == 0:
        raise ResidueSymbolError("Hilbert symbol requires nonzero a, b")
    va, ua = _val_unit_fraction(a, p)
    vb, ub = _val_unit_fraction(b, p)
    return _hilbert_at_2(va, ua, vb, ub) if p == 2 else _hilbert_at_odd(va, ua, vb, ub, p)


def hilbert_symbol_padic(a: PadicNumber, b: PadicNumber) -> int:
    if a.p != b.p:
        raise ResidueSymbolError("primes must match")
    if a.is_zero() or b.is_zero():
        raise ResidueSymbolError("Hilbert symbol requires nonzero arguments")
    p, va, vb = a.p, a.val or 0, b.val or 0
    ua, ub = a.unit.residue, b.unit.residue
    if p == 2:
        if a.unit.prec < 3 or b.unit.prec < 3:
            raise ResidueSymbolError("need precision >= 3 for Hilbert symbol at 2")
        return _hilbert_at_2(va, ua, vb, ub)
    return _hilbert_at_odd(va, ua % p, vb, ub % p, p)


def quadratic_defect(a: "PadicInt") -> int | float:  # noqa: F821
    from padicfield.padic_int import PadicInt
    from padicfield.sqrt import is_square_in_zp, sqrt_padic_int

    assert isinstance(a, PadicInt)
    p = a.p
    if a.is_zero():
        return float("inf")
    if is_square_in_zp(a):
        try:
            diff = a - sqrt_padic_int(a) * sqrt_padic_int(a)
            if diff.is_zero():
                return float(a.prec)
            return float(valuation_int(diff.residue, p))
        except Exception:
            return float(a.prec)
    return float(valuation_int(a.residue, p))


def enumerate_squares_mod_pk(p: int, k: int) -> set[int]:
    require_prime(p)
    if k < 1:
        raise ResidueSymbolError("k must be >= 1")
    mod = p**k
    if mod > 100_000:
        raise ResidueSymbolError("modulus too large to enumerate")
    return {(x * x) % mod for x in range(mod)}


def square_classes_qp(p: int) -> list[str]:
    require_prime(p)
    if p == 2:
        return ["1", "3", "5", "7", "2", "6", "10", "14"]
    u = 2
    while is_quadratic_residue_mod_p(u, p):
        u += 1
    return ["1", str(u), str(p), str(u * p)]


__all__ = [
    "quadratic_residue_symbol",
    "residue_symbol_mod_pk",
    "is_square_mod_pk",
    "hilbert_symbol",
    "hilbert_symbol_padic",
    "quadratic_defect",
    "enumerate_squares_mod_pk",
    "square_classes_qp",
    "legendre_symbol",
    "jacobi_symbol",
    "kronecker_symbol",
]
