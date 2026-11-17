"""Square roots in ``Q_p`` and ``Z_p``."""

from __future__ import annotations

from fractions import Fraction

from padicfield.exceptions import PrecisionError, SquareRootError
from padicfield.hensel import hensel_lift_simple
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.padic_number import PadicNumber
from padicfield.utils import legendre_symbol, require_prime
from padicfield.valuation import remove_factor, valuation_int


def is_quadratic_residue_mod_p(a: int, p: int) -> bool:
    require_prime(p)
    a %= p
    if a == 0 or p == 2:
        return True
    return legendre_symbol(a, p) == 1


def tonelli_shanks(a: int, p: int) -> int:
    """Square root of ``a`` modulo odd prime ``p`` (or ``p=2``)."""
    require_prime(p)
    if p == 2:
        return a % 2
    a %= p
    if a == 0:
        return 0
    if legendre_symbol(a, p) != 1:
        raise SquareRootError(f"{a} is not a quadratic residue modulo {p}")
    Q, S = p - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    if S == 1:
        return pow(a, (p + 1) // 4, p)
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
        if z > p:
            raise SquareRootError("failed to find a quadratic non-residue")
    m, c, r, t = S, pow(z, Q, p), pow(a, (Q + 1) // 2, p), pow(a, Q, p)
    while t != 1:
        t2i, i = t, 0
        while t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1
            if i == m:
                raise SquareRootError("Tonelli–Shanks failed")
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def sqrt_mod_p(a: int, p: int) -> int:
    return tonelli_shanks(a, p)


def is_square_in_zp(a: PadicInt) -> bool:
    p = a.p
    if a.is_zero():
        return True
    v = valuation_int(a.residue, p)
    if v % 2 != 0:
        return False
    if v >= a.prec:
        return True
    unit = a.residue // (p**v)
    if p == 2:
        remaining = a.prec - v
        if remaining < 3:
            mod = 1 << remaining
            return (unit % mod) == (1 % mod)
        return unit % 8 == 1
    return is_quadratic_residue_mod_p(unit % p, p)


def is_square_in_qp(x: PadicNumber) -> bool:
    if x.is_zero():
        return True
    v = x.val
    assert v is not None
    return v % 2 == 0 and is_square_in_zp(x.unit)


def sqrt_padic_int(a: PadicInt, *, prefer_positive: bool = True) -> PadicInt:
    p, prec = a.p, a.prec
    if a.is_zero():
        return PadicInt.zero(p, prec)
    v = valuation_int(a.residue, p)
    if v % 2 != 0:
        raise SquareRootError(f"valuation {v} is odd; not a square in Z_p")
    if v >= prec:
        return PadicInt.zero(p, prec)
    half_v = v // 2
    unit, unit_prec = a.residue // (p**v), prec - v
    sqrt_unit = _sqrt_unit_2adic(unit, unit_prec) if p == 2 else _sqrt_unit_odd(unit, p, unit_prec)
    result_prec = prec - half_v
    residue = (sqrt_unit.residue * (p**half_v)) % (p**result_prec)
    result = PadicInt(residue, p, result_prec)
    if prefer_positive and p != 2 and (sqrt_unit.residue % p) > p - (sqrt_unit.residue % p):
        result = PadicInt((-result).residue, p, result_prec)
    return result


def _sqrt_unit_odd(u: int, p: int, prec: int) -> PadicInt:
    u_mod = u % p
    if not is_quadratic_residue_mod_p(u_mod, p):
        raise SquareRootError(f"unit {u_mod} is not a square mod {p}")
    r0 = tonelli_shanks(u_mod, p)
    return hensel_lift_simple([-(u % (p**prec)), 0, 1], r0, p, prec).root


def _sqrt_unit_2adic(u: int, prec: int) -> PadicInt:
    if prec < 1:
        raise PrecisionError("precision must be >= 1")
    if prec >= 3 and u % 8 != 1:
        raise SquareRootError(f"2-adic unit {u % 8} ≢ 1 (mod 8); not a square")
    if prec < 3:
        if u % 2 == 0:
            raise SquareRootError("not a unit")
        return PadicInt(1, 2, prec)
    x, current = 1, 3
    while current < prec:
        next_prec = min(current + 1, prec)
        mod = 1 << next_prec
        found = None
        for b in (0, 1):
            cand = (x + b * (1 << current)) % mod
            if (cand * cand - u) % mod == 0:
                found = cand
                break
        if found is None and current >= 3:
            half = 1 << (current - 1)
            span = 1 << (next_prec - current + 1)
            for t in range(span):
                cand = (x + t * half) % mod
                if (cand * cand - u) % mod == 0:
                    found = cand
                    break
        if found is None:
            raise SquareRootError(f"failed to lift 2-adic square root at precision {current}")
        x, current = found, next_prec
    return PadicInt(x, 2, prec)


def sqrt_padic_number(x: PadicNumber) -> PadicNumber:
    if x.is_zero():
        return PadicNumber.zero(x.p, x.prec)
    v = x.val
    assert v is not None
    if v % 2 != 0:
        raise SquareRootError(f"valuation {v} is odd; not a square in Q_p")
    return PadicNumber(sqrt_padic_int(x.unit), val=v // 2)


def sqrt_rational(x: int | Fraction, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
    return sqrt_padic_number(PadicNumber.from_rational(x, p, prec))


def all_square_roots_int(a: PadicInt) -> list[PadicInt]:
    r = sqrt_padic_int(a)
    roots = [r]
    neg = -r
    if neg != r:
        roots.append(neg)
    return roots


def square_root_exists_rational(x: int | Fraction, p: int) -> bool:
    require_prime(p)
    frac = Fraction(x)
    if frac == 0:
        return True
    v_num, u_num = remove_factor(frac.numerator, p)
    v_den, u_den = remove_factor(frac.denominator, p)
    if (v_num - v_den) % 2 != 0:
        return False
    if frac < 0 and not _minus_one_is_square(p):
        return False
    unit_prod = abs(u_num * u_den)
    if p == 2:
        return unit_prod % 8 == 1
    return is_quadratic_residue_mod_p(unit_prod % p, p)


def _minus_one_is_square(p: int) -> bool:
    """``-1`` is a square in ``Q_p`` iff ``p = 2`` is false and ``p ≡ 1 (mod 4)``."""
    if p == 2:
        return False
    return legendre_symbol(-1, p) == 1
