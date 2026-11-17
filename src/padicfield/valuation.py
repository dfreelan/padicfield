"""p-adic valuations and absolute values."""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Union

from padicfield.exceptions import ValuationError
from padicfield.utils import factorize, require_prime

NumberLike = Union[int, Fraction]


def valuation_int(n: int, p: int) -> int:
    """Return finite ``v_p(n)`` for nonzero ``n``."""
    require_prime(p)
    if n == 0:
        raise ValuationError("valuation of zero is undefined as a finite integer")
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def valuation(x: NumberLike, p: int) -> float:
    """Extended valuation ``v_p(x) ∈ Z ∪ {+∞}``; zero maps to ``math.inf``."""
    require_prime(p)
    if isinstance(x, Fraction):
        if x == 0:
            return math.inf
        return float(valuation_int(x.numerator, p) - valuation_int(x.denominator, p))
    if isinstance(x, int):
        if x == 0:
            return math.inf
        return float(valuation_int(x, p))
    raise TypeError(f"unsupported type for valuation: {type(x).__name__}")


def unit_part_int(n: int, p: int) -> int:
    """Signed p-free part of nonzero ``n``."""
    require_prime(p)
    if n == 0:
        raise ValuationError("unit part of zero is undefined")
    sign = -1 if n < 0 else 1
    n = abs(n)
    while n % p == 0:
        n //= p
    return sign * n


def remove_factor(n: int, p: int) -> tuple[int, int]:
    """Return ``(v_p(n), n / p^{v_p(n)})`` for nonzero ``n``."""
    require_prime(p)
    if n == 0:
        raise ValuationError("cannot remove factor from zero")
    sign = -1 if n < 0 else 1
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v, sign * n


def padic_norm(x: NumberLike, p: int) -> float:
    """``|x|_p = p^{-v_p(x)}``, with ``|0|_p = 0``."""
    v = valuation(x, p)
    if math.isinf(v):
        return 0.0
    return float(p) ** (-v)


def padic_distance(x: NumberLike, y: NumberLike, p: int) -> float:
    """``|x - y|_p``."""
    diff: NumberLike = Fraction(x) - Fraction(y)
    return padic_norm(diff, p)


def is_p_adic_integer_rational(x: NumberLike, p: int) -> bool:
    """True iff ``v_p(x) >= 0`` (including zero)."""
    v = valuation(x, p)
    return math.isinf(v) or v >= 0


def normalize_rational(x: NumberLike) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    raise TypeError(f"expected int or Fraction, got {type(x).__name__}")


def split_rational(x: NumberLike, p: int) -> tuple[int, Fraction]:
    """Write ``x = p^v * u`` with ``v_p(u) = 0``."""
    require_prime(p)
    frac = normalize_rational(x)
    if frac == 0:
        raise ValuationError("cannot split zero")
    v = valuation_int(frac.numerator, p) - valuation_int(frac.denominator, p)
    u = Fraction(unit_part_int(frac.numerator, p), unit_part_int(frac.denominator, p))
    return v, u


def gauss_valuation(coeffs: list[int], p: int) -> int:
    """``min v_p(c_i)`` over nonzero coefficients."""
    require_prime(p)
    vals = [valuation_int(c, p) for c in coeffs if c != 0]
    if not vals:
        raise ValuationError("Gauss valuation of the zero polynomial is +∞")
    return min(vals)


def valuation_of_factorial(n: int, p: int) -> int:
    """``v_p(n!)`` by de Polignac's formula."""
    require_prime(p)
    if n < 0:
        raise ValueError("n must be non-negative")
    v, pp = 0, p
    while pp <= n:
        v += n // pp
        if pp > n // p:
            break
        pp *= p
    return v


def valuation_of_binomial(n: int, k: int, p: int) -> int:
    """``v_p(C(n, k))`` via Legendre."""
    require_prime(p)
    if k < 0 or k > n:
        raise ValuationError("binomial coefficient is zero; valuation is +∞")
    return (
        valuation_of_factorial(n, p)
        - valuation_of_factorial(k, p)
        - valuation_of_factorial(n - k, p)
    )


def kummer_carry_count(n: int, k: int, p: int) -> int:
    """Kummer: number of carries adding ``k`` and ``n-k`` in base ``p``."""
    require_prime(p)
    if k < 0 or k > n:
        raise ValuationError("binomial coefficient is zero; valuation is +∞")
    a, b, carries, carry = k, n - k, 0, 0
    while a > 0 or b > 0:
        digit_sum = (a % p) + (b % p) + carry
        carry = 1 if digit_sum >= p else 0
        carries += carry
        a //= p
        b //= p
    return carries


def radical_prime_support(n: int) -> set[int]:
    if n == 0:
        raise ValueError("radical of zero is undefined")
    return set(factorize(abs(n)).keys())


def compare_valuations(v1: float, v2: float) -> int:
    if math.isinf(v1) and math.isinf(v2):
        return 0
    if math.isinf(v1):
        return 1
    if math.isinf(v2):
        return -1
    return (v1 > v2) - (v1 < v2)


def min_valuation(*vs: float) -> float:
    if not vs:
        raise ValueError("min_valuation requires at least one argument")
    best = vs[0]
    for v in vs[1:]:
        if compare_valuations(v, best) < 0:
            best = v
    return best


class ValuationRing:
    """Predicate helper for the localization ``Z_{(p)}``."""

    __slots__ = ("p",)

    def __init__(self, p: int) -> None:
        self.p = require_prime(p)

    def contains(self, x: NumberLike) -> bool:
        return is_p_adic_integer_rational(x, self.p)

    def reduce_mod_p(self, x: NumberLike) -> int:
        if not self.contains(x):
            raise ValuationError(f"{x} is not in Z_({self.p})")
        from padicfield.utils import modinv

        frac = normalize_rational(x)
        a = frac.numerator % self.p
        b = frac.denominator % self.p
        return (a * modinv(b, self.p)) % self.p

    def __repr__(self) -> str:
        return f"ValuationRing(p={self.p})"
