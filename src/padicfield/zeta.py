"""Simple local zeta-function and point-counting helpers."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Sequence

from padicfield.exceptions import DomainError
from padicfield.hensel import eval_poly, find_roots_mod_p, multiplicity_mod_p
from padicfield.polynomials import IntPolynomial
from padicfield.utils import require_prime


@dataclass(frozen=True)
class LocalCount:
    """Solution count for a congruence modulo ``p^k``."""

    p: int
    k: int
    count: int
    total: int

    @property
    def density(self) -> Fraction:
        return Fraction(self.count, self.total)


def count_zeros_univariate_mod_pk(coeffs: Sequence[int], p: int, k: int) -> LocalCount:
    """Exhaustive root count modulo ``p^k`` (small moduli only)."""
    require_prime(p)
    if k < 1:
        raise DomainError("k must be >= 1")
    mod = p**k
    if mod > 50_000:
        raise DomainError(f"modulus p^k = {mod} too large for exhaustive count")
    count = sum(1 for x in range(mod) if eval_poly(coeffs, x, mod) % mod == 0)
    return LocalCount(p=p, k=k, count=count, total=mod)


def multiplicative_group_point_count(p: int, k: int = 1) -> int:
    """``|(Z/p^k Z)^×| = φ(p^k)``."""
    require_prime(p)
    if k < 1:
        raise DomainError("k must be >= 1")
    return p - 1 if k == 1 else p**k - p ** (k - 1)


def count_units_mod_pk(p: int, k: int) -> int:
    return multiplicative_group_point_count(p, k)


def projective_line_point_count(p: int) -> int:
    require_prime(p)
    return p + 1


def affine_line_point_count(p: int, k: int = 1) -> int:
    require_prime(p)
    return p**k


def euler_factor_affine_line(p: int) -> Callable[[object], Fraction]:
    """Local factor ``1/(1 - p T)`` of the affine line."""
    require_prime(p)

    def factor(T: object) -> Fraction:
        return Fraction(1) / (1 - p * Fraction(T))  # type: ignore[arg-type]

    return factor


def euler_factor_projective_line(p: int) -> Callable[[object], Fraction]:
    """Local factor ``1/((1-T)(1-p T))`` of ``P^1``."""
    require_prime(p)

    def factor(T: object) -> Fraction:
        t = Fraction(T)  # type: ignore[arg-type]
        return Fraction(1) / ((1 - t) * (1 - p * t))

    return factor


def igusa_naive_poincare(coeffs: Sequence[int], p: int, max_k: int) -> list[Fraction]:
    """Densities ``N_k / p^k`` for zeros of a univariate polynomial."""
    require_prime(p)
    if max_k < 1:
        raise DomainError("max_k must be >= 1")
    return [count_zeros_univariate_mod_pk(coeffs, p, k).density for k in range(1, max_k + 1)]


def local_zeta_motive_guess(poly: IntPolynomial, p: int) -> str:
    """Heuristic reduction-type summary from roots in ``F_p``."""
    require_prime(p)
    if poly.is_zero():
        return "identically zero"
    roots = poly.roots_mod_p(p)
    if not roots:
        return "no roots in F_p (smooth reduction away from infinity)"
    mults = {r: multiplicity_mod_p(poly.coeffs, r, p) for r in roots}
    if all(m == 1 for m in mults.values()):
        return f"square-free with {len(roots)} simple root(s) in F_p: {sorted(mults)}"
    return f"roots with multiplicities in F_p: {dict(sorted(mults.items()))}"


def density_of_units(p: int) -> Fraction:
    require_prime(p)
    return Fraction(p - 1, p)


def density_of_squares_in_units(p: int) -> Fraction:
    require_prime(p)
    return Fraction(1, 8) if p == 2 else Fraction(p - 1, 2 * p)


def count_nonzero_mod_p(coeffs: Sequence[int], p: int) -> int:
    return p - len(find_roots_mod_p(coeffs, p))
