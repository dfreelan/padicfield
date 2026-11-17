"""Polynomials over ``Z`` and ``Z_p`` with Hensel-friendly helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Sequence

from padicfield.exceptions import PolynomialError, PrecisionError
from padicfield.hensel import (
    eval_poly,
    find_roots_mod_p,
    hensel_lift_simple,
    lift_all_simple_roots,
    poly_derivative_coeffs,
)
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.utils import gcd, require_prime
from padicfield.valuation import gauss_valuation


def _normalize(coeffs: Sequence[int]) -> list[int]:
    out = [int(c) for c in coeffs]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [0]


@dataclass
class IntPolynomial:
    """Polynomial with integer coefficients (low degree first)."""

    coeffs: list[int] = field(default_factory=lambda: [0])

    def __post_init__(self) -> None:
        self.coeffs = _normalize(self.coeffs)

    @classmethod
    def from_coeffs(cls, *coeffs: int) -> IntPolynomial:
        return cls(list(coeffs))

    @classmethod
    def monomial(cls, degree: int, coeff: int = 1) -> IntPolynomial:
        if degree < 0:
            raise PolynomialError("degree must be non-negative")
        return cls([0] * degree + [coeff])

    @classmethod
    def zero(cls) -> IntPolynomial:
        return cls([0])

    @classmethod
    def one(cls) -> IntPolynomial:
        return cls([1])

    @property
    def degree(self) -> int:
        return -1 if self.is_zero() else len(self.coeffs) - 1

    def is_zero(self) -> bool:
        return all(c == 0 for c in self.coeffs)

    def content(self) -> int:
        g = 0
        for c in self.coeffs:
            g = gcd(g, abs(c))
        return g

    def primitive_part(self) -> IntPolynomial:
        g = self.content()
        return IntPolynomial.zero() if g == 0 else IntPolynomial([c // g for c in self.coeffs])

    def __add__(self, other: IntPolynomial | int) -> IntPolynomial:
        if isinstance(other, int):
            other = IntPolynomial([other])
        n = max(len(self.coeffs), len(other.coeffs))
        a = self.coeffs + [0] * (n - len(self.coeffs))
        b = other.coeffs + [0] * (n - len(other.coeffs))
        return IntPolynomial([x + y for x, y in zip(a, b)])

    def __radd__(self, other: int) -> IntPolynomial:
        return self + other

    def __sub__(self, other: IntPolynomial | int) -> IntPolynomial:
        return self + (-(other if isinstance(other, IntPolynomial) else IntPolynomial([other])))

    def __rsub__(self, other: int) -> IntPolynomial:
        return IntPolynomial([other]) - self

    def __neg__(self) -> IntPolynomial:
        return IntPolynomial([-c for c in self.coeffs])

    def __mul__(self, other: IntPolynomial | int) -> IntPolynomial:
        if isinstance(other, int):
            return IntPolynomial([c * other for c in self.coeffs])
        if self.is_zero() or other.is_zero():
            return IntPolynomial.zero()
        out = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                out[i + j] += a * b
        return IntPolynomial(out)

    def __rmul__(self, other: int) -> IntPolynomial:
        return self * other

    def __pow__(self, n: int) -> IntPolynomial:
        if n < 0:
            raise PolynomialError("negative powers not supported")
        result, base = IntPolynomial.one(), self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __eq__(self, other: object) -> bool:
        if isinstance(other, int):
            return self.coeffs == [other]
        if not isinstance(other, IntPolynomial):
            return NotImplemented
        return self.coeffs == other.coeffs

    def __hash__(self) -> int:
        return hash(tuple(self.coeffs))

    def evaluate(self, x: int, mod: int | None = None) -> int:
        return eval_poly(self.coeffs, x, mod)

    def derivative(self) -> IntPolynomial:
        return IntPolynomial(poly_derivative_coeffs(self.coeffs))

    def gauss_valuation(self, p: int) -> int:
        return gauss_valuation(self.coeffs, p)

    def reduce_mod(self, m: int) -> IntPolynomial:
        if m <= 0:
            raise ValueError("modulus must be positive")
        return IntPolynomial([c % m for c in self.coeffs])

    def roots_mod_p(self, p: int) -> list[int]:
        return find_roots_mod_p(self.coeffs, p)

    def hensel_lift_root(self, a0: int, p: int, prec: int) -> PadicInt:
        return hensel_lift_simple(self.coeffs, a0, p, prec).root

    def lift_all_simple_roots(self, p: int, prec: int) -> list[PadicInt]:
        return [r.root for r in lift_all_simple_roots(self.coeffs, p, prec)]

    def __repr__(self) -> str:
        return f"IntPolynomial({self.coeffs})"

    def __str__(self) -> str:
        if self.is_zero():
            return "0"
        parts = []
        for i, c in enumerate(self.coeffs):
            if c == 0:
                continue
            if i == 0:
                parts.append(str(c))
            elif i == 1:
                parts.append("x" if c == 1 else ("-x" if c == -1 else f"{c}*x"))
            else:
                parts.append(f"x^{i}" if c == 1 else (f"-x^{i}" if c == -1 else f"{c}*x^{i}"))
        return " + ".join(parts).replace("+ -", "- ")


@dataclass
class PadicPolynomial:
    """Polynomial with uniform-precision ``PadicInt`` coefficients."""

    coeffs: list[PadicInt]
    p: int
    prec: int

    def __post_init__(self) -> None:
        require_prime(self.p)
        if self.prec < 1:
            raise PrecisionError("precision must be >= 1")
        normalized = []
        for c in self.coeffs:
            if c.p != self.p:
                raise PolynomialError(f"coefficient prime {c.p} != {self.p}")
            normalized.append(c.with_precision(self.prec))
        while len(normalized) > 1 and normalized[-1].is_zero():
            normalized.pop()
        self.coeffs = normalized or [PadicInt.zero(self.p, self.prec)]

    @classmethod
    def from_int_poly(cls, poly: IntPolynomial, p: int, prec: int = DEFAULT_PRECISION) -> PadicPolynomial:
        return cls([PadicInt(c, p, prec) for c in poly.coeffs], p, prec)

    @classmethod
    def from_int_coeffs(cls, coeffs: Sequence[int], p: int, prec: int = DEFAULT_PRECISION) -> PadicPolynomial:
        return cls([PadicInt(c, p, prec) for c in coeffs], p, prec)

    @classmethod
    def zero(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicPolynomial:
        return cls([PadicInt.zero(p, prec)], p, prec)

    @classmethod
    def one(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicPolynomial:
        return cls([PadicInt.one(p, prec)], p, prec)

    @property
    def degree(self) -> int:
        return -1 if self.is_zero() else len(self.coeffs) - 1

    def is_zero(self) -> bool:
        return all(c.is_zero() for c in self.coeffs)

    def evaluate(self, x: PadicInt | int) -> PadicInt:
        if isinstance(x, int):
            x = PadicInt(x, self.p, self.prec)
        if x.p != self.p:
            raise PolynomialError("prime mismatch in evaluation")
        acc = PadicInt.zero(self.p, min(self.prec, x.prec))
        for c in reversed(self.coeffs):
            acc = acc * x + c
        return acc

    def derivative(self) -> PadicPolynomial:
        if len(self.coeffs) <= 1:
            return PadicPolynomial.zero(self.p, self.prec)
        return PadicPolynomial([self.coeffs[i] * i for i in range(1, len(self.coeffs))], self.p, self.prec)

    def __add__(self, other: PadicPolynomial) -> PadicPolynomial:
        if self.p != other.p:
            raise PolynomialError(f"prime mismatch: {self.p} vs {other.p}")
        n, prec = max(len(self.coeffs), len(other.coeffs)), min(self.prec, other.prec)
        out = []
        for i in range(n):
            a = self.coeffs[i] if i < len(self.coeffs) else PadicInt.zero(self.p, prec)
            b = other.coeffs[i] if i < len(other.coeffs) else PadicInt.zero(self.p, prec)
            out.append(a.with_precision(prec) + b.with_precision(prec))
        return PadicPolynomial(out, self.p, prec)

    def __mul__(self, other: PadicPolynomial | PadicInt | int) -> PadicPolynomial:
        if isinstance(other, int):
            other = PadicInt(other, self.p, self.prec)
        if isinstance(other, PadicInt):
            return PadicPolynomial(
                [c * other for c in self.coeffs], self.p, min(self.prec, other.prec)
            )
        if self.p != other.p:
            raise PolynomialError(f"prime mismatch: {self.p} vs {other.p}")
        prec = min(self.prec, other.prec)
        out = [PadicInt.zero(self.p, prec) for _ in range(len(self.coeffs) + len(other.coeffs) - 1)]
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                out[i + j] = out[i + j] + a.with_precision(prec) * b.with_precision(prec)
        return PadicPolynomial(out, self.p, prec)

    def __repr__(self) -> str:
        return f"PadicPolynomial({[c.residue for c in self.coeffs]}, p={self.p}, prec={self.prec})"

    def to_int_poly_approx(self) -> IntPolynomial:
        return IntPolynomial([c.residue for c in self.coeffs])


def compose(f: IntPolynomial, g: IntPolynomial) -> IntPolynomial:
    """Return ``f ∘ g``."""
    result, power = IntPolynomial.zero(), IntPolynomial.one()
    for c in f.coeffs:
        result = result + c * power
        power = power * g
    return result


def resultant_univariate(f: IntPolynomial, g: IntPolynomial) -> int:
    """``Res(f, g)`` via Sylvester determinant."""
    if f.is_zero() or g.is_zero():
        return 0
    if f.degree < g.degree:
        sign = -1 if (f.degree * g.degree) % 2 else 1
        return sign * resultant_univariate(g, f)
    n, m = f.degree, g.degree
    if m == 0:
        return g.coeffs[0] ** n
    size = n + m
    S = [[0] * size for _ in range(size)]
    for i in range(m):
        for j, c in enumerate(f.coeffs):
            S[i][i + j] = c
    for i in range(n):
        for j, c in enumerate(g.coeffs):
            S[m + i][i + j] = c
    F = [[Fraction(x) for x in row] for row in S]
    det = Fraction(1)
    for k in range(size):
        pivot = next((i for i in range(k, size) if F[i][k] != 0), None)
        if pivot is None:
            return 0
        if pivot != k:
            F[k], F[pivot] = F[pivot], F[k]
            det = -det
        det *= F[k][k]
        pk = F[k][k]
        for i in range(k + 1, size):
            if F[i][k] == 0:
                continue
            factor = F[i][k] / pk
            for j in range(k, size):
                F[i][j] -= factor * F[k][j]
    return int(det)


def discriminant(f: IntPolynomial) -> int:
    """Discriminant via ``Res(f, f') / lead`` with the usual sign."""
    if f.degree < 1:
        raise PolynomialError("discriminant requires degree >= 1")
    n, a_n = f.degree, f.coeffs[-1]
    res = resultant_univariate(f, f.derivative())
    sign = -1 if (n * (n - 1) // 2) % 2 else 1
    return sign * (res // a_n)
