"""Local-field context objects for ``Q_p`` and finite extensions."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from padicfield.exceptions import DomainError, PrecisionError
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.padic_number import PadicNumber
from padicfield.utils import require_prime


@dataclass(frozen=True)
class ExtensionData:
    """Invariants of a finite extension ``K / Q_p`` with ``n = e f``."""

    p: int
    n: int
    e: int
    f: int

    def __post_init__(self) -> None:
        require_prime(self.p)
        if self.n < 1 or self.e < 1 or self.f < 1:
            raise DomainError("n, e, f must be positive")
        if self.e * self.f != self.n:
            raise DomainError(f"expected e*f = n, got {self.e}*{self.f} != {self.n}")

    @property
    def is_unramified(self) -> bool:
        return self.e == 1

    @property
    def is_totally_ramified(self) -> bool:
        return self.f == 1

    @property
    def is_trivial(self) -> bool:
        return self.n == 1

    def residue_cardinality(self) -> int:
        return self.p**self.f


class LocalField:
    """Working context for ``Q_p`` at a default precision."""

    __slots__ = ("_p", "_prec", "_extension")

    def __init__(
        self,
        p: int,
        prec: int = DEFAULT_PRECISION,
        *,
        extension: ExtensionData | None = None,
    ) -> None:
        self._p = require_prime(p)
        if prec < 1:
            raise PrecisionError("precision must be >= 1", required=1, available=prec)
        self._prec = prec
        if extension is None:
            self._extension = ExtensionData(p=p, n=1, e=1, f=1)
        else:
            if extension.p != p:
                raise DomainError("extension prime must match local field prime")
            self._extension = extension

    @property
    def p(self) -> int:
        return self._p

    @property
    def prec(self) -> int:
        return self._prec

    @property
    def extension(self) -> ExtensionData:
        return self._extension

    def with_precision(self, prec: int) -> LocalField:
        return LocalField(self._p, prec, extension=self._extension)

    def unramified(self, f: int) -> LocalField:
        if f < 1:
            raise DomainError("f must be >= 1")
        return LocalField(self._p, self._prec, extension=ExtensionData(p=self._p, n=f, e=1, f=f))

    def totally_ramified(self, e: int) -> LocalField:
        if e < 1:
            raise DomainError("e must be >= 1")
        return LocalField(self._p, self._prec, extension=ExtensionData(p=self._p, n=e, e=e, f=1))

    def int(self, value: int) -> PadicInt:
        return PadicInt(value, self._p, self._prec)

    def number(self, value: int | Fraction) -> PadicNumber:
        return PadicNumber.from_rational(value, self._p, self._prec)

    def zero(self) -> PadicNumber:
        return PadicNumber.zero(self._p, self._prec)

    def one(self) -> PadicNumber:
        return PadicNumber.one(self._p, self._prec)

    def uniformizer(self) -> PadicNumber:
        return PadicNumber.p_element(self._p, self._prec)

    def teichmuller(self, a: int) -> PadicInt:
        from padicfield.teichmuller import teichmuller_representative

        return teichmuller_representative(a, self._p, self._prec)

    def sqrt(self, x: int | Fraction | PadicNumber) -> PadicNumber:
        from padicfield.sqrt import sqrt_padic_number, sqrt_rational

        if isinstance(x, PadicNumber):
            return sqrt_padic_number(x.with_precision(self._prec))
        return sqrt_rational(x, self._p, self._prec)

    def exp(self, x: PadicNumber) -> PadicNumber:
        from padicfield.exp_log import padic_exp

        return padic_exp(x.with_precision(self._prec))

    def log(self, u: PadicNumber | PadicInt) -> PadicNumber:
        from padicfield.exp_log import padic_log

        return padic_log(u, prec=self._prec)

    def valuation(self, x: int | Fraction | PadicNumber | PadicInt) -> float:
        from padicfield.padic_number import padic_valuation_of

        return padic_valuation_of(x, self._p)

    def __repr__(self) -> str:
        ext = self._extension
        if ext.is_trivial:
            return f"LocalField(Q_{self._p}, prec={self._prec})"
        return (
            f"LocalField(extension/Q_{self._p}, n={ext.n}, e={ext.e}, f={ext.f}, "
            f"prec={self._prec})"
        )


def qp(p: int, prec: int = DEFAULT_PRECISION) -> LocalField:
    return LocalField(p, prec)


def describe_unit_group(p: int, k: int) -> str:
    """Structure of ``(Z/p^k Z)^×`` for small cases."""
    require_prime(p)
    if k < 1:
        raise DomainError("k must be >= 1")
    if p == 2:
        if k == 1:
            return "trivial group"
        if k == 2:
            return "C_2"
        return f"C_2 × C_{2 ** (k - 2)}"
    return f"C_{(p - 1) * (p ** (k - 1))}"


def tame_ramification(p: int, e: int) -> bool:
    require_prime(p)
    return e % p != 0
