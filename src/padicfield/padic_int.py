"""p-adic integers ``Z_p`` at finite precision."""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Iterator, Sequence

from padicfield.exceptions import IncompatiblePrimeError, PrecisionError, ValuationError
from padicfield.utils import digits_base_p, from_digits_base_p, inv_mod_pk, require_prime
from padicfield.valuation import valuation_int

DEFAULT_PRECISION = 20


class PadicInt:
    """Element of ``Z / p^N Z``, viewed as a truncated p-adic integer."""

    __slots__ = ("_p", "_prec", "_value", "_mod")

    def __init__(self, value: int, p: int, prec: int = DEFAULT_PRECISION) -> None:
        self._p = require_prime(p)
        if prec < 1:
            raise PrecisionError("precision must be >= 1", required=1, available=prec)
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"value must be int, got {type(value).__name__}")
        self._prec = prec
        self._mod = p**prec
        self._value = value % self._mod

    @classmethod
    def zero(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
        return cls(0, p, prec)

    @classmethod
    def one(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
        return cls(1, p, prec)

    @classmethod
    def from_digits(cls, digits: Sequence[int], p: int) -> PadicInt:
        require_prime(p)
        if not digits:
            raise ValueError("digits must be non-empty")
        for d in digits:
            if not isinstance(d, int) or d < 0 or d >= p:
                raise ValueError(f"invalid digit {d} for prime {p}")
        return cls(from_digits_base_p(digits, p), p, prec=len(digits))

    @classmethod
    def from_rational(cls, x: int | Fraction, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
        """Embed a p-integral rational into ``Z_p``."""
        require_prime(p)
        frac = Fraction(x)
        if frac == 0:
            return cls.zero(p, prec)
        from padicfield.valuation import valuation_int as vint

        if vint(frac.numerator, p) - vint(frac.denominator, p) < 0:
            raise ValuationError(f"{x} is not a p-adic integer (negative valuation)")
        mod = p**prec
        inv_den = inv_mod_pk(frac.denominator % mod, p, prec)
        return cls((frac.numerator * inv_den) % mod, p, prec)

    @classmethod
    def teichmuller_approx(cls, a: int, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
        from padicfield.teichmuller import teichmuller_representative

        return teichmuller_representative(a, p, prec)

    @property
    def p(self) -> int:
        return self._p

    @property
    def prec(self) -> int:
        return self._prec

    @property
    def modulus(self) -> int:
        return self._mod

    @property
    def residue(self) -> int:
        return self._value

    @property
    def value(self) -> int:
        return self._value

    def digits(self) -> list[int]:
        return digits_base_p(self._value, self._p, length=self._prec)

    def digit(self, i: int) -> int:
        if i < 0 or i >= self._prec:
            raise IndexError(f"digit index {i} out of range [0, {self._prec})")
        return (self._value // (self._p**i)) % self._p

    def valuation(self) -> int | float:
        if self._value == 0:
            return float(self._prec)
        return valuation_int(self._value, self._p)

    def unit_part(self) -> PadicInt:
        if self._value == 0:
            raise ValuationError("unit part of zero is undefined")
        v = valuation_int(self._value, self._p)
        if v >= self._prec:
            raise ValuationError("element is indistinguishable from zero at this precision")
        return PadicInt(self._value // (self._p**v), self._p, self._prec - v)

    def is_unit(self) -> bool:
        return self._value % self._p != 0

    def is_zero(self) -> bool:
        return self._value == 0

    def lift_precision(self, new_prec: int) -> PadicInt:
        if new_prec < self._prec:
            raise PrecisionError(
                f"cannot lift to lower precision {new_prec}",
                required=self._prec,
                available=new_prec,
            )
        return PadicInt(self._value, self._p, new_prec)

    def with_precision(self, new_prec: int) -> PadicInt:
        if new_prec < 1:
            raise PrecisionError("precision must be >= 1", required=1, available=new_prec)
        return PadicInt(self._value, self._p, new_prec)

    def _check_prime(self, other: PadicInt) -> None:
        if self._p != other._p:
            raise IncompatiblePrimeError(self._p, other._p)

    def _align(self, other: PadicInt) -> tuple[PadicInt, PadicInt, int]:
        self._check_prime(other)
        prec = min(self._prec, other._prec)
        return self.with_precision(prec), other.with_precision(prec), prec

    def _coerce(self, other: PadicInt | int) -> PadicInt:
        if isinstance(other, PadicInt):
            return other
        if isinstance(other, int) and not isinstance(other, bool):
            return PadicInt(other, self._p, self._prec)
        raise TypeError(f"cannot coerce {type(other).__name__} to PadicInt")

    def __add__(self, other: PadicInt | int) -> PadicInt:
        other = self._coerce(other)
        a, b, prec = self._align(other)
        return PadicInt(a._value + b._value, self._p, prec)

    def __radd__(self, other: int) -> PadicInt:
        return self.__add__(other)

    def __sub__(self, other: PadicInt | int) -> PadicInt:
        other = self._coerce(other)
        a, b, prec = self._align(other)
        return PadicInt(a._value - b._value, self._p, prec)

    def __rsub__(self, other: int) -> PadicInt:
        return self._coerce(other).__sub__(self)

    def __neg__(self) -> PadicInt:
        return PadicInt(-self._value, self._p, self._prec)

    def __pos__(self) -> PadicInt:
        return PadicInt(self._value, self._p, self._prec)

    def __mul__(self, other: PadicInt | int) -> PadicInt:
        other = self._coerce(other)
        a, b, prec = self._align(other)
        return PadicInt(a._value * b._value, self._p, prec)

    def __rmul__(self, other: int) -> PadicInt:
        return self.__mul__(other)

    def divide(self, other: PadicInt | int) -> PadicInt:
        other = self._coerce(other)
        self._check_prime(other)
        ov = other._value
        if ov == 0:
            raise ZeroDivisionError("p-adic division by zero")
        v = valuation_int(ov, other._p)
        prec = min(self._prec, other._prec)
        if v >= prec:
            raise ZeroDivisionError("divisor is zero at working precision")
        a = self._value % (self._p**prec)
        b = ov % (other._p**prec)
        va = 0 if a == 0 else valuation_int(a, self._p)
        if a != 0 and va < v:
            raise ValuationError("quotient would not be a p-adic integer")
        new_prec = prec - v
        a_unit = a // (self._p**v) if v else a
        b_unit = b // (self._p**v)
        return PadicInt(a_unit * inv_mod_pk(b_unit, self._p, new_prec), self._p, new_prec)

    def __truediv__(self, other: PadicInt | int) -> PadicInt:
        return self.divide(other)

    def __rtruediv__(self, other: int) -> PadicInt:
        return self._coerce(other).divide(self)

    def __floordiv__(self, other: PadicInt | int) -> PadicInt:
        return self.divide(other)

    def __pow__(self, exp: int, modulo: int | None = None) -> PadicInt:
        if modulo is not None:
            raise TypeError("modular pow with three arguments not supported for PadicInt")
        if not isinstance(exp, int) or isinstance(exp, bool):
            raise TypeError("exponent must be int")
        if exp < 0:
            return self.inverse() ** (-exp)
        return PadicInt(pow(self._value, exp, self._mod), self._p, self._prec)

    def inverse(self) -> PadicInt:
        if not self.is_unit():
            raise ValuationError("only units are invertible in Z_p without precision loss")
        return PadicInt(inv_mod_pk(self._value, self._p, self._prec), self._p, self._prec)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, int) and not isinstance(other, bool):
            return self._value == (other % self._mod)
        if not isinstance(other, PadicInt):
            return NotImplemented
        if self._p != other._p:
            return False
        prec = min(self._prec, other._prec)
        mod = self._p**prec
        return (self._value % mod) == (other._value % mod)

    def __hash__(self) -> int:
        return hash((self._p, self._prec, self._value))

    def congruent(self, other: PadicInt | int, k: int) -> bool:
        if k < 0:
            raise ValueError("k must be non-negative")
        if k == 0:
            return True
        if k > self._prec:
            raise PrecisionError(
                f"insufficient precision to test congruence mod p^{k}",
                required=k,
                available=self._prec,
            )
        other_val = other._value if isinstance(other, PadicInt) else int(other)
        if isinstance(other, PadicInt):
            self._check_prime(other)
            if k > other._prec:
                raise PrecisionError(
                    f"insufficient precision on other",
                    required=k,
                    available=other._prec,
                )
        mod = self._p**k
        return (self._value % mod) == (other_val % mod)

    def to_int(self, *, centered: bool = False) -> int:
        if not centered:
            return self._value
        half = self._mod // 2
        return self._value - self._mod if self._value > half else self._value

    def to_rational_approx(self) -> Fraction:
        return Fraction(self._value)

    def __repr__(self) -> str:
        return f"PadicInt({self._value}, p={self._p}, prec={self._prec})"

    def __str__(self) -> str:
        body = " ".join(str(d) for d in reversed(self.digits()))
        return f"... {body} (p={self._p}, O({self._p}^{self._prec}))"

    def latex(self) -> str:
        return rf"{self._value} + O({self._p}^{{{self._prec}}})"

    def expand(self) -> Iterator[tuple[int, int]]:
        for i, d in enumerate(self.digits()):
            yield i, d


def padic_sum(elements: Iterable[PadicInt], *, p: int | None = None, prec: int | None = None) -> PadicInt:
    it = iter(elements)
    try:
        acc = next(it)
    except StopIteration:
        if p is None or prec is None:
            raise ValueError("empty sum requires p and prec") from None
        return PadicInt.zero(p, prec)
    for el in it:
        acc = acc + el
    return acc


def padic_product(
    elements: Iterable[PadicInt], *, p: int | None = None, prec: int | None = None
) -> PadicInt:
    it = iter(elements)
    try:
        acc = next(it)
    except StopIteration:
        if p is None or prec is None:
            raise ValueError("empty product requires p and prec") from None
        return PadicInt.one(p, prec)
    for el in it:
        acc = acc * el
    return acc
