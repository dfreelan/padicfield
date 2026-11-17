"""p-adic numbers ``Q_p`` at finite relative precision."""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Union

from padicfield.exceptions import IncompatiblePrimeError, PrecisionError, ValuationError
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.utils import inv_mod_pk, require_prime
from padicfield.valuation import remove_factor, valuation_int

NumberLike = Union[int, Fraction, "PadicNumber", PadicInt]


class PadicNumber:
    """Element of ``Q_p`` stored as ``p^{val} * unit`` with unit of valuation 0."""

    __slots__ = ("_p", "_val", "_unit", "_prec")

    def __init__(
        self,
        unit: PadicInt | int,
        p: int | None = None,
        *,
        val: int = 0,
        prec: int | None = None,
    ) -> None:
        if isinstance(unit, PadicInt):
            self._p = unit.p
            if p is not None and p != unit.p:
                raise IncompatiblePrimeError(p, unit.p)
            u = unit if prec is None else unit.with_precision(prec)
        else:
            if p is None:
                raise ValueError("prime p is required when unit is an int")
            self._p = require_prime(p)
            u = PadicInt(int(unit), self._p, prec if prec is not None else DEFAULT_PRECISION)

        self._prec = u.prec
        if u.is_zero():
            self._val: int | None = None
            self._unit = PadicInt.zero(self._p, self._prec)
            return
        uv = valuation_int(u.residue, self._p)
        if uv >= u.prec:
            self._val = None
            self._unit = PadicInt.zero(self._p, self._prec)
        elif uv > 0:
            self._val = val + uv
            self._unit = u.unit_part().with_precision(max(1, u.prec - uv))
            self._prec = self._unit.prec
        else:
            self._val = val
            self._unit = u

    @classmethod
    def zero(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
        return cls(PadicInt.zero(p, prec), val=0)

    @classmethod
    def one(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
        return cls(PadicInt.one(p, prec), val=0)

    @classmethod
    def p_element(cls, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
        return cls(PadicInt.one(p, prec), val=1)

    @classmethod
    def from_int(cls, n: int, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
        require_prime(p)
        if n == 0:
            return cls.zero(p, prec)
        v, u = remove_factor(n, p)
        return cls(PadicInt(u, p, prec), val=v)

    @classmethod
    def from_rational(cls, x: int | Fraction, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
        require_prime(p)
        frac = Fraction(x)
        if frac == 0:
            return cls.zero(p, prec)
        v_num, u_num = remove_factor(frac.numerator, p)
        v_den, u_den = remove_factor(frac.denominator, p)
        inv = inv_mod_pk(u_den % (p**prec), p, prec)
        return cls(PadicInt((u_num * inv) % (p**prec), p, prec), val=v_num - v_den)

    @classmethod
    def from_padic_int(cls, x: PadicInt, *, val: int = 0) -> PadicNumber:
        return cls(x, val=val)

    @property
    def p(self) -> int:
        return self._p

    @property
    def prec(self) -> int:
        return self._prec

    @property
    def valuation(self) -> float:
        return math.inf if self._val is None else float(self._val)

    @property
    def val(self) -> int | None:
        return self._val

    @property
    def unit(self) -> PadicInt:
        return self._unit

    def is_zero(self) -> bool:
        return self._val is None

    def is_integral(self) -> bool:
        return self._val is None or self._val >= 0

    def is_unit(self) -> bool:
        return self._val == 0

    def with_precision(self, new_prec: int) -> PadicNumber:
        if new_prec < 1:
            raise PrecisionError("precision must be >= 1", required=1, available=new_prec)
        if self.is_zero():
            return PadicNumber.zero(self._p, new_prec)
        return PadicNumber(self._unit.with_precision(new_prec), val=self._val or 0)

    def absolute_precision(self) -> float:
        return math.inf if self._val is None else float(self._val + self._prec)

    def _coerce(self, other: NumberLike) -> PadicNumber:
        if isinstance(other, PadicNumber):
            if other.p != self._p:
                raise IncompatiblePrimeError(self._p, other.p)
            return other
        if isinstance(other, PadicInt):
            if other.p != self._p:
                raise IncompatiblePrimeError(self._p, other.p)
            return PadicNumber.from_padic_int(other)
        if isinstance(other, Fraction):
            return PadicNumber.from_rational(other, self._p, self._prec)
        if isinstance(other, int) and not isinstance(other, bool):
            return PadicNumber.from_int(other, self._p, self._prec)
        raise TypeError(f"cannot coerce {type(other).__name__} to PadicNumber")

    def __add__(self, other: NumberLike) -> PadicNumber:
        other = self._coerce(other)
        if self.is_zero():
            return other.with_precision(min(other.prec, self._prec) if not other.is_zero() else self._prec)
        if other.is_zero():
            return self.with_precision(min(self._prec, other.prec))
        assert self._val is not None and other._val is not None
        v = min(self._val, other._val)
        abs_result = min(self._val + self._prec, other._val + other._prec)
        rel_result = abs_result - v
        if rel_result < 1:
            return PadicNumber.zero(self._p, 1)
        mod = self._p**rel_result
        shift_s, shift_o = self._val - v, other._val - v
        u_s = self._unit.with_precision(max(self._unit.prec, rel_result - shift_s)).residue
        u_o = other._unit.with_precision(max(other._unit.prec, rel_result - shift_o)).residue
        s = (u_s * (self._p**shift_s) + u_o * (self._p**shift_o)) % mod
        if s == 0:
            return PadicNumber.zero(self._p, max(1, rel_result))
        return PadicNumber(PadicInt(s, self._p, rel_result), val=v)

    def __radd__(self, other: NumberLike) -> PadicNumber:
        return self.__add__(other)

    def __sub__(self, other: NumberLike) -> PadicNumber:
        return self + (-self._coerce(other))

    def __rsub__(self, other: NumberLike) -> PadicNumber:
        return self._coerce(other) - self

    def __neg__(self) -> PadicNumber:
        if self.is_zero():
            return PadicNumber.zero(self._p, self._prec)
        return PadicNumber(-self._unit, val=self._val or 0)

    def __pos__(self) -> PadicNumber:
        return self.with_precision(self._prec)

    def __mul__(self, other: NumberLike) -> PadicNumber:
        other = self._coerce(other)
        if self.is_zero() or other.is_zero():
            return PadicNumber.zero(self._p, min(self._prec, other.prec))
        assert self._val is not None and other._val is not None
        prec = min(self._prec, other.prec)
        unit = self._unit.with_precision(prec) * other._unit.with_precision(prec)
        return PadicNumber(unit, val=self._val + other._val)

    def __rmul__(self, other: NumberLike) -> PadicNumber:
        return self.__mul__(other)

    def __truediv__(self, other: NumberLike) -> PadicNumber:
        other = self._coerce(other)
        if other.is_zero():
            raise ZeroDivisionError("p-adic division by zero")
        if self.is_zero():
            return PadicNumber.zero(self._p, min(self._prec, other.prec))
        assert self._val is not None and other._val is not None
        prec = min(self._prec, other.prec)
        unit = self._unit.with_precision(prec) * other._unit.with_precision(prec).inverse()
        return PadicNumber(unit, val=self._val - other._val)

    def __rtruediv__(self, other: NumberLike) -> PadicNumber:
        return self._coerce(other) / self

    def inverse(self) -> PadicNumber:
        if self.is_zero():
            raise ZeroDivisionError("p-adic inversion of zero")
        assert self._val is not None
        return PadicNumber(self._unit.inverse(), val=-self._val)

    def __pow__(self, exp: int) -> PadicNumber:
        if not isinstance(exp, int) or isinstance(exp, bool):
            raise TypeError("exponent must be int")
        if exp == 0:
            return PadicNumber.one(self._p, self._prec)
        if self.is_zero():
            if exp < 0:
                raise ZeroDivisionError("negative power of zero")
            return PadicNumber.zero(self._p, self._prec)
        assert self._val is not None
        if exp < 0:
            return self.inverse() ** (-exp)
        return PadicNumber(self._unit**exp, val=self._val * exp)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            try:
                other = self._coerce(other)
            except (TypeError, ValuationError):
                return False
        if isinstance(other, PadicInt):
            other = PadicNumber.from_padic_int(other)
        if not isinstance(other, PadicNumber):
            return NotImplemented
        if self._p != other._p:
            return False
        if self.is_zero() and other.is_zero():
            return True
        if self.is_zero() or other.is_zero() or self._val != other._val:
            return False
        prec = min(self._prec, other.prec)
        return self._unit.with_precision(prec) == other._unit.with_precision(prec)

    def __hash__(self) -> int:
        return hash((self._p, self._val, self._unit.residue, self._prec))

    def congruent(self, other: NumberLike, abs_prec: int) -> bool:
        diff = self - other
        if diff.is_zero():
            return True
        assert diff._val is not None
        return diff._val >= abs_prec

    def to_padic_int(self) -> PadicInt:
        if not self.is_integral():
            raise ValuationError("not a p-adic integer")
        if self.is_zero():
            return PadicInt.zero(self._p, self._prec)
        assert self._val is not None
        new_prec = self._val + self._prec
        residue = (self._unit.residue * (self._p**self._val)) % (self._p**new_prec)
        return PadicInt(residue, self._p, new_prec)

    def to_rational_approx(self) -> Fraction:
        if self.is_zero():
            return Fraction(0)
        assert self._val is not None
        if self._val >= 0:
            return Fraction(self._unit.residue * (self._p**self._val))
        return Fraction(self._unit.residue, self._p ** (-self._val))

    def residue_field_image(self) -> int:
        if self.is_zero() or (self._val is not None and self._val > 0):
            return 0
        if self._val is not None and self._val < 0:
            raise ValuationError("element is not integral; cannot reduce to residue field")
        return self._unit.residue % self._p

    def norm(self) -> float:
        if self.is_zero():
            return 0.0
        assert self._val is not None
        return float(self._p) ** (-self._val)

    def __repr__(self) -> str:
        if self.is_zero():
            return f"PadicNumber(0, p={self._p}, prec={self._prec})"
        return f"PadicNumber(val={self._val}, unit={self._unit.residue}, p={self._p}, prec={self._prec})"

    def __str__(self) -> str:
        if self.is_zero():
            return f"0 + O({self._p}^{self._prec})"
        assert self._val is not None
        return (
            f"{self._p}^{self._val} * {self._unit.residue}"
            f" + O({self._p}^{self._val + self._prec})"
        )

    def latex(self) -> str:
        if self.is_zero():
            return rf"0 + O({self._p}^{{{self._prec}}})"
        assert self._val is not None
        return (
            rf"{self._p}^{{{self._val}}}\cdot {self._unit.residue}"
            rf" + O({self._p}^{{{self._val + self._prec}}})"
        )


def embed_rational(x: int | Fraction, p: int, prec: int = DEFAULT_PRECISION) -> PadicNumber:
    return PadicNumber.from_rational(x, p, prec)


def padic_valuation_of(x: NumberLike, p: int) -> float:
    if isinstance(x, PadicNumber):
        if x.p != p:
            raise IncompatiblePrimeError(x.p, p)
        return x.valuation
    if isinstance(x, PadicInt):
        if x.p != p:
            raise IncompatiblePrimeError(x.p, p)
        return float(x.valuation()) if not x.is_zero() else math.inf
    from padicfield.valuation import valuation

    return valuation(x, p)
