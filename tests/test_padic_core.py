"""Tests for PadicInt and PadicNumber."""

from __future__ import annotations

from fractions import Fraction

import pytest

from padicfield.exceptions import IncompatiblePrimeError, ValuationError
from padicfield.padic_int import PadicInt, padic_product, padic_sum
from padicfield.padic_number import PadicNumber, embed_rational, padic_valuation_of

def test_padic_int_arithmetic():
    a = PadicInt(17, 5, 6)
    b = PadicInt(3, 5, 6)
    assert (a + b).residue == (17 + 3) % (5**6)
    assert (a * b).residue == (17 * 3) % (5**6)
    assert (a - b).residue == (17 - 3) % (5**6)
    assert (-a).residue == (-17) % (5**6)

def test_padic_int_digits_and_unit():
    x = PadicInt(42, 5, 5)
    assert x.digits()[0] == 42 % 5
    assert from_digits_roundtrip(x)
    assert PadicInt(10, 5, 4).unit_part().residue % 5 == 2
    assert PadicInt(7, 5, 4).is_unit()
    assert not PadicInt(10, 5, 4).is_unit()

def from_digits_roundtrip(x: PadicInt) -> bool:
    return PadicInt.from_digits(x.digits(), x.p) == x

def test_padic_int_inverse_and_div():
    u = PadicInt(3, 7, 5)
    assert (u * u.inverse()).congruent(1, 5)
    a = PadicInt(15, 5, 6)
    b = PadicInt(5, 5, 6)
    q = a / b
    assert (q * b).congruent(a, q.prec)

def test_padic_int_pow_and_eq():
    a = PadicInt(2, 5, 6)
    assert (a**3).residue == pow(2, 3, 5**6)
    assert a == PadicInt(2 + 5**6, 5, 6)
    assert a == 2
    with pytest.raises(IncompatiblePrimeError):
        _ = PadicInt(1, 5, 3) + PadicInt(1, 7, 3)

def test_from_rational_integral():
    x = PadicInt.from_rational(Fraction(2, 3), 5, 6)
    assert (x * 3).congruent(2, 6)
    with pytest.raises(ValuationError):
        PadicInt.from_rational(Fraction(1, 5), 5, 4)

def test_padic_number_split_and_arith():
    x = PadicNumber.from_int(50, 5, 6)  # 5^2 * 2
    assert x.val == 2
    assert x.unit.residue % 5 == 2
    y = PadicNumber.from_rational(Fraction(2, 5), 5, 6)
    assert y.val == -1
    s = x + y
    assert isinstance(s, PadicNumber)
    assert (x * y).val == 1

def test_padic_number_inverse_pow():
    u = PadicNumber.from_int(3, 7, 8)
    assert (u * u.inverse()).is_unit()
    assert (u**2).unit.residue % 7 == 2
    z = PadicNumber.zero(5, 4)
    assert z.is_zero() and z.is_integral()

def test_embed_and_to_int():
    x = embed_rational(Fraction(7, 3), 5, 6)
    assert padic_valuation_of(x, 5) == 0
    pi = PadicNumber.from_int(25, 5, 5).to_padic_int()
    assert pi.residue % 25 == 0

def test_padic_sum_product():
    els = [PadicInt(i, 5, 4) for i in range(1, 5)]
    assert padic_sum(els).residue == sum(range(1, 5)) % (5**4)
    assert padic_product(els).residue == (24) % (5**4)

def test_precision_align():
    a = PadicInt(1, 5, 8)
    b = PadicInt(1, 5, 3)
    assert (a + b).prec == 3

def test_str_repr():
    x = PadicInt(3, 5, 4)
    assert "PadicInt" in repr(x)
    assert "p=5" in str(x)
    n = PadicNumber.from_int(5, 5, 4)
    assert "PadicNumber" in repr(n)
