"""Tests for utils and valuation."""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from padicfield.exceptions import NotPrimeError, ValuationError
from padicfield.utils import (
    binomial,
    crt,
    digits_base_p,
    factorize,
    from_digits_base_p,
    gcd,
    inv_mod_pk,
    is_prime,
    jacobi_symbol,
    kronecker_symbol,
    legendre_symbol,
    modinv,
    next_prime,
    primes_upto,
    require_prime,
)
from padicfield.valuation import (
    ValuationRing,
    gauss_valuation,
    kummer_carry_count,
    padic_norm,
    remove_factor,
    unit_part_int,
    valuation,
    valuation_int,
    valuation_of_binomial,
    valuation_of_factorial,
)

def test_is_prime_basic():
    assert is_prime(2) and is_prime(3) and is_prime(97)
    assert not is_prime(1) and not is_prime(0) and not is_prime(91)
    assert require_prime(5) == 5
    with pytest.raises(NotPrimeError):
        require_prime(9)

def test_primes_and_factorize():
    assert primes_upto(20) == [2, 3, 5, 7, 11, 13, 17, 19]
    assert next_prime(14) == 17
    assert factorize(360) == {2: 3, 3: 2, 5: 1}
    assert factorize(1) == {}

def test_gcd_modinv_crt():
    assert gcd(24, 36) == 12
    assert modinv(3, 7) == 5
    assert crt([2, 3, 2], [3, 5, 7]) == 23

def test_symbols():
    assert legendre_symbol(2, 7) == 1
    assert legendre_symbol(3, 7) == -1
    assert jacobi_symbol(2, 15) == jacobi_symbol(2, 3) * jacobi_symbol(2, 5)
    assert kronecker_symbol(2, 5) == -1
    assert kronecker_symbol(-1, 5) == 1

def test_digits_and_inv_pk():
    digs = digits_base_p(42, 5, length=4)
    assert from_digits_base_p(digs, 5) == 42
    assert (11 * inv_mod_pk(11, 5, 4)) % (5**4) == 1

def test_binomial():
    assert binomial(10, 3) == 120
    assert binomial(5, 0) == 1

def test_valuation_int_and_rational():
    assert valuation_int(360, 5) == 1
    assert valuation(Fraction(12, 50), 5) == -2  # 12/50 = 6/25
    assert math.isinf(valuation(0, 5))
    assert unit_part_int(-50, 5) == -2
    assert remove_factor(360, 5) == (1, 72)
    with pytest.raises(ValuationError):
        valuation_int(0, 5)

def test_padic_norm_and_gauss():
    assert padic_norm(25, 5) == pytest.approx(0.04)
    assert padic_norm(0, 5) == 0.0
    assert gauss_valuation([0, 25, 125], 5) == 2

def test_factorial_binomial_kummer():
    assert valuation_of_factorial(10, 5) == 2
    assert valuation_of_binomial(10, 3, 5) == kummer_carry_count(10, 3, 5)

def test_valuation_ring():
    R = ValuationRing(5)
    assert R.contains(Fraction(2, 3))
    assert not R.contains(Fraction(1, 5))
    assert R.reduce_mod_p(Fraction(2, 3)) == (2 * modinv(3, 5)) % 5
