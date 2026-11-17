"""Tests for Hensel lifting, Teichmüller representatives, and square roots."""

from __future__ import annotations

import pytest

from padicfield.exceptions import HenselError, SquareRootError, TeichmullerError
from padicfield.hensel import (
    eval_poly,
    find_roots_mod_p,
    hensel_lift,
    hensel_lift_simple,
    hensel_lift_step,
    is_simple_root,
    lift_all_simple_roots,
    multiplicity_mod_p,
    verify_root,
)
from padicfield.padic_int import PadicInt
from padicfield.sqrt import (
    is_square_in_qp,
    is_square_in_zp,
    sqrt_mod_p,
    sqrt_padic_int,
    sqrt_rational,
    square_root_exists_rational,
    tonelli_shanks,
)
from padicfield.teichmuller import (
    all_teichmuller,
    decompose_unit,
    is_teichmuller,
    multiplicative_order_teichmuller,
    teichmuller_power,
    teichmuller_representative,
)


def test_find_roots_and_simple():
    # x^2 - 2 over F_7: roots 3, 4
    coeffs = [-2, 0, 1]
    roots = find_roots_mod_p(coeffs, 7)
    assert set(roots) == {3, 4}
    assert is_simple_root(coeffs, 3, 7)


def test_hensel_lift_sqrt2_mod7():
    coeffs = [-2, 0, 1]
    result = hensel_lift_simple(coeffs, 3, 7, 6)
    assert verify_root(coeffs, result.root)
    assert result.simple
    r = result.root.residue
    assert (r * r - 2) % (7**6) == 0


def test_hensel_lift_step_and_all():
    coeffs = [-2, 0, 1]
    a1 = hensel_lift_step(coeffs, 3, 7, 1)
    assert eval_poly(coeffs, a1, 49) % 49 == 0
    lifted = lift_all_simple_roots(coeffs, 7, 5)
    assert len(lifted) == 2


def test_hensel_rejects_nonsimple():
    # (x-1)^2 = x^2 - 2x + 1 over F_5
    coeffs = [1, -2, 1]
    with pytest.raises(HenselError):
        hensel_lift(coeffs, 1, 5, 4)


def test_multiplicity():
    assert multiplicity_mod_p([1, -2, 1], 1, 5) == 2


def test_teichmuller_basic():
    w = teichmuller_representative(3, 7, 6)
    assert w.residue % 7 == 3
    assert (w ** 6).congruent(1, 6)
    assert is_teichmuller(w)
    assert teichmuller_representative(0, 5, 4).is_zero()
    assert teichmuller_representative(1, 2, 5).residue == 1


def test_teichmuller_power_and_decompose():
    assert teichmuller_power(3, 2, 7, 5) == teichmuller_representative(2, 7, 5)  # 3^2=9≡2
    u = PadicInt(17, 5, 6)
    omega, principal = decompose_unit(u)
    assert is_teichmuller(omega)
    assert principal.residue % 5 == 1
    assert multiplicative_order_teichmuller(3, 7) == 6
    assert len(all_teichmuller(5, 4)) == 5


def test_tonelli_and_sqrt_zp():
    assert tonelli_shanks(2, 7) in (3, 4)
    assert sqrt_mod_p(2, 7) ** 2 % 7 == 2
    a = PadicInt(2, 7, 8)
    assert is_square_in_zp(a)
    r = sqrt_padic_int(a)
    assert (r * r).congruent(a, r.prec)


def test_sqrt_qp_and_rational():
    s = sqrt_rational(2, 7, 8)
    assert is_square_in_qp(s * s)
    assert square_root_exists_rational(2, 7)
    assert not square_root_exists_rational(3, 7)
    with pytest.raises(SquareRootError):
        sqrt_rational(3, 7, 6)


def test_sqrt_2adic():
    # 17 ≡ 1 mod 8, so square in Z_2
    a = PadicInt(17, 2, 8)
    assert is_square_in_zp(a)
    r = sqrt_padic_int(a)
    assert (r * r).congruent(a, min(r.prec, a.prec))
