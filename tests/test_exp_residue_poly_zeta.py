"""Tests for exp/log, residue symbols, polynomials, zeta, local fields, CLI."""

from __future__ import annotations

from fractions import Fraction

import pytest

from padicfield.cli import main
from padicfield.exceptions import ConvergenceError, DomainError
from padicfield.exp_log import (
    exp_converges,
    exp_log_roundtrip_ok,
    log_converges,
    padic_exp,
    padic_log,
)
from padicfield.local_field import ExtensionData, LocalField, describe_unit_group, qp
from padicfield.padic_int import PadicInt
from padicfield.padic_number import PadicNumber
from padicfield.polynomials import (
    IntPolynomial,
    PadicPolynomial,
    compose,
    discriminant,
    resultant_univariate,
)
from padicfield.residue import (
    hilbert_symbol,
    is_square_mod_pk,
    quadratic_residue_symbol,
    residue_symbol_mod_pk,
    square_classes_qp,
)
from padicfield.zeta import (
    count_units_mod_pk,
    count_zeros_univariate_mod_pk,
    density_of_units,
    euler_factor_affine_line,
    igusa_naive_poincare,
    local_zeta_motive_guess,
    projective_line_point_count,
)


def test_exp_log_principal_unit():
    p, prec = 5, 10
    x = PadicNumber.from_int(p, p, prec) * 3  # valuation 1 > 1/(p-1)
    assert exp_converges(x)
    e = padic_exp(x)
    assert log_converges(e)
    y = padic_log(e)
    assert (y - x).is_zero() or (y - x).val >= prec - 3  # type: ignore[operator]
    assert exp_log_roundtrip_ok(x)


def test_log_of_padic_int():
    u = PadicInt(1 + 2 * 5, 5, 8)
    assert log_converges(u)
    L = padic_log(u)
    assert not L.is_zero() or u.congruent(1, 8)


def test_exp_rejects_bad_disk():
    x = PadicNumber.from_int(1, 5, 6)  # valuation 0
    assert not exp_converges(x)
    with pytest.raises(ConvergenceError):
        padic_exp(x)


def test_legendre_and_residue_pk():
    assert quadratic_residue_symbol(2, 7) == 1
    assert residue_symbol_mod_pk(2, 7, 3) == 1
    assert residue_symbol_mod_pk(3, 7, 3) == -1
    assert is_square_mod_pk(2, 7, 4)
    assert not is_square_mod_pk(3, 7, 4)
    assert is_square_mod_pk(0, 5, 3)


def test_hilbert_symbol():
    assert hilbert_symbol(1, 1, 5) == 1
    assert hilbert_symbol(2, 5, 7) in (-1, 1)
    # (2,2)_2 = 1? Actually (2,2)_2 = -1? Check known value: (2,2)_2 = 1? 
    # Serre: (2,2)_2 = (-1)^{ω(2)ω(2)?} wait units... 2 = 2^1 * 1, so α=1,u=1
    # (2,2)_2 = (-1)^{ε(1)ε(1) + 1·ω(1) + 1·ω(1)} = 1
    assert hilbert_symbol(2, 2, 2) == 1
    assert len(square_classes_qp(7)) == 4
    assert len(square_classes_qp(2)) == 8


def test_int_polynomial_ops():
    f = IntPolynomial([-2, 0, 1])  # x^2 - 2
    g = IntPolynomial([-1, 1])  # x - 1
    assert (f + g).coeffs[0] == -3
    assert (f * g).degree == 3
    assert f.evaluate(3, 7) == 0
    assert f.derivative().coeffs == [0, 2]
    assert set(f.roots_mod_p(7)) == {3, 4}
    root = f.hensel_lift_root(3, 7, 5)
    assert (root * root).congruent(2, 5)


def test_padic_polynomial_and_resultant():
    f = IntPolynomial.from_coeffs(-2, 0, 1)
    pf = PadicPolynomial.from_int_poly(f, 7, 5)
    assert pf.evaluate(3).residue % 7 == 0
    # Res(x-1, x-2) = 2-1 = 1? Res(x-a,x-b)=b-a
    a = IntPolynomial([-1, 1])
    b = IntPolynomial([-2, 1])
    assert abs(resultant_univariate(a, b)) == 1
    # disc(x^2-2) = 8
    assert discriminant(f) == 8
    h = compose(IntPolynomial([0, 0, 1]), IntPolynomial([1, 1]))  # (1+x)^2
    assert h.coeffs == [1, 2, 1]


def test_zeta_counts():
    lc = count_zeros_univariate_mod_pk([-2, 0, 1], 7, 2)
    assert lc.count >= 2
    assert count_units_mod_pk(5, 3) == 5**3 - 5**2
    assert projective_line_point_count(5) == 6
    assert density_of_units(5) == Fraction(4, 5)
    dens = igusa_naive_poincare([-2, 0, 1], 7, 2)
    assert len(dens) == 2
    guess = local_zeta_motive_guess(IntPolynomial([-2, 0, 1]), 7)
    assert "simple" in guess
    fac = euler_factor_affine_line(5)
    assert fac(Fraction(1, 25)) == Fraction(1) / (1 - Fraction(5, 25))


def test_local_field():
    K = qp(5, 8)
    assert K.int(3).p == 5
    assert K.number(Fraction(1, 2)).is_unit()
    assert K.uniformizer().val == 1
    w = K.teichmuller(2)
    assert w.residue % 5 == 2
    L = K.unramified(3)
    assert L.extension.f == 3 and L.extension.e == 1
    R = K.totally_ramified(2)
    assert R.extension.e == 2
    ext = ExtensionData(p=5, n=6, e=2, f=3)
    assert ext.residue_cardinality() == 125
    assert "C_" in describe_unit_group(5, 3)


def test_cli_val_and_legendre(capsys):
    assert main(["val", "360", "5"]) == 0
    assert "1" in capsys.readouterr().out
    assert main(["legendre", "2", "7"]) == 0
    assert main(["expand", "42", "5", "--prec", "4"]) == 0
    assert main(["teich", "3", "7", "--prec", "4"]) == 0
    assert main(["info", "5"]) == 0
    assert main(["hilbert", "2", "5", "7"]) == 0
    assert main(["arith", "add", "3", "4", "5", "--prec", "4"]) == 0
    assert main(["units", "5", "2"]) == 0


def test_cli_hensel_sqrt(capsys):
    assert main(["hensel", "-2", "0", "1", "--p", "7", "--root", "3", "--prec", "5"]) == 0
    assert main(["sqrt", "2", "7", "--prec", "5"]) == 0
    assert main(["norm", "25", "5"]) == 0
