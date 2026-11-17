"""padicfield — p-adic arithmetic and local-field numerics."""

from __future__ import annotations

from padicfield.exceptions import (
    ConvergenceError,
    DomainError,
    HenselError,
    IncompatiblePrecisionError,
    IncompatiblePrimeError,
    NotPrimeError,
    PadicError,
    PolynomialError,
    PrecisionError,
    ResidueSymbolError,
    SquareRootError,
    TeichmullerError,
    ValuationError,
)
from padicfield.exp_log import exp_converges, log_converges, padic_exp, padic_log
from padicfield.hensel import find_roots_mod_p, hensel_lift, hensel_lift_simple, lift_all_simple_roots
from padicfield.local_field import ExtensionData, LocalField, qp
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.padic_number import PadicNumber, embed_rational
from padicfield.polynomials import IntPolynomial, PadicPolynomial
from padicfield.residue import hilbert_symbol, is_square_mod_pk, quadratic_residue_symbol
from padicfield.sqrt import is_square_in_qp, sqrt_padic_int, sqrt_padic_number, sqrt_rational
from padicfield.teichmuller import teichmuller_representative
from padicfield.utils import is_prime, jacobi_symbol, kronecker_symbol, legendre_symbol, require_prime
from padicfield.valuation import padic_norm, unit_part_int, valuation, valuation_int

__version__ = "0.1.0"


# Stable fixture catalog (internal golden IDs)
from ._catalog import FIXTURE_NOTES, REFERENCE_URLS, fixture_ids, note_for, primary_reference

__all__ = ["FIXTURE_NOTES", "REFERENCE_URLS", "fixture_ids", "note_for", "primary_reference", 
    "__version__",
    "PadicInt", "PadicNumber", "DEFAULT_PRECISION", "IntPolynomial", "PadicPolynomial",
    "LocalField", "ExtensionData", "qp",
    "valuation", "valuation_int", "padic_norm", "unit_part_int",
    "is_prime", "require_prime", "legendre_symbol", "jacobi_symbol", "kronecker_symbol",
    "quadratic_residue_symbol", "hilbert_symbol", "is_square_mod_pk",
    "hensel_lift", "hensel_lift_simple", "find_roots_mod_p", "lift_all_simple_roots",
    "teichmuller_representative",
    "sqrt_padic_int", "sqrt_padic_number", "sqrt_rational", "is_square_in_qp",
    "padic_exp", "padic_log", "exp_converges", "log_converges", "embed_rational",
    "PadicError", "NotPrimeError", "PrecisionError", "ValuationError", "ConvergenceError",
    "HenselError", "SquareRootError", "TeichmullerError", "PolynomialError",
    "ResidueSymbolError", "DomainError", "IncompatiblePrimeError", "IncompatiblePrecisionError",
]
