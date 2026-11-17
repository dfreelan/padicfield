"""Exception hierarchy for padicfield.

All library-specific errors inherit from :class:`PadicError` so callers can
catch a single base type when desired.
"""

from __future__ import annotations


class PadicError(Exception):
    """Base exception for all padicfield errors."""


class NotPrimeError(PadicError, ValueError):
    """Raised when a prime modulus or prime base is required but not provided."""

    def __init__(self, value: object, message: str | None = None) -> None:
        self.value = value
        msg = message or f"expected a prime, got {value!r}"
        super().__init__(msg)


class PrecisionError(PadicError, ValueError):
    """Raised when arithmetic would lose required p-adic precision."""

    def __init__(self, message: str, *, required: int | None = None, available: int | None = None) -> None:
        self.required = required
        self.available = available
        super().__init__(message)


class ValuationError(PadicError, ValueError):
    """Raised for invalid valuation operations (e.g. valuation of zero)."""


class ConvergenceError(PadicError, ValueError):
    """Raised when a p-adic series or iteration fails to converge."""

    def __init__(self, message: str, *, radius: object | None = None) -> None:
        self.radius = radius
        super().__init__(message)


class HenselError(PadicError, ValueError):
    """Raised when Hensel lifting cannot proceed (singular derivative, etc.)."""


class SquareRootError(PadicError, ValueError):
    """Raised when a p-adic square root does not exist."""


class TeichmullerError(PadicError, ValueError):
    """Raised when a Teichmüller representative cannot be constructed."""


class PolynomialError(PadicError, ValueError):
    """Raised for invalid polynomial operations over local rings."""


class ResidueSymbolError(PadicError, ValueError):
    """Raised when a quadratic residue symbol is undefined."""


class DomainError(PadicError, ValueError):
    """Raised when an argument lies outside the domain of a p-adic function."""


class IncompatiblePrimeError(PadicError, ValueError):
    """Raised when two p-adic objects use different primes."""

    def __init__(self, p1: int, p2: int) -> None:
        self.p1 = p1
        self.p2 = p2
        super().__init__(f"incompatible primes: {p1} vs {p2}")


class IncompatiblePrecisionError(PadicError, ValueError):
    """Raised when precision policies conflict and cannot be reconciled."""

    def __init__(self, prec1: int, prec2: int) -> None:
        self.prec1 = prec1
        self.prec2 = prec2
        super().__init__(f"incompatible precisions: {prec1} vs {prec2}")
