"""Teichmüller representatives in ``Z_p``."""

from __future__ import annotations

from padicfield.exceptions import PrecisionError, TeichmullerError
from padicfield.padic_int import DEFAULT_PRECISION, PadicInt
from padicfield.utils import factorize, inv_mod_pk, require_prime


def teichmuller_representative(a: int, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
    """Return the Teichmüller lift ``ω(a)`` of ``a mod p`` to precision ``prec``.

    For ``a ≡ 0`` returns zero. For odd ``p``, Newton-lifts ``x^{p-1} - 1 = 0``.
    For ``p = 2``, returns the idempotent ``0`` or ``1``.
    """
    require_prime(p)
    if prec < 1:
        raise PrecisionError("precision must be >= 1", required=1, available=prec)
    a %= p
    if a == 0:
        return PadicInt.zero(p, prec)
    if p == 2:
        return PadicInt.one(2, prec)
    return _teichmuller_newton(a, p, prec)


def _teichmuller_newton(a: int, p: int, prec: int) -> PadicInt:
    x = a % p
    if pow(x, p - 1, p) != 1:
        raise TeichmullerError(f"{a} is not a unit residue mod {p}")
    current = 1
    while current < prec:
        next_prec = min(2 * current, prec)
        mod = p**next_prec
        xp2 = pow(x, p - 2, mod)
        fx = (xp2 * x - 1) % mod
        fpx = ((p - 1) * xp2) % mod
        try:
            inv = inv_mod_pk(fpx, p, next_prec)
        except ValueError as exc:
            raise TeichmullerError(f"Newton step failed: {exc}") from exc
        x = (x - fx * inv) % mod
        current = next_prec
    return PadicInt(x, p, prec)


def teichmuller_character_value(a: int, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
    return teichmuller_representative(a, p, prec)


def teichmuller_power(a: int, n: int, p: int, prec: int = DEFAULT_PRECISION) -> PadicInt:
    """``ω(a)^n = ω(a^n)`` via the character property on residues."""
    require_prime(p)
    if a % p == 0:
        if n > 0:
            return PadicInt.zero(p, prec)
        if n == 0:
            return PadicInt.one(p, prec)
        raise TeichmullerError("negative power of ω(0) is undefined")
    if p == 2:
        return PadicInt.one(2, prec) if a % 2 else PadicInt.zero(2, prec)
    target = pow(a % p, n % (p - 1), p)
    return teichmuller_representative(target, p, prec)


def all_teichmuller(p: int, prec: int = DEFAULT_PRECISION) -> list[PadicInt]:
    require_prime(p)
    if p > 500:
        raise TeichmullerError(f"refusing to enumerate Teichmüller set for large p={p}")
    return [teichmuller_representative(a, p, prec) for a in range(p)]


def teichmuller_log_table(p: int, prec: int = DEFAULT_PRECISION) -> dict[int, PadicInt]:
    return {a: teichmuller_representative(a, p, prec) for a in range(p)}


def is_teichmuller(x: PadicInt) -> bool:
    return x == teichmuller_representative(x.residue % x.p, x.p, x.prec)


def multiplicative_order_teichmuller(a: int, p: int) -> int:
    """Order of ``a`` in ``F_p^×`` (equals the order of ``ω(a)``)."""
    require_prime(p)
    a %= p
    if a == 0:
        raise TeichmullerError("ω(0) is not a unit")
    if p == 2:
        return 1
    n = p - 1
    for q, e in factorize(n).items():
        for _ in range(e):
            if pow(a, n // q, p) == 1:
                n //= q
            else:
                break
    return n


def decompose_unit(u: PadicInt) -> tuple[PadicInt, PadicInt]:
    """Write unit ``u = ω(ū) · ⟨u⟩`` with ``⟨u⟩ ∈ 1 + p Z_p``."""
    if not u.is_unit():
        raise TeichmullerError("decompose_unit requires a unit")
    omega = teichmuller_representative(u.residue % u.p, u.p, u.prec)
    return omega, u.divide(omega)


def frobenius_on_teichmuller(x: PadicInt) -> PadicInt:
    """Teichmüller representatives are fixed by Frobenius; verify and return."""
    if not is_teichmuller(x):
        raise TeichmullerError("expected a Teichmüller representative")
    return x


def teichmuller_from_idempotent(e: PadicInt) -> bool:
    """True if ``e^2 = e``."""
    return e * e == e
