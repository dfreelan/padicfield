"""p-adic exponential and logarithm on their disks of convergence."""

from __future__ import annotations

import math
from math import factorial

from padicfield.exceptions import ConvergenceError, DomainError, PrecisionError
from padicfield.padic_int import PadicInt
from padicfield.padic_number import PadicNumber
from padicfield.utils import require_prime
from padicfield.valuation import valuation_of_factorial


def exp_convergence_radius_valuation(p: int) -> float:
    """Series converges for ``v_p(x) > 1/(p-1)``."""
    require_prime(p)
    return 1.0 / (p - 1)


def min_integral_valuation_for_exp(p: int) -> int:
    require_prime(p)
    return 2 if p == 2 else 1


def log_converges(x: PadicNumber | PadicInt) -> bool:
    """True if ``x`` is a principal unit (``1+pZ_p``, or ``1+4Z_2``)."""
    if isinstance(x, PadicInt):
        p = x.p
        if p == 2:
            return x.residue % 2 == 1 if x.prec < 2 else x.residue % 4 == 1
        return x.residue % p == 1
    p = x.p
    diff = x - PadicNumber.one(p, x.prec)
    if diff.is_zero():
        return True
    assert diff.val is not None
    return diff.val >= (2 if p == 2 else 1)


def exp_converges(x: PadicNumber) -> bool:
    return x.is_zero() or x.valuation > exp_convergence_radius_valuation(x.p)


def _v_p_of_int(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def _series_terms_needed(p: int, v_x: float, target_prec: int) -> int:
    if math.isinf(v_x):
        return 1
    threshold = exp_convergence_radius_valuation(p)
    if v_x <= threshold:
        raise ConvergenceError(
            f"v_p(x)={v_x} does not exceed 1/(p-1)={threshold}",
            radius=threshold,
        )
    N = 1
    while N < 10_000:
        if N * v_x - valuation_of_factorial(N, p) >= target_prec:
            if all(
                m * v_x - valuation_of_factorial(m, p) >= target_prec
                for m in range(N, N + p + 2)
            ):
                return N + p + 2
        N += 1
    raise ConvergenceError("unable to bound exp/log series truncation")


def _log_terms_needed(p: int, v_x: float, target_prec: int) -> int:
    if v_x <= 0:
        raise ConvergenceError("log series requires v_p(x) > 0", radius=0)
    N = 1
    while N < 10_000:
        if N * v_x - _v_p_of_int(N, p) >= target_prec:
            if all(m * v_x - _v_p_of_int(m, p) >= target_prec for m in range(N, N + p + 2)):
                return N + p + 2
        N += 1
    raise ConvergenceError("unable to bound log series truncation")


def padic_exp(x: PadicNumber, *, prec: int | None = None) -> PadicNumber:
    """``exp(x)`` via truncated series; requires ``v_p(x) > 1/(p-1)``."""
    p = x.p
    working_prec = prec if prec is not None else x.prec
    if working_prec < 1:
        raise PrecisionError("precision must be >= 1")
    if x.is_zero():
        return PadicNumber.one(p, working_prec)
    if not exp_converges(x):
        raise ConvergenceError(
            f"exp does not converge at valuation {x.valuation}",
            radius=exp_convergence_radius_valuation(p),
        )
    assert x.val is not None
    n_terms = _series_terms_needed(p, x.valuation, working_prec + max(0, x.val))
    acc = term = PadicNumber.one(p, working_prec)
    for n in range(1, n_terms + 1):
        term = term * x / n
        acc = acc + term
        if term.is_zero():
            break
    return acc.with_precision(working_prec)


def padic_log(u: PadicNumber | PadicInt, *, prec: int | None = None) -> PadicNumber:
    """``log(u)`` for a principal unit ``u``."""
    if isinstance(u, PadicInt):
        p, working = u.p, (prec if prec is not None else u.prec)
        u_num = PadicNumber.from_padic_int(u.with_precision(working))
    else:
        p, working = u.p, (prec if prec is not None else u.prec)
        u_num = u.with_precision(working)
    if not log_converges(u_num):
        raise DomainError(f"log requires a principal unit; got {u_num}")
    x = u_num - PadicNumber.one(p, working)
    if x.is_zero():
        return PadicNumber.zero(p, working)
    assert x.val is not None
    n_terms = _log_terms_needed(p, x.valuation, working + max(0, -x.val))
    acc = PadicNumber.zero(p, working)
    xpow = PadicNumber.one(p, working)
    for n in range(1, n_terms + 1):
        xpow = xpow * x
        term = xpow / n
        acc = acc - term if n % 2 == 0 else acc + term
        if term.is_zero():
            break
    return acc.with_precision(working)


def padic_log_int(u: PadicInt, *, prec: int | None = None) -> PadicNumber:
    return padic_log(u, prec=prec)


def padic_exp_int(x: PadicInt, *, prec: int | None = None) -> PadicNumber:
    return padic_exp(PadicNumber.from_padic_int(x), prec=prec)


def exp_log_roundtrip_ok(x: PadicNumber, *, tol_val: int | None = None) -> bool:
    if not exp_converges(x):
        return False
    e = padic_exp(x)
    if not log_converges(e):
        return False
    diff = padic_log(e) - x
    if diff.is_zero():
        return True
    assert diff.val is not None
    return diff.val >= (tol_val if tol_val is not None else max(1, x.prec - 2))


def formal_exp_series(p: int, max_n: int) -> list[tuple[int, int]]:
    require_prime(p)
    return [(1, factorial(n)) for n in range(max_n + 1)]


def formal_log_series(max_n: int) -> list[tuple[int, int]]:
    coeffs: list[tuple[int, int]] = [(0, 1)]
    for n in range(1, max_n + 1):
        coeffs.append((1 if n % 2 == 1 else -1, n))
    return coeffs


def padic_pow_principal(u: PadicInt, alpha: PadicNumber) -> PadicNumber:
    if not log_converges(u):
        raise DomainError("base must be a principal unit")
    return padic_exp(alpha * padic_log(u))


def binomial_series_unit(
    x: PadicNumber, alpha: int | PadicNumber, *, prec: int | None = None
) -> PadicNumber:
    """``(1+x)^α`` via direct power or exp-log when convergent."""
    p = x.p
    working = prec if prec is not None else x.prec
    if x.is_zero():
        return PadicNumber.one(p, working)
    base = PadicNumber.one(p, working) + x
    if isinstance(alpha, int) and 0 <= alpha < 50:
        return base**alpha
    if not log_converges(base):
        raise DomainError("binomial series domain not satisfied")
    a = alpha if isinstance(alpha, PadicNumber) else PadicNumber.from_int(alpha, p, working)
    return padic_exp(a * padic_log(base), prec=working)
