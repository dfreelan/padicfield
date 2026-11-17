"""Hensel's lemma and root lifting from ``F_p`` to ``Z_p``."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from padicfield.exceptions import HenselError, PrecisionError
from padicfield.padic_int import PadicInt
from padicfield.utils import inv_mod_pk, modinv, require_prime


@dataclass(frozen=True)
class HenselLiftResult:
    """Outcome of a Hensel lift."""

    root: PadicInt
    precision: int
    simple: bool

    @property
    def p(self) -> int:
        return self.root.p


def eval_poly(coeffs: Sequence[int], x: int, mod: int | None = None) -> int:
    """Horner evaluation of ``Σ c_i x^i``."""
    acc = 0
    for c in reversed(coeffs):
        acc = acc * x + c
        if mod is not None:
            acc %= mod
    return acc


def eval_poly_derivative(coeffs: Sequence[int], x: int, mod: int | None = None) -> int:
    """Evaluate formal derivative at ``x``."""
    if len(coeffs) <= 1:
        return 0
    return eval_poly([i * coeffs[i] for i in range(1, len(coeffs))], x, mod)


def poly_derivative_coeffs(coeffs: Sequence[int]) -> list[int]:
    if len(coeffs) <= 1:
        return [0]
    return [i * coeffs[i] for i in range(1, len(coeffs))]


def find_roots_mod_p(coeffs: Sequence[int], p: int) -> list[int]:
    """Brute-force roots in ``F_p`` (small primes only)."""
    require_prime(p)
    if p > 10_000:
        raise HenselError(f"brute-force root search unsuitable for p={p}")
    return [a for a in range(p) if eval_poly(coeffs, a, p) == 0]


def is_simple_root(coeffs: Sequence[int], a: int, p: int) -> bool:
    require_prime(p)
    return eval_poly(coeffs, a, p) == 0 and eval_poly_derivative(coeffs, a, p) % p != 0


def hensel_lift_simple(coeffs: Sequence[int], a0: int, p: int, prec: int) -> HenselLiftResult:
    """Newton-lift a simple root from ``mod p`` to ``mod p^{prec}``."""
    require_prime(p)
    if prec < 1:
        raise PrecisionError("precision must be >= 1", required=1, available=prec)
    a0 %= p
    if eval_poly(coeffs, a0, p) != 0:
        raise HenselError(f"{a0} is not a root of f modulo {p}")
    if eval_poly_derivative(coeffs, a0, p) % p == 0:
        raise HenselError(f"derivative vanishes at {a0} mod {p}; root is not simple")
    a, current = a0, 1
    while current < prec:
        next_prec = min(2 * current, prec)
        mod = p**next_prec
        fa = eval_poly(coeffs, a, mod)
        fpa = eval_poly_derivative(coeffs, a, mod)
        try:
            inv_fp = inv_mod_pk(fpa, p, next_prec)
        except ValueError as exc:
            raise HenselError(f"derivative became non-unit during lifting: {exc}") from exc
        a = (a - fa * inv_fp) % mod
        current = next_prec
    return HenselLiftResult(root=PadicInt(a, p, prec), precision=prec, simple=True)


def hensel_lift_step(coeffs: Sequence[int], a: int, p: int, k: int) -> int:
    """Lift a simple root from ``mod p^k`` to ``mod p^{k+1}``."""
    require_prime(p)
    if k < 1:
        raise ValueError("k must be >= 1")
    mod_k, mod_next = p**k, p ** (k + 1)
    a %= mod_k
    fa = eval_poly(coeffs, a, mod_next)
    if fa % mod_k != 0:
        raise HenselError(f"f(a) is not 0 mod p^{k}")
    t = (fa // mod_k) % p
    fpa = eval_poly_derivative(coeffs, a, p) % p
    if fpa == 0:
        raise HenselError("derivative vanishes mod p; cannot take a simple Hensel step")
    b = (-t * modinv(fpa, p)) % p
    return (a + b * mod_k) % mod_next


def hensel_lift(
    coeffs: Sequence[int],
    a0: int,
    p: int,
    prec: int,
    *,
    allow_multiple: bool = False,
) -> HenselLiftResult:
    """Lift a root from ``mod p`` to ``mod p^{prec}``."""
    require_prime(p)
    a0 %= p
    if eval_poly(coeffs, a0, p) != 0:
        raise HenselError(f"{a0} is not a root modulo {p}")
    if is_simple_root(coeffs, a0, p):
        return hensel_lift_simple(coeffs, a0, p, prec)
    if not allow_multiple:
        raise HenselError(
            f"root {a0} is not simple mod {p}; pass allow_multiple=True to attempt generalized lifting"
        )
    return _hensel_lift_multiple(coeffs, a0, p, prec)


def _hensel_lift_multiple(coeffs: Sequence[int], a0: int, p: int, prec: int) -> HenselLiftResult:
    """Best-effort generalized Hensel lifting for nonsimple roots."""
    a = a0 % p
    for k in range(1, prec):
        mod_k, mod_next = p**k, p ** (k + 1)
        fa = eval_poly(coeffs, a, mod_next)
        if fa % mod_k != 0:
            raise HenselError(f"lost root property at precision {k}")
        fpa = eval_poly_derivative(coeffs, a, mod_next)
        m = 0 if fa == 0 else 0
        tmp = fa
        if fa == 0:
            m = k + 1
        else:
            while tmp % p == 0:
                tmp //= p
                m += 1
        n = 0
        tmp = fpa
        if fpa == 0:
            n = k + 1
        else:
            while tmp % p == 0:
                tmp //= p
                n += 1
        if m < 2 * n:
            raise HenselError(f"no lift exists at step {k}: m={m}, n={n}")
        if m > k:
            a %= mod_next
            continue
        found = None
        for b in range(p):
            candidate = (a + b * mod_k) % mod_next
            if eval_poly(coeffs, candidate, mod_next) % mod_next == 0:
                found = candidate
                break
        if found is None and fpa % p != 0:
            t = (fa // mod_k) % p
            bb = (-t * modinv(fpa % p, p)) % p
            found = (a + bb * mod_k) % mod_next
        if found is None:
            raise HenselError(f"failed to lift multiple root at precision {k}")
        a = found
    return HenselLiftResult(root=PadicInt(a, p, prec), precision=prec, simple=False)


def lift_all_simple_roots(coeffs: Sequence[int], p: int, prec: int) -> list[HenselLiftResult]:
    return [
        hensel_lift_simple(coeffs, a0, p, prec)
        for a0 in find_roots_mod_p(coeffs, p)
        if is_simple_root(coeffs, a0, p)
    ]


def newton_lift(
    f: Callable[[int, int], int],
    fprime: Callable[[int, int], int],
    a0: int,
    p: int,
    prec: int,
) -> PadicInt:
    """Abstract Newton lift for callables ``f`` / ``fprime``."""
    require_prime(p)
    a = a0 % p
    if f(a, p) % p != 0:
        raise HenselError("initial value is not a root mod p")
    if fprime(a, p) % p == 0:
        raise HenselError("derivative vanishes; Newton lift requires a simple root")
    current = 1
    while current < prec:
        next_prec = min(2 * current, prec)
        mod = p**next_prec
        a = (a - f(a, mod) * inv_mod_pk(fprime(a, mod), p, next_prec)) % mod
        current = next_prec
    return PadicInt(a, p, prec)


def verify_root(coeffs: Sequence[int], root: PadicInt) -> bool:
    return eval_poly(coeffs, root.residue, root.modulus) % root.modulus == 0


def multiplicity_mod_p(coeffs: Sequence[int], a: int, p: int) -> int:
    """Multiplicity of root ``a`` in ``F_p[x]`` via successive derivatives."""
    require_prime(p)
    a %= p
    c = list(coeffs)
    m = 0
    while eval_poly(c, a, p) % p == 0:
        m += 1
        c = poly_derivative_coeffs(c)
        if c == [0] or not any(x % p for x in c):
            return m
    return m
