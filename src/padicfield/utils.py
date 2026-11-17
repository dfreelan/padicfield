"""Number-theoretic utilities for padicfield (exact integer arithmetic)."""

from __future__ import annotations

import math
import random
from typing import Iterable, Sequence

from padicfield.exceptions import NotPrimeError


def is_prime(n: int) -> bool:
    """Return True if ``n`` is prime (deterministic Miller–Rabin for ordinary sizes)."""
    if n < 2:
        return False
    if n in (2, 3, 5, 7):
        return True
    if n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0:
        return False
    if n < 1_000_000:
        limit = int(math.isqrt(n)) + 1
        i = 11
        while i < limit:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True
    return _miller_rabin(n)


def require_prime(p: int, *, name: str = "p") -> int:
    """Validate that ``p`` is an int prime."""
    if not isinstance(p, int) or isinstance(p, bool):
        raise NotPrimeError(p, f"{name} must be an int prime, got {type(p).__name__}")
    if not is_prime(p):
        raise NotPrimeError(p, f"{name} must be prime, got {p}")
    return p


def _miller_rabin(n: int, *, rounds: int = 16) -> bool:
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    if n < 3_317_044_064_679_887_385_961_981:
        witnesses = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    else:
        witnesses = tuple(random.randrange(2, n - 1) for _ in range(rounds))
    for a in witnesses:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n: int) -> int:
    """Least prime strictly greater than ``n``."""
    if n < 2:
        return 2
    candidate = n + 1 + (n % 2 == 0)
    if candidate == 2:
        return 2
    if candidate % 2 == 0:
        candidate += 1
    while not is_prime(candidate):
        candidate += 2
    return candidate


def primes_upto(limit: int) -> list[int]:
    """All primes ``<= limit``."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            start = i * i
            sieve[start : limit + 1 : i] = b"\x00" * (((limit - start) // i) + 1)
    return [i for i in range(2, limit + 1) if sieve[i]]


def factorize(n: int) -> dict[int, int]:
    """Prime factorization of ``|n|`` as ``{p: e}``."""
    if n == 0:
        raise ValueError("cannot factorize 0")
    n = abs(n)
    if n == 1:
        return {}
    factors: dict[int, int] = {}
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    f = 3
    while f * f <= n:
        while n % f == 0:
            factors[f] = factors.get(f, 0) + 1
            n //= f
        f += 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def gcd(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def egcd(a: int, b: int) -> tuple[int, int, int]:
    """Return ``(g, x, y)`` with ``a x + b y = g = gcd(a, b)``."""
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    g, x1, y1 = egcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def modinv(a: int, m: int) -> int:
    """Inverse of ``a`` modulo ``m``."""
    if m <= 0:
        raise ValueError(f"modulus must be positive, got {m}")
    g, x, _ = egcd(a % m, m)
    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {m}")
    return x % m


def crt(residues: Sequence[int], moduli: Sequence[int]) -> int:
    """CRT for pairwise coprime moduli; result in ``[0, M)``."""
    if len(residues) != len(moduli) or not residues:
        raise ValueError("bad CRT system")
    M = 1
    for m in moduli:
        if m <= 0:
            raise ValueError(f"modulus must be positive, got {m}")
        M *= m
    result = 0
    for a_i, m_i in zip(residues, moduli):
        M_i = M // m_i
        result = (result + a_i * M_i * modinv(M_i, m_i)) % M
    return result


def legendre_symbol(a: int, p: int) -> int:
    """Legendre symbol ``(a/p) ∈ {-1, 0, 1}`` for odd prime ``p``."""
    require_prime(p)
    if p == 2:
        raise ValueError("Legendre symbol requires an odd prime")
    a %= p
    if a == 0:
        return 0
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def jacobi_symbol(a: int, n: int) -> int:
    """Jacobi symbol ``(a/n)`` for odd positive ``n``."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("Jacobi symbol requires odd positive n")
    a %= n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def kronecker_symbol(a: int, n: int) -> int:
    """Kronecker symbol ``(a/n)`` for ``n != 0``."""
    if n == 0:
        return 1 if abs(a) == 1 else 0
    result = 1
    if n < 0:
        if a < 0:
            result = -1
        n = -n
    v2 = 0
    while n % 2 == 0:
        n //= 2
        v2 += 1
    if v2:
        if a % 2 == 0:
            return 0
        two_sym = 1 if ((a * a - 1) // 8) % 2 == 0 else -1
        if v2 % 2 == 1:
            result *= two_sym
    if n == 1:
        return result
    return result * jacobi_symbol(a, n)


def binomial(n: int, k: int) -> int:
    """Binomial coefficient ``C(n, k)``."""
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def digits_base_p(n: int, p: int, length: int | None = None) -> list[int]:
    """Base-``p`` digits of ``n`` (LSB first), optionally padded/truncated."""
    if p < 2 or n < 0:
        raise ValueError("invalid base or negative n")
    digs: list[int] = []
    x = n
    if x == 0:
        digs = [0]
    else:
        while x:
            digs.append(x % p)
            x //= p
    if length is not None:
        if length < 0:
            raise ValueError("length must be non-negative")
        digs = (digs + [0] * length)[:length]
    return digs


def from_digits_base_p(digits: Iterable[int], p: int) -> int:
    """Integer from base-``p`` digits (LSB first)."""
    if p < 2:
        raise ValueError("base must be >= 2")
    result, power = 0, 1
    for d in digits:
        if d < 0 or d >= p:
            raise ValueError(f"digit {d} out of range for base {p}")
        result += d * power
        power *= p
    return result


def inv_mod_pk(a: int, p: int, k: int) -> int:
    """Inverse of a unit ``a`` modulo ``p**k`` via Hensel lifting."""
    require_prime(p)
    if k < 1:
        raise ValueError("k must be >= 1")
    mod = p**k
    a %= mod
    if gcd(a, p) != 1:
        raise ValueError(f"{a} is not a unit modulo {p}^{k}")
    inv = modinv(a % p, p)
    precision = 1
    while precision < k:
        precision = min(2 * precision, k)
        inv = (inv * (2 - a * inv)) % (p**precision)
    return inv % mod
