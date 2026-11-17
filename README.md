# padicfield

**p-adic arithmetic and local-field numerics in pure Python.**

`padicfield` provides exact valuations, finite-precision elements of the rings
of p-adic integers \(\mathbf{Z}_p\) and the fields \(\mathbf{Q}_p\), together
with classical constructions from local number theory: Hensel lifting,
Teichmüller representatives, square roots, the p-adic exponential and
logarithm on their disks of convergence, quadratic residue and Hilbert
symbols, and simple local zeta helpers.

The runtime depends only on the Python standard library (exact integers and
`fractions.Fraction`). Precision is always finite and explicit: every
computation reports the modulus \(p^N\) or relative precision at which the
result is known.

## Installation

```bash
pip install -e ".[dev]"
```

Requires Python 3.10+.

## Quick start

```python
from padicfield import (
    PadicInt,
    PadicNumber,
    valuation,
    teichmuller_representative,
    hensel_lift_simple,
    sqrt_rational,
    padic_exp,
    padic_log,
    hilbert_symbol,
    LocalField,
)

# Valuation and norm
assert valuation(360, 5) == 1          # 360 = 5 * 72

# Z_p arithmetic at precision 8
a = PadicInt(17, p=5, prec=8)
b = PadicInt(3, p=5, prec=8)
print(a * b)                           # product mod 5^8
print(a.digits())                      # base-5 expansion

# Embed rationals into Q_p
x = PadicNumber.from_rational(2 / 7, p=5, prec=10)

# Teichmüller representative ω(3) in Z_7
w = teichmuller_representative(3, p=7, prec=6)
assert (w ** 6).congruent(1, 6)        # ω(a)^{p-1} = 1

# Hensel: lift the root x ≡ 3 (mod 7) of x^2 - 2
root = hensel_lift_simple([-2, 0, 1], a0=3, p=7, prec=8).root
assert (root * root).congruent(2, 8)

# Square root in Q_7
s = sqrt_rational(2, p=7, prec=8)

# Exp / log on principal units
K = LocalField(5, prec=12)
u = K.int(1 + 5 * 3)                   # ≡ 1 (mod 5)
y = padic_log(u)
restored = padic_exp(y)
assert restored.congruent(PadicNumber.from_padic_int(u), abs_prec=8)

# Hilbert symbol
assert hilbert_symbol(2, 5, 7) in (-1, 1)
```

## Command-line interface

The package installs a `padicfield` console script:

```bash
padicfield val 360 5
padicfield expand 42 5 --prec 8
padicfield teich 3 7 --prec 6
padicfield sqrt 2 7 --prec 8
padicfield hensel -2 0 1 --p 7 --root 3 --prec 6
padicfield legendre 10 13
padicfield hilbert 2 5 7
padicfield info 5
```

## Mathematical scope

| Topic | Module | Notes |
| --- | --- | --- |
| Valuation \(v_p\), \(\lvert\,\cdot\,\rvert_p\) | `valuation` | Exact on \(\mathbf{Q}\) |
| \(\mathbf{Z}_p\) at precision \(N\) | `padic_int` | Digits, units, inversion |
| \(\mathbf{Q}_p\) | `padic_number` | \(p^v\cdot u\) normal form |
| Hensel's lemma | `hensel` | Simple roots; Newton doubling |
| Teichmüller lifts | `teichmuller` | \(\omega\colon\mathbf{F}_p^\times\to\mathbf{Z}_p^\times\) |
| Square roots | `sqrt` | Tonelli–Shanks + Hensel; \(p=2\) via mod 8 |
| \(\exp\), \(\log\) | `exp_log` | Disks \(v_p(x)>1/(p-1)\) and principal units |
| Residue / Hilbert symbols | `residue` | Legendre, Kronecker, \((a,b)_p\) |
| Polynomials | `polynomials` | \(\mathbf{Z}[x]\), \((\mathbf{Z}/p^N)[x]\), resultants |
| Local zeta helpers | `zeta` | Point counts mod \(p^k\), Euler factors |
| Local-field context | `local_field` | \(e,f,n\) bookkeeping |

## Precision conventions

- A `PadicInt(value, p, prec=N)` represents a class in \(\mathbf{Z}/p^N\mathbf{Z}\).
- A `PadicNumber` stores \(p^{v}\cdot u\) with relative precision equal to the
  precision of the unit \(u\). Absolute precision is \(v + \mathrm{prec}(u)\).
- Arithmetic on mismatched precisions truncates to the minimum.

## Development

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/pytest -q
```

## License

MIT — see [LICENSE](LICENSE).

## References

- J.-P. Serre, *A Course in Arithmetic*
- N. Koblitz, *p-adic Numbers, p-adic Analysis, and Zeta-Functions*
- F. Q. Gouvêa, *p-adic Numbers: An Introduction*
