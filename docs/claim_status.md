# Claim status

Every claim in this project carries exactly one of four labels. The point of
the labels is to make it impossible to quote a number from `verify.py` as
though it were a theorem.

| Label | Meaning |
|---|---|
| **Proved here** | A full derivation appears in the notebook. You can check it line by line without leaving the document. |
| **Classical theorem** | A published result, stated precisely and attributed. Not re-derived here. |
| **Numerical evidence** | Computed. True at the precision used, for the range tested, and nowhere else. Not a proof. |
| **Open** | No claim is made. |

There is a fifth category this project does **not** currently occupy:

| **Certified** | Computed with rigorous error control (interval arithmetic, validated quadrature, explicit bounds on every truncation) so that the output is a machine-checked mathematical statement. |

Nothing in this repository is Certified. `mpmath` computes at high precision
but does not carry directed rounding or error bounds, so every number below is
Numerical evidence, not Certified. See [`certification.md`](certification.md)
for what closing that gap would take, and
[`numerical_methods.md`](numerical_methods.md) for where the error in each
method comes from.

This is also why no function in the package is named `certify`. The comparison
routine is `counting.experimentally_check_up_to`; the name states the claim.

---

## Notebook stages

| Stage | Content | Status |
|---|---|---|
| 0 | Notation, assumed background | — |
| 1 | Dirichlet series, convergence for Re(s) > 1 | Proved here |
| 2 | Euler product; divergence of the sum of 1/p | Proved here |
| 3 | Analytic continuation; the pole at s = 1 | Proved here |
| 4 | Theta function and Poisson summation | Proved here, using Poisson summation as **Classical theorem** |
| 5 | The functional equation; xi(s) = xi(1-s) | Proved here, using Gamma reflection/duplication and Stirling as **Classical theorem** |
| 6 | Riemann–von Mangoldt formula for N(T) | Proved here (main term); error term **Classical theorem** |
| 7 | Hardy's theorem: infinitely many zeros on the line; Selberg's positive proportion | **Classical theorem** (Hardy 1914; Selberg 1942; Conrey 1989 for 2/5) |
| 8 | The explicit formula for psi(x) | Proved here, using Perron's formula as **Classical theorem** |
| 9 | Equivalent formulations: Mertens (false), Li's criterion, Nyman–Beurling, de Bruijn–Newman | **Classical theorem** (Odlyzko–te Riele 1985; Bombieri–Lagarias 1999; Nyman 1950/Beurling 1955; Rodgers–Tao 2020) |
| 10 | Hilbert–Pólya, GUE statistics, Berry–Keating | **Classical theorem** where stated as theorems; **Open** where stated as conjecture (Montgomery's pair correlation is conjectural) |
| 11 | Ledger of what was and was not established | — |
| — | **The Riemann Hypothesis itself** | **Open** |

## Numerical checks in `verify.py`

| Check | What it computes | Status |
|---|---|---|
| Functional equation | \|zeta(s) − 2^s pi^(s−1) sin(pi s/2) Gamma(1−s) zeta(1−s)\| ~ 1e−26 | Numerical evidence |
| xi symmetry | \|xi(s) − xi(1−s)\| ~ 1e−26 | Numerical evidence |
| `mp.nzeros` circularity | nzeros(T) − theta(T)/pi − 1 − S(T) = 0 to working precision (7.3e−29 at 30 digits) | Numerical evidence (see correction below) |
| Strip count vs line count | argument principle vs Z sign changes, T ≤ 200 | Numerical evidence |
| Zero finder cross-check | grid+bisection vs `mp.zetazero`, agreement to 1e−22 | Numerical evidence |
| Gram's law | no failures for n < 40; first failure at n = 126 | Numerical evidence |
| Spacing statistics | unfolded gaps, GUE and Poisson overlay | Numerical evidence |
| Explicit formula | partial sums converge toward psi(100) = 94.045311 | Numerical evidence |
| Montgomery pair correlation | — | **Open** (conjecture; would not imply RH) |

---

## A correction to an earlier version of this project

An earlier `verify.py` compared `mp.nzeros(T)` against the Riemann–von Mangoldt
main term and described the agreement as an *independent* confirmation, with
`mp.nzeros` characterised as "an independent count of zeros directly on the
critical line." It also reported a "200 = 200" match between `mp.nzeros(T)` and
the number of ordinates returned by `mp.zetazero(1..200)`.

Both framings were wrong, for the same reason.

Reading `mpmath/functions/zetazeros.py`: `nzeros(t)` is built from
`gram_index(t) = floor(theta(t)/pi)` plus Rosser-block bookkeeping of `Z(t)`
sign changes, and `zetazero(n)` calls the same `find_rosser_block_zero` and
`separate_zeros_in_block` routines. So:

- `nzeros` vs the counting formula is not two computations agreeing. Since
  N(T) = theta(T)/pi + 1 + S(T) *by definition of S*, and `nzeros` is built
  from theta, the "discrepancy" being reported was simply S(T). The check
  reproduces an identity to whatever working precision is requested — measured
  worst case 1.7e−24 at 25 digits, 7.3e−29 at 30, 6.6e−49 at 50 — and carries
  no information about the location of any zero. (This is now asserted as a test, in
  `tests/test_counting.py::test_nzeros_is_not_independent_of_the_formula`.)
- `nzeros` vs counting `zetazero` outputs is circular: same code path.

The replacement is `counting.count_zeros_argument_principle`, which evaluates

    (1 / 2 pi i) * contour integral of zeta'(s)/zeta(s) ds

around the rectangle Re(s) in [−1, 2], Im(s) in [0.5, T]. That rectangle
contains the whole critical strip; the pole at s = 1 and all trivial zeros sit
on the real axis and are excluded. The integrand mentions only zeta and zeta'.
No Gram points, no theta, no Z, no Rosser blocks.

Comparing that count against the number of sign changes of Z on (0, T] is then
a real comparison: the first counts zeros *anywhere in the strip*, the second
finds zeros *on the critical line*. Equality means every zero below height T is
on the line and simple. That is what "RH verified up to height T" means, and it
is the check the project now runs.

Two caveats survive even so, and they are why the label is Numerical evidence
rather than Certified:

1. `mp.quad` returns no rigorous error bound. The reported quadrature residual
   (distance from the nearest integer) is a diagnostic, not a proof of
   convergence. It degrades as T grows, and the degradation is not monotone in
   T alone — a contour whose top edge passes near a zero loses accuracy
   sharply. At `maxdegree = 8`, `dps = 30`: 3.9e−31 at T = 30, 7.2e−17 at
   T = 100, 3.9e−11 at T = 120, and 1.1e−8 at T = 150, where `Z(150) = −0.09`
   puts a zero almost exactly on the contour. Raising to `maxdegree = 10`,
   `dps = 35` recovers T = 150 (4.8e−35) and reaches T = 200 (3.1e−27).

   Because of this, `count_zeros_argument_principle_int` refuses to round when
   the residual exceeds an explicit tolerance, raising `ResidualToleranceError`
   instead of returning a confidently wrong integer.
2. The sign-change grid could in principle step over a close pair of zeros.
   A coarse grid would show *fewer* line zeros than strip zeros, so this
   failure mode is detectable rather than silent — but only because the
   independent strip count exists to detect it.

The published record (van de Lune–te Riele–Winter, Gourdon, and the
distributed `ZetaGrid` and successor efforts) has checked the first 10^13 zeros
and beyond. That is far past anything here, and it is still not a proof.
