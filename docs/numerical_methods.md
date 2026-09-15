# Numerical methods

How each computation in this project works, where its error comes from, and
why none of the results is certified.

## Arbitrary-precision arithmetic

Everything runs on [`mpmath`](https://mpmath.org/), which represents real
numbers as arbitrary-precision binary floats and lets the working precision be
set at runtime (`mp.mp.dps`, in decimal digits).

Precision is **global mutable state**, which is a footgun. A function that
sets it and forgets to restore it silently changes the accuracy of everything
computed afterwards. This package routes every precision change through
`zeta_tools.working_precision`, a context manager that restores the previous
value on exit, including on exception.

The same footgun has a subtler form that bit this project's own test suite.
Consider:

```python
g = gram_point(5, dps=30)                  # computed at 30 digits
assert abs(theta(g) - 5 * mp.pi) < 1e-20   # compared at the ambient precision
```

The comparison runs at whatever precision is current *outside* the call — 15
digits by default — so `mp.pi` is only accurate to 15 digits and the residual
cannot fall below about `1e-16` no matter how good the root is. The assertion
fails, and it looks like the root finder is broken when the arithmetic around
it is at fault. The fix is to do the comparison inside a matching context:

```python
g = gram_point(5, dps=30)
with mp.workdps(30):
    assert abs(theta(g) - 5 * mp.pi) < 1e-20   # residual is 0.0
```

Several tests in `tests/` carry this pattern with a comment, because it is the
kind of mistake that recurs.

## Grid-and-bisection zero searches

`zeros.find_on_line_zeros` samples Hardy's `Z(t)` on a uniform grid, brackets
each sign change, and bisects.

Hardy's `Z` is real for real `t` and satisfies `|Z(t)| = |zeta(1/2 + it)|`, so
its sign changes are exactly the odd-order zeros of zeta on the critical line.

Three failure modes, all of them real:

1. **A close pair can be stepped over.** If two zeros are separated by less
   than the grid step, `Z` changes sign twice inside one interval and the scan
   sees nothing. This is the Lehmer phenomenon, and it gets worse with height
   as the zeros crowd together. Default step `0.05` is comfortable at the
   heights this project reaches and would be inadequate near `t = 7005`, where
   Lehmer's famous pair sits.
2. **An even-order zero produces no sign change at all.** No such zero is
   known — all computed zeros are simple — but nothing here rules one out.
3. **Bisection terminates on width, not on a proven enclosure.** The returned
   ordinate is where the bracket collapsed, not a point proved to contain a
   zero. `zeros.zero_residual` reports `|zeta(1/2 + i*gamma)|` at the returned
   value so this is visible rather than assumed; it is small, never zero.

None of these is silent *if* an independent total count is available to
contradict an undercount. That is the argument for the next section.

## The argument principle

`counting.count_zeros_argument_principle` evaluates

    (1 / 2 pi i) * contour integral of zeta'(s)/zeta(s) ds

around the rectangle `Re(s) in [-1, 2]`, `Im(s) in [0.5, T]`.

By the argument principle this equals the number of zeros of zeta inside,
counted with multiplicity, minus the number of poles. The rectangle contains
the entire critical strip; the pole at `s = 1` and every trivial zero lie on
the real axis and are excluded by the lower edge at `Im(s) = 0.5`, which also
sits below the first zero ordinate `14.1347...`.

What makes this the *independent* count: the integrand mentions only `zeta`
and `zeta'`. No Riemann-Siegel theta, no Gram points, no Hardy `Z`, no Rosser
blocks. It also counts zeros **wherever they are in the strip**, which is the
property that gives the comparison against on-line zeros its content.

Three practical constraints:

**Cost grows with `T`.** The integrand oscillates faster along the taller
edges, so the quadrature needs more nodes, each one a zeta evaluation. This is
why the method is used to heights of order hundreds here rather than millions,
and why Turing's method exists.

**`T` must avoid zero ordinates.** The top edge of the contour runs through the
strip at height `T`. If a zero sits near it, `zeta` is nearly zero on the path,
the integrand is large, and accuracy collapses. Measured, at
`maxdegree = 8`, `dps = 30`:

| `T` | residual | note |
|---|---|---|
| 30 | 3.9e-31 | |
| 50 | 1.8e-27 | |
| 100 | 7.2e-17 | |
| 120 | 3.9e-11 | |
| 150 | 1.1e-08 | `Z(150) = -0.09`, nearly a zero |

Raising to `maxdegree = 10`, `dps = 35` recovers `T = 150` (residual 4.8e-35)
and reaches `T = 200` (residual 3.1e-27), at roughly 2.5-4 minutes per
evaluation.

**Quadrature parameters are not self-tuning.** `mp.quad` with a given
`maxdegree` returns an answer whether or not it has converged.

## Residuals and conditioning

In exact arithmetic the contour integral is an integer. The computed value is
not, and the distance to the nearest integer is the **residual**.

`count_zeros_argument_principle_int` rounds only when the residual is within an
explicit `residual_tolerance` (default `1e-8`), and raises
`ResidualToleranceError` otherwise. It does not round and warn. Rounding an
unconverged result manufactures a confident integer out of a failed
computation, and the failure then propagates into whatever conclusion depends
on it.

The asymmetry worth internalising:

- **A large residual is strong evidence the count is wrong.** It is a genuine
  detector of failure.
- **A small residual is a weak indication the count is right.** It says the
  quadrature is self-consistent at this refinement level. It does **not** bound
  the true error, because `mp.quad` returns no enclosure. A method can converge
  smoothly to the wrong answer.

So a small residual is a numerical diagnostic, not a theorem.

## Why `mpmath` is not interval-certified

`mpmath` computes with high precision, not with *guaranteed* precision. It
tracks how many digits it is carrying; it does not track how many of those
digits are correct. There is no directed rounding, no error propagation, no
enclosure returned alongside a value.

The practical consequence: when `hardy_z(t)` returns a small negative number,
that is mpmath's best estimate of `Z(t)`. If the true value is a small
*positive* number, the sign is wrong, a sign change is spurious or missed, and
the zero count is off — with nothing in the output to indicate it.

Interval arithmetic libraries (`arb`, via `python-flint`) return an enclosure
`[a, b]` guaranteed to contain the true value. If `0` is not in the enclosure,
the sign is *proved*. That is the difference between evidence and certificate,
and it is why no function in this package is named `certify`.

The comparison routine is `counting.experimentally_check_up_to`. The name is
the claim.

## Why a small numerical residual is evidence, not a theorem

Collecting the above into the rule the project runs on.

Every number produced here is the output of a finite-precision computation with
no rigorous error bound. The argument-principle strip count and the Hardy-`Z`
on-line count arrive at their answers by genuinely different routes — one
integrates `zeta'/zeta` around a contour, the other tracks sign changes of a
real function along a line — so a bug in one is unlikely to produce the same
wrong answer as a bug in the other. That is why their agreement is worth
something.

It is worth less than "independent", though, and the difference matters.
**These use different numerical counting procedures, but both rely on
finite-precision numerical evaluation and therefore provide experimental
evidence rather than a certified result.**

Concretely, the two share a foundation: every value either one needs comes out
of `mp.zeta`, at finite precision, with no enclosure. A systematic error in
mpmath's zeta evaluation — or simply a loss of precision severe enough to flip
a sign — could corrupt both counts at once, and their agreement would not
reveal it. The procedures differ; the arithmetic underneath them does not.

To turn evidence into certificate, four things are needed, and all four are
missing:

1. Interval arithmetic for every `zeta` and `Z` evaluation, so each sign is
   proved.
2. Validated quadrature returning an enclosure, so the contour count is proved.
3. Explicit truncation bounds on the Riemann-Siegel expansion used to evaluate
   `Z`.
4. Verified constants in any `S_1` bound, traced to a specific paper and
   matched against its stated range of validity.

The roadmap is in [`certification.md`](certification.md).
