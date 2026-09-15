# From checking to certifying: Turing's method

This document explains the gap between what this repository does and what a
rigorous verification of RH up to a given height would require. It is a
tutorial and a roadmap. **None of it is implemented here yet.**

## The problem with the argument principle

`counting.count_zeros_argument_principle` gives a genuinely independent count
of zeros in the critical strip. It also gets expensive fast: the integrand
zeta'/zeta oscillates more as the contour gets taller, so the quadrature needs
more nodes, and each node is a zeta evaluation. At `maxdegree = 10`,
`dps = 35`, a single evaluation at T = 200 takes minutes; the largest height
this project has checked is T = 200.

Nobody verifies RH to height 10^13 by contour integration. Turing's method is
what makes large-scale verification possible, and it is worth understanding
precisely because it shows what a *certified* computation looks like.

## The setup

The counting function splits exactly:

    N(T) = theta(T)/pi + 1 + S(T)

where theta is the Riemann–Siegel theta function — smooth, explicit, cheap to
evaluate to any precision — and S(T) = (1/pi) arg zeta(1/2 + iT) is everything
else. This is an identity, not an approximation. It is also the source of the
circularity trap described in [`claim_status.md`](claim_status.md): you cannot
"check" N(T) against theta(T)/pi + 1 and learn anything, because their
difference is S(T) by definition.

Now observe three things.

1. **N(T) is an integer.**
2. **theta(T)/pi + 1 is computable to arbitrary precision.**
3. Therefore *determining N(T) exactly is equivalent to determining S(T) to
   within an interval of length less than 1* — you need only enough resolution
   to pick out which integer N(T) is.

So the whole problem reduces to bounding S(T).

## Why bounding S(T) pointwise doesn't work

Unconditionally, the best known pointwise bounds are of the shape
S(T) = O(log T), with explicit versions giving something like
|S(T)| ≤ 0.111 log T + 0.275 log log T + 2.450 for T sufficiently large
(Trudgian 2014, refining Backlund; **check the constants against the paper
before relying on them** — they have been improved several times).

At T = 10^13 that bound is around 5 — far too weak. It permits eleven different
integer values of N(T). Pointwise bounds on S are hopeless for this purpose.

## Turing's insight

Turing's 1953 observation is that while S(T) itself is badly behaved, its
**antiderivative** is not. Write

    S_1(T) = integral from 0 to T of S(t) dt.

S_1 is bounded far better than S, because S oscillates around zero and the
integral averages the oscillation away. Explicit bounds of the form

    | integral from t1 to t2 of S(t) dt |  <=  a + b log(t2)

are available, with small explicit constants (Lehman 1970; Brent 1979;
Trudgian 2011 gives a = 2.067, b = 0.059 for t1, t2 ≥ 168 pi — again, **verify
against the source**).

Combine that with the structure of S: between consecutive zeros on the line,
S decreases smoothly and jumps up by 1 at each zero. So if you have *found* a
certain number of sign changes of Z below T, you get a candidate value for
N(T); the integral bound then rules out every other integer, provided you
evaluate Z a little way past T and the Gram points there behave.

That is Turing's method: **you never bound S(T) directly; you bound its
integral over a window past your target height, and use the integer-valuedness
of N to snap to the answer.**

## The algorithm, concretely

To certify N(T) = m:

1. Find sign changes of Z on (0, T] by evaluating Z on a grid. Suppose you
   find m of them. Each is a genuine on-line zero, so N(T) ≥ m.
2. Continue evaluating Z past T, at Gram points g_n, g_{n+1}, ..., until you
   find a run where Gram's law holds (Z alternates in sign).
3. Apply the explicit bound on the integral of S over that run.
4. Conclude that S(T) lies in an interval containing exactly one admissible
   integer, hence N(T) = m exactly.
5. Since you already exhibited m zeros on the line, every zero below T is on
   the line — and each is simple, since m distinct sign changes account for
   all m zeros counted with multiplicity.

Step 5 is the payoff and it is worth restating: RH up to height T follows from
an *equality of two counts*, never from the on-line count alone.

## What "certified" additionally requires

Turing's method makes the counting tractable. Making it *rigorous* needs more:

- **Interval arithmetic** for every zeta and Z evaluation, so each computed
  sign is provably correct rather than probably correct. `arb` / `python-flint`
  is the standard tool; `mpmath` cannot do this.
- **Validated quadrature** with an enclosure, if any integral is evaluated.
- **Explicit truncation bounds** on the Riemann–Siegel expansion used to
  evaluate Z, not just the empirical accuracy mpmath targets.
- **Verified constants** in the S_1 bound, traced to a specific paper and
  matched to the hypotheses (the bounds have ranges of validity in t).

Only when all four hold does the output become a mathematical statement rather
than strong evidence.

## Roadmap for this repository

**Milestone 1 — done.** Reproduce and visualize critical-line zeros with
documented precision and tests, and an independent strip count at modest
heights. Delivered as an installable package with 118 tests, figure generation,
and the strip-vs-line comparison agreeing up to T = 200. Labelled Numerical
evidence throughout; the comparison routine is named
`experimentally_check_up_to`, not `certify_up_to`.

**Milestone 2 — next.** An educational Turing-method workflow: Gram points
(done, in `critical_line.gram`), Rosser blocks, Gram's law failures past
n = 126, and the S_1 bound applied with constants taken from Trudgian (2011)
and checked against the paper's stated range of validity. To be marked
experimental throughout, because `mpmath` arithmetic is not validated — the
workflow would show *how* certification works without producing a certificate.

**Milestone 3 — later.** Port the sign evaluations to interval arithmetic
(`python-flint` / `arb`) so that a run produces a genuine certificate for a
stated height, with every constant traced to a source. Only at this point could
any label in this project honestly change from Numerical evidence to
Certified.

And even then: for a specific finite height, and no further. The Riemann
Hypothesis would still be open. It is a statement about infinitely many zeros,
and no computation reaches them. See
[`mathematical_scope.md`](mathematical_scope.md).

## References

Full citations, including the numerical-verification records, are in
[`references.md`](references.md). The core sources for this document are
Turing (1953), Lehman (1970), Brent (1979), Trudgian (2011) for the explicit
`S_1` constants, and Platt (2017) for a modern rigorous-interval
implementation.
