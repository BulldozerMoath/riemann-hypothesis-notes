# Mathematical scope

What this project is about, what it actually computes, and the precise sense
in which those are different things.

## What the Riemann Hypothesis says

The Riemann zeta function is defined for `Re(s) > 1` by

    zeta(s) = sum over n >= 1 of n^(-s)

and extended to the rest of the complex plane by analytic continuation. The
continued function has a single simple pole at `s = 1` and vanishes at the
negative even integers `-2, -4, -6, ...`; those are the *trivial* zeros, and
they come from the poles of the Gamma factor in the functional equation.

Every other zero — the *nontrivial* zeros — lies in the critical strip
`0 < Re(s) < 1`. There are infinitely many, and they are symmetric about both
the real axis and the line `Re(s) = 1/2`.

**The Riemann Hypothesis:** every nontrivial zero has `Re(s) = 1/2`.

That is the entire statement. It is open. It has been open since 1859, it is
one of the seven Clay Millennium Prize Problems, and nothing in this
repository moves it.

Its importance is that the zeros control the distribution of primes. Von
Mangoldt's explicit formula makes this exact: the prime-counting function is a
smooth main term minus an oscillating sum over the zeros, with each zero
contributing a wave whose amplitude is governed by its real part. Zeros on the
critical line give the smallest possible error term in the Prime Number
Theorem. A single zero off the line would mean primes clump more irregularly
than anyone currently believes.

## What this project computes

Five things, all finite and all numerical.

**Identities.** That the functional equation and the symmetry `xi(s) = xi(1-s)`
hold at sample points, to the working precision. These are theorems; checking
them tests the implementation, not the mathematics.

**Zeros on the critical line.** Located by scanning Hardy's `Z(t)` for sign
changes and bisecting. `Z` is real-valued and vanishes exactly where
`zeta(1/2 + it)` does, so this is a reliable way to find on-line zeros — and
it is *only* able to find on-line zeros. It cannot see off the line at all.

**Zeros in the whole strip.** Counted by the argument principle: a contour
integral of `zeta'/zeta` around a rectangle enclosing the strip. This counts
every zero inside, wherever it sits, and is the only computation here that
could in principle detect a zero off the line.

**Spacing statistics.** Gaps between consecutive on-line zeros, unfolded to
unit mean and compared against the GUE prediction from random matrix theory.

**The explicit formula.** Partial sums over zeros, converging toward `psi(x)`
computed independently from prime powers.

## Why checking finitely many zeros cannot prove RH

This is the point the whole project is organised around, so it is worth
stating without hedging.

RH is a universally quantified statement about an infinite set. Any
computation terminates. Between those two facts there is no bridge.

Concretely: suppose every zero below height `T` is verified to lie on the
critical line, for some enormous `T`. What has been established is

> for all zeros `rho` with `0 < Im(rho) <= T`: `Re(rho) = 1/2`.

RH asserts the same with no upper bound. The gap between those is not
quantitative — it is not that a larger `T` gets you closer in any meaningful
sense. There are infinitely many zeros above any height you reach, and the
verified ones are a set of measure zero within them.

History is unkind here. The Mertens conjecture — that `|M(x)| < sqrt(x)` —
held for every value anyone could compute, was checked far beyond what
intuition suggested was necessary, and is **false**. Odlyzko and te Riele
disproved it in 1985 by showing

    limsup M(x) x^(-1/2) >  1.06        liminf M(x) x^(-1/2) < -1.009

which contradicts the conjecture. Their proof exhibits no counterexample at
all; in their own words, it "does not produce any single value of `x` for
which `|M(x)| > x^(1/2)`". None has been found since. They expected none below
`10^20`, "and maybe not even for `x < 10^30`" — comfortably past anything
computable then or now.

Had the Mertens conjecture been true it would have implied RH. Every
computation ever run on it agreed with it, and all of that agreement was
worthless as evidence.

So: numerical verification of RH to great heights is genuinely interesting,
genuinely hard, and genuinely not evidence that a proof is close.

## Three distinct things

Keeping these separate is what the project's status labels enforce.

**Numerical experimentation.** Computing zeros, spacings, partial sums, and
looking at them. Produces intuition and pictures. Establishes nothing.

**Finite-height verification.** Proving, for a specific `T`, that all zeros
below `T` are on the line. This is a real mathematical statement, and it has a
real proof technique behind it: compare a count of *all* zeros in the strip
below `T` against a count of zeros found *on the line* below `T`, and show they
are equal. Large-scale versions use Turing's method and interval arithmetic
so that every computed sign is provably correct. The published record — van de
Lune, te Riele and Winter for the first 1.5 billion zeros, and later
distributed efforts reaching much further — is of this kind.

**A proof of RH.** A finite argument covering all infinitely many zeros. Does
not exist.

This project does the first. It does the second in an *experimental* sense
only: the comparison of counts is implemented and does agree, but `mpmath`
gives no rigorous error bounds, so the output is strong evidence rather than a
certificate. The path from here to genuine finite-height certification is laid
out in [`certification.md`](certification.md). It does not attempt the third,
and neither does any amount of work along these lines.

## What would actually be needed

Nobody knows. The known equivalent formulations — Li's criterion, the
Nyman-Beurling criterion, the de Bruijn-Newman constant being `<= 0`, bounds on
the Mertens function — are reformulations, not footholds; each is exactly as
hard as RH itself. The Hilbert-Pólya program, which hopes to find a
self-adjoint operator whose eigenvalues are the zero ordinates, would give
reality of the eigenvalues for free, and has produced beautiful mathematics
and no such operator in a century of trying.

The honest summary is that RH is open, the approaches are known, and none of
them is close.
