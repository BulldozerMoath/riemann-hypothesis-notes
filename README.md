# The Critical Line

**The first version of this project contained a validation check that could
not have failed.** Finding that, and working out what it implies, is what the
project is now about.

The check compared a computed count of the zeros of the Riemann zeta function
against the Riemann–von Mangoldt counting formula. The two agreed to
twenty-eight decimal places. That looked like confirmation. It was not: the
routine producing the count, `mpmath`'s `nzeros(T)`, is *built from* the same
formula — via `gram_index = floor(theta(t)/pi)` — so the difference it reports
is the function

$$S(T) = N(T) - \frac{\theta(T)}{\pi} - 1$$

which is *defined* as that difference. The check returns an identity. It would
have agreed to twenty-eight decimal places no matter where the zeros were,
including if every one of them had been off the critical line.

Worse, it looks *better* the harder you compute. The residual tracks the
arithmetic precision and nothing else — 7.3e-29 at 30 digits, 6.6e-49 at 50 —
which is the signature of this failure and is easy to mistake for rigour.

That error is now pinned down as a test
(`tests/test_counting.py::test_nzeros_is_not_independent_of_the_formula`), and
the rest of the repository is the work of doing it properly.

## The diagnostic

The general question, which is the portable part:

> For any check of the form *A agrees with B* — was B derived from A, or from
> A's assumptions? What value of A would have produced disagreement? If none
> exists, the check is reporting an identity, not evidence.

This applies far beyond zeta functions. Anywhere a model is validated against
a quantity that inherits its assumptions — modeled exposure in a risk
denominator, empirical-Bayes counts checked against the model that generated
them, random cross-validation on spatially autocorrelated data — the same
failure is available, usually without being visible in the diagnostics.

The zeta function is worth using as the teaching case for one reason: it is a
domain where the circularity is **provable** rather than merely suspected.
Usually, establishing that a check was circular needs exactly the ground truth
you were missing in the first place.

## Doing it properly

The honest version of a numerical check of the Riemann Hypothesis is an
**equality between two counts**:

- the number of zeros in the critical **strip** below height `T`, and
- the number of zeros found **on the critical line** below height `T`.

If they agree, every zero below `T` is on the line and simple. The first count
has to be computed without assuming where the zeros are, which
`count_zeros_argument_principle` does by contour-integrating `zeta'/zeta`
around a rectangle enclosing the strip — an integrand mentioning only zeta and
its derivative.

Both counts agree up to **T = 200** (79 in the strip, 79 on the line).

**NON-CERTIFIED NUMERICAL EVIDENCE ONLY.** These use different numerical
counting procedures, but both rely on finite-precision numerical evaluation
and therefore provide experimental evidence rather than a certified result.
They are not fully independent either: every value either count needs comes
out of `mp.zeta`, so a systematic error there could corrupt both at once
without the agreement revealing it.

For scale, Platt and Trudgian verified RH rigorously, using interval
arithmetic, to height 3·10¹². This project is ten orders of magnitude below
that and not rigorous. **It does not prove the Riemann Hypothesis and does not
claim to.** RH is open and nothing here bears on it.

No function in this package is named `certify`, because `mpmath` provides no
interval-certified arithmetic. The comparison routine is
`counting.experimentally_check_up_to`, and the name is the claim.

## Install

```bash
git clone https://github.com/BulldozerMoath/riemann-hypothesis-notes
cd riemann-hypothesis-notes

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .

pytest -q
```

118 tests, about 85 seconds.

## Run

```bash
python scripts/plot_hardy_z.py --start 10 --end 40           #  ~2 s
python scripts/verify_core_claims.py                         #  ~80 s
python scripts/explore_zeros.py --count 100                  #  ~27 s
python scripts/explicit_formula_demo.py --x 100 --zeros 100  #  ~27 s
python scripts/analyze_spacings.py --count 1000              #  ~10 min
```

Section 2 of `verify_core_claims.py` demonstrates the circularity directly —
it is the fastest way to see the point.

Every script takes `--dps` and `--output-dir`, labels its output
`[EXPERIMENTAL]`, and exits nonzero if a check misses its stated tolerance.
Figures land in `outputs/`.

Reaching T = 200 needs a longer run, and `--maxdegree 8` does not converge
that high — the script raises rather than rounding:

```bash
python scripts/verify_core_claims.py --max-height 200 --maxdegree 10   # ~9 min
```

## The exposition

Underneath all of this is an 11-stage derivation of the classical theory —
Dirichlet series through the Euler product, analytic continuation, the
functional equation, the Riemann–von Mangoldt formula, the explicit formula,
equivalent formulations of RH, and the modern research program. Two standalone
HTML pages in `site/`; open either in a browser, no server needed.

Each stage marks what is derived in full versus cited from a named classical
result. That labelling discipline is what made the original error findable.

## Layout

```
src/critical_line/     the package
  zeta_tools.py        zeta, functional equation, xi, precision management
  hardy_z.py           Riemann-Siegel theta, Hardy's Z
  zeros.py             grid-and-bisection zero finding
  counting.py          N(T) three ways — and which comparisons carry information
  gram.py              Gram points, Gram's law as a diagnostic
  spacing.py           gaps, unfolding, GUE comparison, Lehmer-like pairs
  explicit_formula.py  von Mangoldt Lambda, psi(x), primes from zeros
  plots.py             figures, each carrying the standing disclaimer
  constants.py         defaults, disclaimer text, reference data
scripts/               five runnable entry points
tests/                 118 tests
docs/                  scope, methods, claim status, certification, references
site/                  the two HTML exposition pages
notebooks/             guided walkthrough
paper.md               submission paper (see Status)
outputs/               generated figures (gitignored)
```

## Status labels

Every claim carries exactly one.

| Label | Meaning |
|---|---|
| Proved here | Full derivation in the exposition |
| Classical theorem | Published result, cited and attributed |
| Numerical evidence | Computed; true at the precision and range tested |
| Certified | Rigorous error bounds — **nothing here has this label** |
| Open | No claim made |

`mpmath` tracks how many digits it carries, not how many are correct. Reaching
*Certified* would need interval arithmetic (`python-flint` / `arb`) and
Turing's method with verified constants — see
[`docs/certification.md`](docs/certification.md), which is a roadmap, not an
implementation.

## Documentation

- [`docs/learning_objectives.md`](docs/learning_objectives.md) — objectives,
  prerequisites, suggested path, exercises, notes for instructors
- [`docs/mathematical_scope.md`](docs/mathematical_scope.md) — what RH says,
  and why checking finitely many zeros cannot settle it
- [`docs/numerical_methods.md`](docs/numerical_methods.md) — precision, zero
  searches, the argument principle, and why a small residual is evidence
  rather than a theorem
- [`docs/claim_status.md`](docs/claim_status.md) — per-claim status and the
  `nzeros` correction in full
- [`docs/certification.md`](docs/certification.md) — Turing's method; checking
  versus certifying
- [`docs/references.md`](docs/references.md) — sources for every classical
  result cited

## Other errors kept on the record

Three more defects were found during development and are documented where they
occurred rather than quietly fixed:

- A **precision-context error** that recurred three times: comparisons
  performed at the ambient precision floor near `1e-16` regardless of how
  accurate their inputs are. It once made a passing check report a spurious
  failure, and once made a correct identity look twelve orders of magnitude
  worse than it is.
- A **false independence claim** between `hardy_z` and `hardy_z_from_zeta`,
  presented as a cross-check. `mp.siegelz` only uses the Riemann–Siegel
  expansion above `|t| > 500 * prec` — far above any height used here — so
  below that it computes the same formula as the function it was being checked
  against.
- A **factual error** about where the first Mertens counterexample lies,
  stated in the wrong direction.

Each is the same species as the headline error, which is the point.

## Status

An earlier version of this work was prepared for submission to the Journal of
Open Source Education; `paper.md` and `paper.bib` remain in the repository.
JOSE has since paused submissions pending board deliberation on eligibility,
and the module has not yet been used in a taught course. It is complete and
usable for self-study in the meantime.

## What's mine and what isn't

The mathematics is classical or drawn from published work — Riemann, Hadamard,
von Mangoldt, Hardy, Selberg, Turing, Lehmer, Montgomery, Odlyzko–te Riele,
Rodgers–Tao — cited in [`docs/references.md`](docs/references.md). **None of it
originates here.**

What is mine is the exposition, the verification code, and the correction:
working the derivations through stage by stage, building numerics so every
claim has something checkable behind it, and then finding that one of those
checks was worthless.

I used Claude (Anthropic) as a derivation and writing partner throughout —
checking algebra, structuring exposition, writing code. The `nzeros`
circularity was identified during that collaboration by reading `mpmath`'s
source, after an external reviewer flagged the original wording as overstated.
Stated plainly rather than presented as unassisted work.

## Why this exists

I'm a civil engineering PhD student — transportation safety and pedestrian
research, not analytic number theory. I wanted to understand RH rather than
quote it. What I did not expect was that the most useful thing to come out of
it would be a lesson about validation that transfers directly back to my own
field, where models are routinely checked against quantities that inherit
their assumptions.

## License

Dual-licensed, as open educational materials should be: code under the MIT
License ([`LICENSE-CODE`](LICENSE-CODE)), prose and figures under CC BY 4.0
([`LICENSE-CONTENT`](LICENSE-CONTENT)). See [`LICENSE`](LICENSE) for which
applies where.
