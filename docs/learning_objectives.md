# Learning objectives and instructional design

## The idea the module is built around

Numerical agreement feels like confirmation. Often it is not. The specific
failure this module teaches is **circular validation**: checking a computed
quantity against something that was derived from the same assumptions, and
reading the inevitable agreement as independent support.

This failure is hard to teach in the abstract, because in most real settings
you cannot prove your check was circular — you can only suspect it. The
Riemann zeta function is unusual in that the circularity is *demonstrable*.
The module walks learners into the trap, lets them watch a check produce
agreement to 29 decimal places, and then shows — from library source code and
from a one-line identity — that the check could not have failed regardless of
where the zeros were.

The zeta function is the vehicle. The transferable skill is the destination.

## Audience and prerequisites

**Audience.** Scientists and engineers who validate computational models —
graduate students and researchers who write code that produces numbers they
then have to defend. No number theory background is assumed or needed.

**Prerequisites.**

- Python: comfortable installing a package and running a script.
- Mathematics: complex numbers, integration, and the idea of a function of a
  complex variable. Calculus is enough; complex analysis is not required to
  follow the argument, only to follow every proof in the optional exposition.
- No prior exposure to the Riemann Hypothesis.

**Time.** The core path (objectives 1–4) is about 90 minutes, including the
runtime of the scripts. The full path with the exposition is 6–8 hours and is
better spread over a week.

## Learning objectives

After working through the module, a learner should be able to:

1. **State what a numerical verification of a mathematical claim establishes,
   and what it cannot.** Specifically, articulate why verifying a property for
   every case below a finite bound gives no information about cases above it,
   using the Mertens conjecture as the concrete cautionary case — a conjecture
   that held for every computable value, implied the Riemann Hypothesis, and
   is false.

2. **Identify circular validation in a numerical check.** Given a check of the
   form "quantity A agrees with quantity B," determine whether B was derived
   from A or from A's assumptions. Learners do this on a worked case where the
   answer is provable: `mp.nzeros(T)` compared against the Riemann–von Mangoldt
   formula, where the reported discrepancy is `S(T)` *by definition*.

3. **Construct a check that is not circular.** Explain why counting zeros in a
   region by contour integration is independent of counting them along a line
   by sign changes, and — equally important — state precisely how far that
   independence extends (both computations bottom out in the same `mp.zeta`).

4. **Distinguish a numerical diagnostic from a proof of convergence.** Explain
   why a residual near zero is weak evidence that a quadrature converged while
   a large residual is strong evidence that it did not, and why the asymmetry
   matters when deciding whether to trust a result.

5. **Recognise when finite-precision arithmetic has silently limited a
   result.** Reproduce the precision-context failure documented in
   `numerical_methods.md`, in which a comparison performed at the ambient
   precision floors at `1e-16` no matter how accurate its inputs — a bug that
   appeared three separate times during this project's own development and was
   caught each time only because a check failed loudly.

6. **Assign honest status labels to computational claims.** Classify a result
   as proved, cited, numerical evidence, certified, or open, and justify why
   nothing produced by non-interval arithmetic can carry the fourth label.

## Suggested path

| Step | Activity | Time |
|---|---|---|
| 1 | Read `docs/mathematical_scope.md` — what RH says, why finite checks cannot settle it | 20 min |
| 2 | Run `scripts/verify_core_claims.py` and read section 2 of its output closely | 15 min |
| 3 | **Before reading further:** predict what section 2's residual column means. Then read `docs/claim_status.md` | 20 min |
| 4 | Work `notebooks/critical_line_exploration.ipynb` sections 4 and 5 | 30 min |
| 5 | Read `docs/numerical_methods.md` on residuals and precision contexts | 25 min |
| 6 | *Optional:* the full exposition in `site/critical-line-notebook.html` | 4–6 h |
| 7 | *Optional:* `docs/certification.md` — what rigour would actually require | 40 min |

Step 3 is the pedagogical hinge. Learners who read the explanation before
forming a prediction reliably report that the circularity was obvious; those
who predict first usually get it wrong. Instructors should enforce the
ordering.

## Exercises

1. `count_sign_changes_hardy_z` scans on a grid of step 0.05. Find a pair of
   zeros it would miss, using `spacing.lehmer_like_pairs` on the first 1000
   ordinates. What step would be needed? What would the resulting undercount
   look like in `experimentally_check_up_to` — silent, or detected?

2. `count_zeros_argument_principle_int` raises rather than rounding when the
   residual exceeds tolerance. Find a `T` where the default `maxdegree=8`
   fails and `maxdegree=10` succeeds. Explain what property of that particular
   `T` caused the failure. (Hint: evaluate `hardy_z` at it.)

3. Reproduce the precision-context bug deliberately: compute `xi(s)` and
   `xi(1-s)` at 30 digits but form `1-s` outside the precision context.
   Explain the floor you observe.

4. Take a validation step from your own field. Write down what it compares
   against what, and trace whether the comparison target was derived from the
   thing being validated. Most are; the interesting question is how you would
   tell.

Exercise 4 is the one that matters. The rest are practice for it.

## Assessment

There is no autograder, and a correct numerical answer is not the learning
outcome. Suggested assessment is a one-page written argument for exercise 4:
a validation procedure from the learner's own work, an honest analysis of
whether it is circular, and — if it is — a description of what a
non-circular alternative would require, including whether it is affordable.
