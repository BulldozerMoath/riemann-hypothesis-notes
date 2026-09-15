---
title: 'The Critical Line: teaching circular validation through a case where the circularity is provable'
tags:
  - Python
  - numerical verification
  - reproducibility
  - validation
  - arbitrary-precision arithmetic
  - Riemann zeta function
authors:
  - name: Moath Alshannaq
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Civil Engineering, University of Texas at Arlington, USA
    index: 1
date: 15 September 2026
bibliography: paper.bib
---

# Summary

`The Critical Line` is a self-contained learning module about a specific way
numerical validation fails. A researcher checks a computed quantity against a
reference, obtains agreement to many decimal places, and concludes the
computation is sound — when the reference was derived from the same
assumptions, so agreement was guaranteed regardless of whether the computation
was right.

This failure is difficult to teach because in most settings it cannot be
demonstrated, only suspected. The module uses a case where it can be proved.
Learners run a check comparing `mpmath`'s zero-counting function `nzeros(T)`
against the Riemann–von Mangoldt counting formula, observe agreement to
twenty-eight decimal places or better, and then establish from library source
code and a one-line
identity that the check could not have failed for *any* configuration of
zeros. They then build a comparison that is not circular, and examine exactly
how far its independence extends.

The module comprises an installable Python package (118 tests), five
command-line scripts, a Jupyter walkthrough, figure generation, and
documentation stating the status of every claim. The Riemann zeta function is
the vehicle; validation practice is the subject.

# Statement of need

Computational researchers are trained to validate models and rarely trained to
audit the validation itself. Circular validation — comparing a result against a
quantity derived from the same assumptions — is common and largely invisible,
because it produces exactly the agreement that signals success.

Existing material does not address this well. Numerical analysis texts cover
error and conditioning; reproducibility training covers environments and
version control. Neither teaches a learner to ask whether a passing check
*could have failed*. The obstacle is pedagogical: convincing examples are hard
to construct, since demonstrating circularity usually requires access to a
ground truth that is unavailable precisely when it matters.

The Riemann zeta function resolves this. Its zero-counting function satisfies
the exact identity $N(T) = \theta(T)/\pi + 1 + S(T)$, where $S(T)$ is *defined*
as that difference. Any check comparing a $\theta$-derived count against the
formula therefore reports $S(T)$ and nothing else — provably, not arguably.
Learners can verify this themselves at any precision they choose.

A second property makes the domain unusually suitable: verifying the Riemann
Hypothesis below a finite height is a real, published research activity
[@lune:1986; @platt:2021], while the hypothesis itself remains open. The gap
between "checked to $3 \cdot 10^{12}$" and "proved" is concrete and
unbridgeable, so learners encounter the limits of finite verification as
current practice rather than as an abstraction. The Mertens conjecture
supplies the cautionary case: it agreed with all available computation,
implied the Riemann Hypothesis, and is false [@odlyzko:1985].

# Learning objectives and instructional design

Learners who complete the module can state what finite verification
establishes and what it cannot; identify circular validation in a check;
construct a non-circular comparison and state the limits of its independence;
distinguish a numerical diagnostic from a proof of convergence; recognise
where finite-precision arithmetic has silently capped a result; and assign
honest status labels to computational claims. Full objectives, prerequisites,
timings and exercises are in `docs/learning_objectives.md`.

The core path takes roughly 90 minutes and assumes no number theory. Its
pedagogical hinge is ordering: learners predict what the residual column of
`verify_core_claims.py` means *before* reading the explanation. Those who read
first report the circularity as obvious; those who predict first usually get
it wrong.

The design commitment is that the module's own claims meet the standard it
teaches. No function is named `certify`, because `mpmath` provides no
interval-certified arithmetic [@mpmath]; the comparison routine is
`experimentally_check_up_to`. The argument-principle counter raises an
exception rather than rounding a result whose quadrature residual exceeds
tolerance. Where the module reaches agreement between a contour-integral count
and a sign-change count at height $T = 200$, it labels the result
non-certified numerical evidence and notes that both computations rest on the
same underlying `mp.zeta`. A genuinely certified result would require interval
arithmetic and Turing's method with verified constants [@turing:1953;
@trudgian:2011; @platt:2017; @johansson:2017]; `docs/certification.md` sets
out what that would take and does not pretend it is implemented.

Three defects found during development are retained as teaching material,
including a precision-context error that recurred three times and a false
independence claim between two functions that turned out to share a code path.
Each is documented where it occurred.

# Experience of use

The module has not yet been used in a taught course. It is complete and
immediately usable for self-study, and has been through an independent
reproduction audit — clean-environment installation, full test suite, all
scripts end to end, and verification of every cited source. Classroom piloting
is planned; the author welcomes reports from instructors who adopt it.

# Acknowledgement of AI assistance

The exposition, package and this paper were developed with Claude (Anthropic)
as a derivation and writing partner: checking algebra, structuring
explanations, and drafting code. The author directed the work, verified the
mathematics against the cited sources, and is responsible for all content. The
`mp.nzeros` circularity central to the module was identified during that
collaboration by reading `mpmath`'s source. This is stated because the module
teaches honest reporting of how results were obtained.

# Availability

Source, documentation and exposition: <https://github.com/BulldozerMoath/riemann-hypothesis-notes>.
Code is MIT-licensed; documentation and figures are CC BY 4.0.

# References
