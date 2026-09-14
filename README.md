# The Critical Line — expository notes on the Riemann Hypothesis

This repository contains a set of study notes I built while learning the
classical theory behind the Riemann Hypothesis (RH), plus a companion page
that checks the notebook's claims against real numerical computation.

**This is not a proof of RH, and it does not claim to be.** The Riemann
Hypothesis is an open problem — one of the seven Millennium Prize Problems —
and nothing here resolves it. What this repo *is*: a from-scratch,
step-by-step derivation of the standard mathematics that everyone studying
RH has to learn first (the analytic continuation of ζ(s), the functional
equation, the Euler product, the explicit formula connecting primes to
zeros, the equivalent formulations of RH, and a tour of the modern research
program), written out with full proofs where a full proof is feasible at
this level, and honest citations to the literature where it isn't.

## Contents

- **`critical-line-notebook.html`** — the main notebook, "The Critical
  Line." Open it in any browser (no server needed). It's organized into 11
  stages, from the basic definition of ζ(s) through the explicit formula,
  equivalent forms of RH (Mertens, Li's criterion, Nyman–Beurling, the
  de Bruijn–Newman constant), and the Hilbert–Pólya / random-matrix-theory
  research program. Each stage marks clearly what is derived in full versus
  cited from a named classical result (e.g. Poisson summation, Stirling's
  formula, Perron's formula).
- **`ten-thousand-zeros.html`** — a numerical companion, "Ten Thousand
  Zeros." It uses `mpmath` to compute real zeros of ζ(s) (out past the
  10,000th one) and checks them against the notebook's claims: the
  functional equation, the Riemann–Siegel Z function, the zero-counting
  formula N(T), the GUE level-repulsion statistics, and the explicit
  formula's convergence to ψ(x).
- **`verify.py`** — a standalone Python script with the same three core
  numerical checks (functional equation, N(T) cross-check, explicit-formula
  convergence), so the claims can be reproduced outside the browser. Needs
  `mpmath` (`pip install mpmath`).

## What's original here and what isn't

The mathematics is all classical or, for the later stages, drawn from
20th- and 21st-century published results (Riemann, Hadamard, von Mangoldt,
Hardy, Selberg, Montgomery, Odlyzko–te Riele, Rodgers–Tao, and others — the
notebook cites them by name at the relevant stage). What's original is the
exposition: working through the derivations myself, stage by stage, and
building the numerical verification to check every claim against direct
computation.

I used Claude (Anthropic) as a step-by-step derivation and writing partner
while building this — checking my algebra, catching malformed equations,
and helping structure the exposition and the verification code. I'm
disclosing that plainly rather than presenting it as unassisted work.

## Why this exists

I'm a civil engineering PhD student (transportation safety / pedestrian
research) — pure analytic number theory isn't my field. I wanted to
actually understand RH rather than just quote it, so I worked through the
standard path into the problem and kept an honest record of how far genuine
derivation gets you versus where you have to start citing research papers.
If it's useful to anyone else learning the same material, that's a bonus.

## License

Feel free to reuse, adapt, or build on this for your own study. Attribution
appreciated but not required.
