# References

Sources for every classical result this project states, cited and not
re-derived. Where a constant or a numerical record is quoted, the entry says so
and says to check it against the source rather than against this file.

## Reference works

- **NIST Digital Library of Mathematical Functions**, Chapter 25 (Zeta and
  Related Functions). <https://dlmf.nist.gov/25>
  Definitions, the functional equation, the Riemann-Siegel formula
  (§25.10), and tabulated zero ordinates. The authority used for the
  ordinate values in `critical_line.constants.FIRST_ZERO_ORDINATES`.
- H. M. Edwards, *Riemann's Zeta Function*, Academic Press, 1974
  (Dover reprint 2001). The standard route through Riemann's 1859 memoir,
  including the explicit formula and the Riemann-Siegel formula.
- E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd edition
  revised by D. R. Heath-Brown, Oxford University Press, 1986.
- A. Ivić, *The Riemann Zeta-Function: Theory and Applications*, Wiley 1985
  (Dover reprint 2003).

## Riemann's memoir

- B. Riemann, "Ueber die Anzahl der Primzahlen unter einer gegebenen Grösse",
  *Monatsberichte der Berliner Akademie*, November 1859.
  The functional equation, the completed `xi`, the explicit formula, and the
  hypothesis itself.

## The Riemann-Siegel formula

- C. L. Siegel, "Über Riemanns Nachlaß zur analytischen Zahlentheorie",
  *Quellen und Studien zur Geschichte der Mathematik, Astronomie und Physik*
  2 (1932), 45-80.
  Siegel's recovery of the asymptotic formula from Riemann's unpublished
  papers — the reason `Z(t)` is computable at large height at all.
- NIST DLMF §25.10 for the modern statement.

## The Riemann-von Mangoldt formula

- H. von Mangoldt, "Zu Riemanns Abhandlung 'Ueber die Anzahl der Primzahlen
  unter einer gegebenen Grösse'", *Journal für die reine und angewandte
  Mathematik* 114 (1895), 255-305.
  The counting function `N(T) = theta(T)/pi + 1 + S(T)` and the proof of
  Riemann's explicit formula.
- R. J. Backlund, "Über die Nullstellen der Riemannschen Zetafunktion",
  *Acta Mathematica* 41 (1918), 345-375.
  Explicit bounds on `S(T)`; `mpmath`'s `backlunds` is named for this work.

## Hardy's theorem

- G. H. Hardy, "Sur les zéros de la fonction `zeta(s)` de Riemann",
  *Comptes Rendus de l'Académie des Sciences* 158 (1914), 1012-1014.
  Infinitely many zeros lie on the critical line. Note what this does **not**
  say: it is consistent with infinitely many lying off it.
- A. Selberg, "On the zeros of Riemann's zeta-function",
  *Skrifter Norske Videnskaps-Akademi Oslo* I (1942), no. 10.
  A positive proportion of zeros lies on the line.
- J. B. Conrey, "More than two fifths of the zeros of the Riemann zeta
  function are on the critical line", *Journal für die reine und angewandte
  Mathematik* 399 (1989), 1-26.

## Gram's law

- J.-P. Gram, "Sur les zéros de la fonction `zeta(s)` de Riemann",
  *Acta Mathematica* 27 (1903), 289-304.
- J. I. Hutchinson, "On the roots of the Riemann zeta-function",
  *Transactions of the American Mathematical Society* 27 (1925), 49-60.
  Records the failure of Gram's law at `n = 126`.
- Gram's law is known to fail infinitely often; see Titchmarsh §10 and the
  discussion in Edwards §6.5. It is an empirical tendency, not a theorem, and
  not a proof technique.

## Turing's method

- A. M. Turing, "Some calculations of the Riemann zeta-function",
  *Proceedings of the London Mathematical Society* (3) 3 (1953), 99-117.
  The method of bounding the integral of `S(t)` rather than `S(t)` itself.
- R. S. Lehman, "On the distribution of zeros of the Riemann zeta-function",
  *Proceedings of the London Mathematical Society* (3) 20 (1970), 303-320.
- R. P. Brent, "On the zeros of the Riemann zeta function in the critical
  strip", *Mathematics of Computation* 33 (1979), 1361-1372.
- T. S. Trudgian, "Improvements to Turing's method",
  *Mathematics of Computation* 80 (2011), 2259-2279.
  Sharpened explicit constants in the `S_1` bound. **The constants have been
  revised more than once; take them from the paper, not from a summary,** and
  check the stated range of validity in `t` before applying them.

## Lehmer's phenomenon

- D. H. Lehmer, "On the roots of the Riemann zeta-function",
  *Acta Mathematica* 95 (1956), 291-298.
- D. H. Lehmer, "Extended computation of the Riemann zeta-function",
  *Mathematika* 3 (1956), 102-108.
  The close pair near `t = 7005`. Lehmer pairs are unusually close
  consecutive zeros on the critical line. They are not counterexamples to the
  Riemann Hypothesis, but they make naive zero searches more difficult.

## Numerical verification of zeros

Quoted here only where the source states the figure; do not extrapolate.

- J. van de Lune, H. J. J. te Riele, D. T. Winter, "On the zeros of the
  Riemann zeta function in the critical strip. IV",
  *Mathematics of Computation* 46 (1986), 667-681.
  The first 1,500,000,001 zeros verified to lie on the critical line.
- X. Gourdon, "The `10^13` first zeros of the Riemann zeta function, and zeros
  computation at very large height", 2004. Unpublished preprint, widely
  circulated; it has not been through peer review, so treat the figure as a
  reported computation rather than a certified record.
- D. J. Platt, "Isolating some non-trivial zeros of zeta",
  *Mathematics of Computation* 86 (2017), 2449-2467.
  A rigorous interval-arithmetic implementation — the model for what
  "certified" means in this context.
- D. J. Platt, T. S. Trudgian, "The Riemann hypothesis is true up to
  `3 * 10^12`", *Bulletin of the London Mathematical Society* 53 (2021),
  792-797. A peer-reviewed rigorous verification.

## Equivalent formulations and related results

- A. M. Odlyzko, H. J. J. te Riele, "Disproof of the Mertens conjecture",
  *Journal für die reine und angewandte Mathematik* 357 (1985), 138-160.
  The cautionary tale: a conjecture implying RH, supported by all available
  numerical evidence, and false.
- E. Bombieri, J. C. Lagarias, "Complements to Li's criterion for the Riemann
  hypothesis", *Journal of Number Theory* 77 (1999), 274-287.
- B. Nyman, *On the One-Dimensional Translation Group and Semi-Group in
  Certain Function Spaces*, thesis, Uppsala, 1950; A. Beurling, "A closure
  problem related to the Riemann zeta-function",
  *Proceedings of the National Academy of Sciences* 41 (1955), 312-314.
- B. Rodgers, T. Tao, "The de Bruijn-Newman constant is non-negative",
  *Forum of Mathematics, Pi* 8 (2020), e6.

## Random matrix theory and zero statistics

- H. L. Montgomery, "The pair correlation of zeros of the zeta function",
  in *Analytic Number Theory*, Proceedings of Symposia in Pure Mathematics 24,
  American Mathematical Society, 1973, 181-193.
  The pair correlation conjecture — a **conjecture**, and one that would not
  imply RH if proved.
- A. M. Odlyzko, "On the distribution of spacings between zeros of the zeta
  function", *Mathematics of Computation* 48 (1987), 273-308.
  The numerical agreement with GUE statistics.
- M. V. Berry, J. P. Keating, "The Riemann zeros and eigenvalue asymptotics",
  *SIAM Review* 41 (1999), 236-266.

## Software

- F. Johansson and others, *mpmath: a Python library for arbitrary-precision
  floating-point arithmetic*. <https://mpmath.org/>
  The `nzeros` and `zetazero` implementations discussed in
  [`claim_status.md`](claim_status.md) are in
  `mpmath/functions/zetazeros.py`.
- F. Johansson, "Arb: efficient arbitrary-precision midpoint-radius interval
  arithmetic", *IEEE Transactions on Computers* 66 (2017), 1281-1292.
  The library a certified version of this project would need.
