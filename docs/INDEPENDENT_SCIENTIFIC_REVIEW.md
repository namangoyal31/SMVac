# Independent Scientific Review of SMVac

**Subject:** SMVac, branch `formalization`, HEAD `3d514e2` (repository by Naman Goyal, IISc Bangalore).

**Nature of this document.** This is a read-only review written as a skeptical referee would before
publication. It was not written to defend the project. No source code, dataset, figure or git object was
changed. The review's own calculations were run in a temporary directory outside the repository and are
described in enough detail here to be repeated from the stored data.

**Notation for evidence.**

* **[E]**: verified directly. This covers a computation run for this review, a quotation checked in the
  primary source, or a repository fact.
* **[I]**: inference or judgement by the reviewer.

Equation numbers of other papers refer to their arXiv versions.

> **Addendum: status after the final correction stage.** This review is kept as the record of the assessment of
> HEAD `3d514e2`. Its recommendations 1 and 2 in §8 were incorporated into the repository afterwards:
>
> * attribution to Wipf (1985), Espinosa (2020) and Di Luzio et al. (2016), with a notation map in
>   `docs/SCIENTIFIC_METHODOLOGY.md` §6.5a;
> * the unresummed first- and second-order formulas added to a tested, reproducible comparison,
>   `validation/compare_approximations.py`, with output in `validation/results/approximation_comparison.json`;
> * the README narrative rewritten;
> * the "only numerics are reliable" and "every finite-action configuration" statements corrected;
> * the provenance wording replaced by the version suggested in §1.3.
>
> The new script reproduces the numbers quoted in §2.4 from the committed data:
>
> * maximum first-order error 1.54% over the 21,626 standard-tier controlled points ($q=0.338$; `tree`: 1.41%);
> * 2.65% on the parabolic model for $0.1\le q<1$, and 10.2% at $q=0.05$;
> * 1.7% on the SM lines for $q\ge0.097$, and 14.9% at $q=0.044$.
>
> One detail is added. Under a consistent metric, in which an infinite estimate counts as an infinite error,
> $S_{\rm curv}$ and $S_{\rm CGM}$ are infinite at one `tree` point, $(M_h,M_t)=(120.1,168.5)$ GeV with $q=0.173$.
> The 42% and 46% quoted in §2.4 are the largest *finite* errors.

---

## 1. Executive assessment

### 1.1 Bottom line

SMVac is a careful, well-validated **numerical benchmark** of the tree-level bounce action of the RG-improved
Standard Model potential, $V=\lambda_P(\phi)\phi^4/4$, against closed-form approximations. Its numerical
reference is credible. One independent check is new in this review: the reference reproduces the published
**second-order** analytical result of Espinosa (2020) to about $10^{-3}$ in the coefficient (§2.2). The
software engineering, provenance and reproducibility are of professional quality.

The project's **analytical results are largely already in the literature, or are elementary consequences of
known methods.** The literature search for this review found:

* The thick-wall coefficient $\frac{13}{72}$ is Espinosa (2020), eq. (3.21), in identical form. The same paper
  has the Fubini-weighted average $6\int_0^1x(1-x)f(x\phi_0)\,dx$, the optimal centre $\phi_0=\mu_0e^{5/6}$,
  the scaling (Affleck) constraint, and the full second-order term, eq. (6.5).
* The lower bound $S\ge 8\pi^2/(3|\lambda_{\min}|)$ follows directly from a known method: Wipf (1985) combines
  the Coleman–Glaser–Martin inf-max characterisation with the sharp Sobolev constant. A weaker general bound,
  $24/\lambda_{\rm cr}$, is Sato & Takimoto (2018). Di Luzio, Isidori & Ridolfi (2016) remark that the
  constant-$\lambda$ formula acts as a lower bound, and they had already solved the RG-improved tree-level SM
  bounce numerically.
* The variational ("CGM") estimate with a Fubini trial is a resummed form of Espinosa's first-order formula. The
  use of trial functions in the inf-max characterisation to obtain upper bounds is Wipf (1985).

The **finiteness criterion** $q>\tfrac12{\rm Var}(\ln g)$ and the two-parameter threshold $3/32$ are correct.
This review proves them rigorously, with a closed form, in §2.3. I did not find them in the literature, but they
are elementary and of limited scope.

The review also contains a **finding that changes one of the project's main conclusions**. The plain,
*unresummed* first-order formula $S\simeq S_{\rm FL}\,[1+13/(72q)]$, Espinosa's eq. (3.21), is far more
accurate near the stability boundary than any of SMVac's "improved" estimates:

* about 1.5% at worst over all standard-tier points of both controlled scans;
* within 2.7% on the parabolic model for every $q\ge0.1$.

This was found by evaluating that formula on SMVac's own stored data (§2.4). SMVac's improved estimates
($S_{\rm curv}$, $S_{\rm CGM}$, $S_{\rm CGM2}$) are resummations that diverge. The README statement that only the
numerical solution is reliable for $0.05\lesssim q\lesssim1$ is therefore not supported as written.

**Is there a defensible research contribution beyond good software?** A modest one, if it is reframed
honestly. It would be a numerical-benchmark paper on the accuracy of analytical approximations to the SM bounce
from thick to thin walls, with:

1. confirmation of the known first- and second-order results;
2. the exact explanation of why fixed-shape variational and resummed estimates fail;
3. the empirical observation that the unresummed first-order formula survives deep into the near-critical
   region;
4. the thin-wall crossover and the Planck-plateau phenomenology, as a property of the frozen-coupling model.

It would **not** support claims of new analytical theory for the thick-wall regime, nor of improving the
Buttazzo et al. calculation. SMVac addresses one component only, the tree-level action of the RG-improved
potential, and not their one-loop treatment (§5).

### 1.2 Summary classification of the principal results (details in §6)

| result | classification |
|---|---|
| Numerical reference (accuracy $\sim10^{-12}$–$10^{-14}$, standard tier) | 1. established, independently supported |
| $S\ge S_{\rm FL}$ for every finite-action bounce | 2. correct; follows directly from known methods (Wipf 1985) |
| $(S-S_{\rm FL})/S_{\rm FL}\to\frac{13}{72q}$ | 2. correct; already published (Espinosa 2020, eq. 3.21) |
| Second-order term | 2. published (Espinosa 2020, eq. 6.5); SMVac reproduces it, a strong numerical check |
| $S_{\rm curv}$, $S_{\rm CGM}$ as "new" estimates | 2. resummed and variational forms of known results; not improvements near criticality |
| Finiteness criterion $q>\frac12{\rm Var}(\ln g)$; $13/72$; $3/32$ | 3. potentially new, but elementary; proved here |
| Identity $S=\frac{\pi^2}{2}\lambda''M\,(q-\frac12{\rm Var}_b\ln\phi)$ | 3. elementary corollary of the Derrick and Affleck constraints |
| Universality of $S/S_{\rm FL}-1$ in $q$ (2–3%) | 4. empirical |
| Error thresholds and crossovers in $q$ and $\delta$ | 4. empirical |
| "Only numerics are reliable for $0.05\lesssim q\lesssim1$" | 5. not supported (see §2.4) |
| "Plateau points metastable since $S\ge S_{\rm FL}$ for every finite-action configuration" | 5. incorrectly stated; a correct formulation exists (§2.1) |
| No-finite-bounce region and Planck-sensitive bounces | 6. outside the validity of the physical model (frozen couplings) |
| FL error is small compared with RG-improvement ambiguities in the decay criterion | 1/2. consistent with the literature, within the tree-level framework |

### 1.3 Repository state and provenance

**Git state [E].**

* Branch `formalization`, HEAD `3d514e2fc0a0252ae36348a9c11319c2dab730a8`, working tree clean at the start of
  this review. This document is the only new file.
* History: `a7c070a` (original; `main` = `origin/main` = `a7c070a`, untouched), then `6b7ae5a`, `1d0d435`,
  `ad63002`, `3d514e2`.
* The message-only rewrite maps `f3d2e25`→`6b7ae5a`, `bb4b5b7`→`1d0d435` and `6fdb8aa`→`ad63002`. For each pair
  the tree hashes are identical (checked with `git rev-parse <c>^{tree}`).
* The backup ref `refs/original/refs/heads/formalization` still points to `6fdb8aa`.
* No commit message contains a co-author trailer.
* `formalization` has no remote-tracking branch.

**Scan provenance.** The committed scan files carry `git_commit = bb4b5b70e9…-dirty`. Four states need to be
kept apart.

1. **Clean committed base at scan time:** `bb4b5b7`, now `1d0d435` with an identical tree. It has
   `BounceOptions::max_nodes = 2048`. **[E]**
2. **Uncommitted modifications compiled into the scan binary.** The development session record shows the
   following sequence. **[E]**
   * At 00:52:53 IST (19:22:53 UTC) a single `sed` changed `max_nodes` from 2048 to 1024, followed by
     `cmake --build build`.
   * `build/smvac_scan.exe` has mtime 00:53:04 IST.
   * No other compiled source was modified between the commit (00:39:40 IST) and that build.
   * The three committed scans ran 19:25:07–19:58:06 UTC (scan headers).
   * During the scans, two accessor methods were added to headers (01:15 IST). A stagnation-threshold
     experiment in `src/bounce.cpp` was compiled only into a scratch directory and then reverted (01:16–01:17).
   * The binary mtime shows that `build/` was not rebuilt afterwards.
   * Earlier scan runs started at 00:31 and 00:38 IST were killed, and their outputs discarded.
3. **The final committed source `ad63002`.** Relative to `1d0d435`, the scan-relevant C++ differs only by the
   `max_nodes = 1024` line and the two unused accessors. The CMake file also gained an optional target.
   `git diff --ignore-cr-at-eol 1d0d435 ad63002 -- src include apps python/smvac/run_scan.py configs`. **[E]**
4. **Reproduction.** A clean CMake build of `ad63002` recomputed 233 seeded, stratified points. All matched
   (`validation/results/reproduce_sample.csv`). **[E]**
   * 227 scan points (83 `cw1`, 36 `tree`, 108 exploratory; about 0.5% of 42,622), on status, $S_{\rm num}$,
     $S_{\rm FL}$ and $S_{\rm CGM}$ at $10^{-12}$ relative.
   * 4 shooting-fallback points at $10^{-10}$.
   * 2 near-critical model points at $10^{-9}$.

**Assessment [I].** The evidence firmly supports the statement that the scans were produced by code
functionally identical to the scan-relevant code in `ad63002`. Two caveats:

* The evidence for item 2 is the session record and filesystem timestamps, **neither of which is part of the
  repository**.
* The reproduction is a 0.5% sample. It is evidence for, not proof of, reproducibility of every point.

The current README sentence "The scans therefore come from the source of `ad63002`" is slightly too strong.
Suggested wording:

> *Provenance.* The committed scans (`results/*.csv.gz`) were produced on 2026-11-08, 19:25–19:58 UTC, by a binary
> built at 19:23 UTC from commit `bb4b5b7` (renamed `1d0d435` by a later message-only edit; identical tree) plus one
> uncommitted change, `BounceOptions::max_nodes = 1024` (`include/smvac/bounce.hpp`), which is included in `ad63002`.
> The scan-relevant source of `ad63002` differs from that build only by two accessor methods that the scan code
> does not call. The build state is reconstructed from the development record and file timestamps; the scan
> headers themselves record only `bb4b5b7…-dirty`. A seeded, stratified sample of 227 scan points recomputed from a
> clean build of `ad63002` reproduced the stored values to $10^{-12}$; this supports, but does not prove,
> reproducibility of all 42,622 points.

The other stored datasets were produced from the working tree that became `ad63002` **[E]**:

* The fallback CSV was generated 01:43–02:03 IST. The Python code it uses is unchanged in `ad63002`.
* The near-critical CSVs were generated 02:05–02:17 IST with `analytic.py` in its final state, except for a
  later refactoring that extracted the `cgm2_trial` function. The refactor was functionally identical (equal
  values at $q=2.5$ to all printed digits), and 2 model points were recomputed with the final code.
* The SM-line $S_{\rm CGM2}$ values were not recomputed with the final code.

---

## 2. Mathematical verification

All derivations below were reconstructed from the equations, independently of the project's documents.
Conventions:

* $V=\lambda(\phi)\phi^4/4$, false vacuum at $\phi=0$, $x=\ln\phi$, $\lambda''=d^2\lambda/dx^2$.
* With $\mu=\phi$, $\lambda''=\beta_\lambda'\equiv d\beta_\lambda/d\ln\mu$. In $t=\ln\mu^2$ units,
  $\lambda''=4\,d^2\lambda/dt^2$.
* $S=T+V$, with $T=\tfrac12\int|\nabla\phi|^2d^4x$ and $V=\int\tfrac14\lambda\phi^4d^4x$.

### 2.1 The lower bound $S_{\rm FL}=8\pi^2/(3|\lambda_{\min}|)$

**Derivation (verified).**

1. *Pohozaev/Derrick identity.* For a finite-action critical point in $d=4$, $2T+4V=0$, so $T=-2V$ and
   $S=T/2=-V$.
2. *Pointwise bound.* If $\lambda(\phi)\ge\lambda_{\min}$ on the range of the bounce, then
   $-V\le\frac{|\lambda_{\min}|}{4}\int\phi^4$.
3. *Sharp Sobolev inequality in $\mathbb R^4$* (Aubin 1976; Talenti 1976): $\int\phi^4\le K(\int|\nabla\phi|^2)^2$
   with $K=3/(32\pi^2)$. Equality holds only for $\phi=a/(b+|x-x_0|^2)$. The constant was re-verified on the
   Fubini profile: $K=0.009498860966469$ against $3/(32\pi^2)=0.009498860966469$. **[E]**
4. Combining these:
   $S=-V\le\frac{|\lambda_{\min}|}{4}K\big(\int|\nabla\phi|^2\big)^2=\frac{|\lambda_{\min}|}{4}K(2T)^2=|\lambda_{\min}|K(2S)^2=4K|\lambda_{\min}|S^2$.
   Hence
   $$S\ge\frac{1}{4K|\lambda_{\min}|}=\frac{8\pi^2}{3|\lambda_{\min}|}.$$

**Exact assumptions [I].**

* $\phi\in D^{1,2}(\mathbb R^4)$ (finite $T$, $\phi\to0$ at infinity), $V(\phi)\in L^1$, and $\phi$ a weak
  solution.
* Decay $\phi=O(r^{-2})$ and $|\nabla\phi|=O(r^{-3})$. These hold for a massless false vacuum by the linearised
  equation, and they make the boundary terms of the Pohozaev multiplier argument vanish.
* The naive scaling argument needs care because $x\cdot\nabla\phi\notin L^2$ for $r^{-2}$ tails. The multiplier
  argument on balls avoids this.
* $V$ need only be continuous with a piecewise-$C^1$ derivative, so the derivative kink at $M_{\rm Pl}$ is
  harmless. O(4) symmetry is not required.
* For the Higgs doublet, the bound holds with $|\Phi|$ via Kato's inequality $|\nabla|\Phi||\le|\nabla\Phi|$.
* $\lambda_{\min}$ may be the minimum over the range of the bounce. The global minimum, including the frozen
  regions, gives a weaker but valid bound. The implemented potential satisfies this: $\lambda_P$ is frozen below
  the matching scale and above $M_{\rm Pl}$.

**Equality [I].** Equality requires both the Aubin–Talenti profile and $\lambda(\phi(x))=\lambda_{\min}$ wherever
$\phi\ne0$, i.e. $\lambda\equiv\lambda_{\min}$ on $(0,\phi_0]$. For any running $\lambda$ the inequality is
therefore strict. The project's statement is correct.

**Applicability to the implemented problem [I].** The bound is a theorem about $V=\lambda_P(\phi)\phi^4/4$ with
its false vacuum at $\phi=0$, i.e. with the electroweak mass term and the vev dropped. This is standard; Di Luzio
et al. (2016) discuss its validity. Nothing is lost in going from constant to running $\lambda$, because only
$\lambda\ge\lambda_{\min}$ is used. The bound says nothing about the one-loop decay exponent.

**A misstatement found [E].** `docs/VALIDATION_REPORT.md` §5.2 says "$S\ge S_{\rm FL}$ for every finite-action
configuration". That is false: the action is unbounded below on general configurations. The correct statement,
which still supports the classification of the plateau points, is as follows.

* For any configuration with $V<0$, the dilatation maximum $\max_aS(\phi_a)=T^2/(4|V|)\ge S_{\rm FL}$, by the
  same Sobolev step.
* Any path from $\phi=0$ to a configuration with $S<0$ must cross the Pohozaev manifold $\{T+2V=0\}$. Near
  $\phi=0$, Sobolev gives $T+2V\ge T-2K|\lambda_{\min}|T^2>0$. At the endpoint, $T+2V<T+V<0$.
* On that manifold, $S=T^2/(4|V|)\ge S_{\rm FL}$.
* Hence the mountain-pass barrier is at least $S_{\rm FL}$.

The physics caveat remains: $\Lambda_B\le M_{\rm Pl}$ and the frozen-coupling model are both assumed.

**Prior art [E].**

* Wipf (1985) derives lower bounds on the bounce action from the inf-max (CGM) characterisation combined with
  Hölder–Young and Sobolev inequalities "using the sharpest constants" (his eqs. 15–19). The SMVac bound is a
  direct special case of this method for the critical quartic.
* Sato & Takimoto (2018) give the general bound $S\ge24/\lambda_{\rm cr}$,
  $\lambda_{\rm cr}=\max[-4V/|\phi|^4]$. They note that the Fubini value $8\pi^2/(3\lambda)\simeq26.3/\lambda$
  exceeds it. The sharp-Sobolev argument gives the stronger constant $8\pi^2/3$ for exactly this class.
* Di Luzio et al. (2016) note that $8\pi^2/(3|\lambda(1/R_{\rm SM})|)$ "provides a lower bound" (heuristically,
  by expanding about the minimum).
* I did not find the sharp statement written explicitly. It should be presented as an immediate application of
  known tools, not as a new theorem.

### 2.2 The thick-wall formula $S/S_{\rm FL}-1\sim13/(72q)$

**Derivation (verified).**

* The Fubini family $\phi_h=h\,w(r/R)$, with $w=1/(1+s^2)$ and $R^2=8/(|\lambda_{\min}|h^2)$, is degenerate under
  dilatations at $\lambda=\lambda_{\min}$.
* At first order, the correction is the perturbation evaluated on the family. With $\phi^4r^3dr\propto
  w(1-w)\,dw$:
  $$\delta S(h)=\frac{8\pi^2}{3\lambda_{\min}^2}\,\big\langle\lambda(\ln h+\ln w)-\lambda_{\min}\big\rangle,\qquad
  \langle f\rangle=6\int_0^1w(1-w)f\,dw .$$
* The dilatation zero mode is fixed by stationarity in $h$.
* For quadratic $\lambda$, $\langle\cdot\rangle=\frac{\lambda''}{2}\langle(\ln h+\ln w-x_0)^2\rangle$. This is
  stationary (a minimum) at $\ln h=x_0-\langle\ln w\rangle=x_0+\frac56$, where it equals
  $\frac{\lambda''}{2}{\rm Var}(\ln w)$.
* Since $w\sim{\rm Beta}(2,2)$, ${\rm Var}(\ln w)=\psi'(2)-\psi'(4)=\frac{13}{36}$. Hence
  $$\frac{S}{S_{\rm FL}}=1+\frac{13}{72}\frac{\lambda''}{|\lambda_{\min}|}+O(q^{-2}).$$
* Checks: $\langle\ln w\rangle=-0.8333333$ and $\frac12{\rm Var}=0.1805556$ by quadrature. **[E]**

**Status [I].** This is a *formal* first-order asymptotic result: degenerate perturbation theory about the Fubini
family. It is not an empirical fit, and it is not proved with error bounds. It is confirmed numerically: the ratio
of $q(S/S_{\rm FL}-1)$ to $13/72$ is 0.99992 at $q=1000$ and 0.99293 at $q=10$
(`near_critical_summary.json`). **[E]**

For non-quadratic $\lambda$, as in the SM, the general first-order correction is the full average
$\langle\lambda-\lambda_{\min}\rangle$. The coefficient $13/72$ captures only its quadratic part. The higher
Taylor terms of $\lambda$ are of the same order in $1/|\lambda_{\min}|$, because the Fubini profile spans
$O(1)$ e-folds (standard deviation of $\ln w$ = 0.60).

**Prior art [E].** Espinosa (2020), §3, takes exactly $\lambda=-\lambda_0+\frac12\beta'_\lambda\ln^2(\phi/\mu_0)$ and
finds $S=\frac{8\pi^2}{3\lambda_0}+\frac{13\pi^2}{27\lambda_0^2}\beta'_\lambda+O(\beta'^2)$ (eq. 3.21). Since
$\frac{13\pi^2/27}{8\pi^2/3}=\frac{13}{72}$, this is identical to SMVac's formula.

The same paper gives, in matching form:

* $\phi_0=\mu_0e^{5/6}$ (eq. 3.17);
* the Fubini average $\langle f\rangle_{\phi_0}=6\int_0^1x(1-x)f(x\phi_0)\,dx$ (eq. 3.25);
* $S\simeq\langle8\pi^2/(3|\lambda(\phi)|)\rangle_{\phi_0}$ (eq. 3.27);
* the scaling constraint $\int\phi_B^4\beta_\lambda(\phi_B)r^3dr=0$ (eq. 3.23), which is SMVac's dilatation
  identity;
* the second-order term (eq. 6.5).

**Second order, a new cross-check of SMVac's numerics [E].**

* Evaluating Espinosa's eq. (6.5) at $\kappa=0$ and $L_0=5/6$ gives the second-order term. The $O(\beta')$ shift
  of $L_0$ enters only at third order, because $S$ is stationary in $L_0$. With
  $K_0=\int_0^1(1-6x+6x^2){\rm Li}_2(x)^2dx=0.1187177$, the bracket is $B(5/6)=-0.32746$. In SMVac's convention,
  $S/S_{\rm FL}=1+\frac{13}{72q}+\frac{c_2}{q^2}$ with $c_2=B/24=-0.013644$.
* SMVac's stored model data give $[(S/S_{\rm FL}-1)-\frac{13}{72q}]\,q^2=-0.01363$ (q = 1000), $-0.01361$ (300),
  $-0.01355$ (100), $-0.01333$ (30). This converges to $-0.013644$ with an $O(1/q)$ trend.
* The fitted $S_{\rm CGM}$ coefficient (0.0462) equals the analytic value $(13/72)^2-c_2=0.04624$.

The residual left after the first-order term agrees with the analytic second-order term to about $10^{-3}$. At
$q=1000$ that is about $10^{-11}$ in $S$. This independently validates SMVac's reference through second order.
It also shows that the thick-wall expansion is fully known.

### 2.3 The variational finiteness criterion

**Statement.** For a trial $\phi=h\,g(r/R)$ of fixed shape, maximising over $R$ gives $T^2/(4|V|)$, with
$V=\frac{\pi^2}{2}h^4R^4\int\lambda(\ln h+\ln g)\,g^4s^3ds$. For $\lambda=\lambda''[\frac12(x-x_0)^2-q]$, write
$a=\ln h-x_0$ and average with respect to $g^4s^3ds$:
$$V\propto\tfrac12\big\langle(a+\ln g)^2\big\rangle-q,\qquad \min_a\tfrac12\langle(a+\ln g)^2\rangle=\tfrac12{\rm Var}(\ln g).$$
Hence **$V<0$ for some $h$ if and only if $q>q_c[g]\equiv\frac12{\rm Var}(\ln g)$** **[E]** (algebra re-derived).

This is necessary and sufficient, provided $T$, $\int g^4s^3ds$ and the first two log-moments are finite. It
holds exactly only for quadratic $\lambda$; in the SM the poles shift slightly with the cubic terms, as the data
confirm.

**Fubini shape ($n=1$):** $q_c=\frac{13}{72}$, as in §2.2. Then $S_{\rm CGM}=\frac{8\pi^2}{3\lambda''(q-13/72)}$,
which is exactly the resummed first-order formula $S_{\rm FL}/(1-\frac{13}{72q})$.

**Two-parameter family $w_n=(1+s^{2n})^{-1/n}$: closed form (new in this review) [E].**

* Substitute $t=s^{2n}$. The weight $w^4s^3ds\propto t^{2/n-1}(1+t)^{-4/n}dt$ is a beta-prime distribution
  $\beta'(\frac2n,\frac2n)$. So $1/(1+t)\sim{\rm Beta}(\frac2n,\frac2n)$, and $\ln w=\frac1n\ln\frac1{1+t}$. Hence
  $$q_c(n)=\frac{\psi'(2/n)-\psi'(4/n)}{2n^2}.$$
* Values: $q_c(1)=\frac{13}{72}$, $q_c(2)=\frac18$. These match direct quadrature to 12 digits for
  $0.25\le n\le100$, including $n<1$. At $n=10^4$ they match to $2\times10^{-9}$, limited by the quadrature.
* *Monotonicity, proved.* With $z=2/n$ and $\psi'(z)=\sum_k(z+k)^{-2}$:
  $$q_c=\frac18\Big[\tfrac34+\sum_{k\ge1}g_k(z)\Big],\qquad g_k=\frac{z^2}{(z+k)^2}-\frac{z^2}{(2z+k)^2},\qquad g_k'(z)=2zk\big[(z+k)^{-3}-(2z+k)^{-3}\big]>0.$$
  So $q_c$ strictly decreases in $n$. As $n\to\infty$, $g_k\to0$ and $q_c\to\frac{3}{32}$. As $n\to0$,
  $q_c\to\infty$.
* *Limit member.* $w=\min(1,s^{-2})$ has half its weight at $\ln w=0$ and half distributed as $-{\rm Exp}(2)$, so
  ${\rm Var}=\frac3{16}$. In closed form, $S_\infty=4\pi^2/[\lambda''(q-3/32)]$ at $a=\frac14$. This was
  re-derived: $\min_a(a^2-\frac a2+\frac14)=\frac3{16}$ at $a=\frac14$.

**Conclusion.** For exactly quadratic $\lambda$, $S_{\rm CGM2}<\infty$ iff $q>3/32$, because $q_c(n)>3/32$ for
every finite $n$. **The divergence at $3/32$ is therefore mathematically established, not only observed through
the optimiser.** **[E]**

The implementation restricts $0.5\le n\le1000$. Since $q_c(1000)-3/32=2.4\times10^{-9}$, this is immaterial.

**Exact bounce identity.** For quadratic $\lambda$, the dilatation (Affleck) constraint gives
$\langle\ln\phi\rangle_b=x_0$. With Derrick's identity, $V=-T/2$, this yields
$$S=\tfrac{\pi^2}{2}\lambda''\Big(\int f^4d\tau\Big)\big(q-\tfrac12{\rm Var}_b\ln\phi\big)>0 .$$
This is correct, and SMVac's check 8 confirms $\langle\ln\phi\rangle_b-x_0\sim10^{-14}$. **[E]** It is a one-line
corollary of Espinosa's eq. (3.5), the same constraint, and Derrick's identity. **[I]**

**Significance [I].** For the FL family, the criterion simply locates the pole of the resummed first-order
formula, which is evident from $S_{\rm FL}/(1-13/(72q))$. Its general form, finiteness iff the trial's
log-field spread is below $2q$, is a clean and correct way to explain why *any* fixed-shape or Fubini-tailed
family must fail. The $3/32$ result belongs to an ad hoc family. I did not find the criterion in the physics
literature reviewed (§4). Given its elementary nature, a prior appearance in some form cannot be excluded.

### 2.4 Assessment of the "improved estimates"

The stated statuses are correct:

* $S_{\rm FL}$ is a rigorous lower bound (§2.1).
* $S_{\rm CGM}$ and $S_{\rm CGM2}$ are upper estimates, rigorous only if the CGM infimum is attained.
  Attainment was not verified for the RG-improved potential. With $\lambda_{\rm Pl}<0$ the potential has
  *critical* growth $-|\lambda_{\rm Pl}|\phi^4$ at large field, where compactness can fail (Brezis–Nirenberg
  1983). The inequality $S_{\rm num}\le S_{\rm CGM2}\le S_{\rm CGM}$ holds at all stored points where both are
  finite. **[E]**
* $S_{\rm curv}$ is not a bound. **[E]** $S_{\rm curv}<S_{\rm num}$ at the SM central point; it is larger on the
  parabolic model.

**The reviewer's comparison with the published unresummed formulas [E].** These use the stored
`near_critical_model.csv` and `near_critical_sm.csv`, the controlled scans, and $c_2=-0.013644$.

*Relative errors on the parabolic model:*

| $q$ | $S_{\rm FL}$ | 1st order $S_{\rm FL}(1+\frac{13}{72q})$ | 2nd order | $S_{\rm CGM}$ | $S_{\rm CGM2}$ |
|---|---|---|---|---|---|
| 10 | $-1.8\times10^{-2}$ | $+1.3\times10^{-4}$ | $-9\times10^{-6}$ | $+4.5\times10^{-4}$ | $+1.0\times10^{-4}$ |
| 1 | $-0.15$ | $+6.7\times10^{-3}$ | $-4.9\times10^{-3}$ | $+4.1\times10^{-2}$ | $+8.7\times10^{-3}$ |
| 0.35 | $-0.33$ | $+2.2\times10^{-2}$ | $-5.4\times10^{-2}$ | $+0.39$ | $+7.2\times10^{-2}$ |
| 0.2 | $-0.46$ | $+2.7\times10^{-2}$ | $-0.16$ | $+4.5$ | $+0.26$ |
| 0.1 | $-0.64$ | $+8.8\times10^{-4}$ | $-0.49$ | $\infty$ | $+7.1$ |
| 0.07 | $-0.73$ | $-4.0\times10^{-2}$ | $-0.79$ | $\infty$ | $\infty$ |
| 0.05 | $-0.81$ | $-0.10$ | $-1.16$ | $\infty$ | $\infty$ |

*Standard-tier controlled-scan points (both prescriptions, 21,626 points):* the unresummed first-order formula
has median $|{\rm err}|=9\times10^{-6}$ and maximum 1.5%. For comparison, the maxima are 34% (`cw1`) and 50%
(`tree`) for $S_{\rm FL}$, 42% for $S_{\rm curv}$ and 46% for $S_{\rm CGM}$ (`cw1`).

*Near-critical SM lines:* the first-order formula stays within 1.7% for $q\ge0.097$. It is −3.2% at $q=0.084$
(`cw1`, 120) and −15% at $q=0.044$.

**Interpretation [I].**

1. In the bulk, $S_{\rm CGM}$ is the most accurate, with median error $2\times10^{-7}$.
2. Near criticality, the simplest *published* formula beats every estimate the project introduced, by one to
   three orders of magnitude. It stays finite and few-per-cent accurate far outside its formal domain. At small
   $q$ it behaves like $q^{-2}$, close to the observed local exponent of $-1.6$ to $-1.95$.
3. There is no theoretical justification for this accuracy. It is an empirical observation. The second-order
   truncation is *worse* near the boundary, as expected for an asymptotic series.
4. The README statements "both improved estimates break down… only the numerical solution is reliable" and "no
   closed-form asymptotics cover $0.05\lesssim q\lesssim1$" are true *for the project's resummed and variational
   estimates*. They are not true for closed forms in general and need revision.

---

## 3. Numerical-method assessment

**Formulation [I].** The variables $\tau=\ln r$ and $f=r\phi$ turn the scale-covariant bounce into a homoclinic
problem on a finite interval with asymptotic Robin conditions $f'=\pm f$. This choice is natural and well
suited to Chebyshev collocation. The neglected term in the boundary conditions is
$O(e^{-2|\tau-\tau_c|})$. Doubling the half-width from 10 to 25 changes $S$ by less than $10^{-14}$. **[E]**

**Resolution and convergence [E].**

* Exponential convergence in $N$ for thick walls (plateau for $N\ge192$).
* Thin walls need more nodes: on the parabolic model at $q=0.1$, $6\times10^{-5}$ at $N=256$ falls to
  $8\times10^{-13}$ at $N=768$.
* Acceptance requires both exact identities, Derrick and dilatation, below $10^{-9}$. These are monitored, not
  imposed.
* Off-node residuals lie at the expected round-off floor. Clenshaw–Curtis and Simpson quadrature agree.
* Long-double re-solution differs by $\le4\times10^{-15}$. One-ulp perturbations of the potential coefficients
  change $S$ by $\le3\times10^{-13}$.

These checks are independent of solver agreement and adequately support the quoted accuracy. The second-order
agreement with Espinosa (§2.2) is a further, fully independent analytic check.

**Uncertainty estimates [I].** The production $\delta S$ (two-resolution difference plus identity residual times
$S$) is a heuristic but conservative estimator. It is consistent with the tests above.

**Planck kink [I].** Freezing $\lambda$ above $M_{\rm Pl}$ makes $G=\lambda+\lambda'/4$ discontinuous, and
spectral convergence becomes algebraic. The "reduced-accuracy" tier (identities $<10^{-5}$, $\delta S/S\sim10^{-6}$)
rests on the following support:

* one direct cross-check with $N=1536$ against shooting, agreeing to $10^{-7}$;
* agreement between the adaptive shooting solvers on the 78 fallback bounces, to $5\times10^{-10}$.

This is adequate. However, the *physical* error at these points (the arbitrary freezing) is uncontrolled and
orders of magnitude larger.

**Solver agreement [I].** The documents state correctly that agreement between solvers is not proof. The
C++ and Python implementations share no code, and their RG inputs agree to $10^{-10}$. A common error in the
problem definition (the potential) is excluded only by the input checks against Buttazzo et al. eq. (61).

**Near-critical region [E].** The deep thin-wall regime ($q\lesssim0.05$, $S\gtrsim5\times10^7$) is unresolved
on the SM lines. This is documented. The model points at $q=0.07$ and $0.05$ rest on M4 alone, with M2 failing,
but the identities hold to $10^{-9}$.

**No-finite-bounce classification [I].** The criteria (every shot of both solvers undershoots, and the FL
variational optimum runs to the edge of the search range) are numerical evidence, not a proof of non-existence.
A rigorous framework exists:

* For critical growth at large field, a least-action bounce exists when the variational infimum lies
  *strictly below* the bubble level $8\pi^2/(3|\lambda_{\rm Pl}|)$. Otherwise minimising sequences concentrate
  (Brezis–Nirenberg 1983; concentration–compactness). The plateau points sit at
  $S_{\rm inf}/S_{\rm FL}\le1.024$, consistent with this picture.
* Espinosa (2020) notes that $\beta_\lambda$ must change sign over the bounce for a bounce to exist.
* Andreassen et al. (2017) note that the rate is infinite if $\lambda(\mu)$ has no minimum.
* Buttazzo et al. (2013) prescribe $\lambda(\Lambda_B)\to\lambda(M_{\rm Pl})$ for $\Lambda_B>M_{\rm Pl}$, which
  matches SMVac's plateau infimum.

**Precision versus accuracy [I].** Numerical precision ($10^{-12}$) far exceeds physical accuracy:

* the prescription difference is 2.7% of $S$;
* the omitted one-loop corrections are $O(1)$ units of action;
* gravitational corrections are $\delta S_G=256\pi^3\Lambda_B^2/(45|\lambda|M_{\rm Pl}^2)$ (Buttazzo et al., §6.3).

The high precision is meaningful only for comparing approximations to the *same* mathematical problem. That is
exactly how SMVac uses it.

---

## 4. Literature and prior-art review

The search covered physics journals and arXiv, using web searches, with full text read where marked. It was not
exhaustive. In particular, the mathematical literature on variational bounds for critical elliptic problems was
sampled, not surveyed.

| # | reference | what it establishes | relation to SMVac |
|---|---|---|---|
| 1 | J. R. Espinosa, *Vacuum decay in the Standard Model: analytical results with running and gravity*, JCAP 06 (2020) 052, [arXiv:2003.06219](https://arxiv.org/abs/2003.06219) (full text read) | §3: model $\lambda=-\lambda_0+\frac12\beta'\ln^2(\phi/\mu_0)$; perturbed Fubini bounce; eq. (3.17) $\phi_0=\mu_0e^{5/6}$; eq. (3.21) $S=\frac{8\pi^2}{3\lambda_0}+\frac{13\pi^2}{27\lambda_0^2}\beta'$; eqs. (3.25)–(3.29) Fubini average; eq. (3.23) scaling constraint; footnote 3: near-critical bounce is thin-walled and not Fubini-like; §6, eq. (6.5): second order, including gravity | **Contains SMVac's 13/72 result, the averaging construction and the second-order term.** SMVac's numerics confirm eq. (6.5) |
| 2 | A. Wipf, *Upper and lower bounds for the bounce action*, J. Phys. A 18 (1985) 2521 ([PDF](https://www.tpi.uni-jena.de/qfphysics/homepage/wipf/publications/papers/bounceaction.pdf), read) | inf-max characterisation of the bounce action; lower bounds via Hölder–Young and Sobolev inequalities with sharpest constants (eqs. 15–19); upper bounds from trial functions | method behind SMVac's lower bound and variational upper estimates |
| 3 | R. Sato, M. Takimoto, *Absolute lower bound on the bounce action*, PRL 120 (2018) 091802, [arXiv:1707.01099](https://arxiv.org/abs/1707.01099) | $S\ge24/\lambda_{\rm cr}$ (eq. 17) for general multi-field potentials; notes the Fubini value $26.3/\lambda$ | weaker general bound; SMVac's sharp-Sobolev bound improves the constant for this class |
| 4 | L. Di Luzio, G. Isidori, G. Ridolfi, *Stability of the electroweak ground state in the Standard Model and its extensions*, PLB 753 (2016) 150, [arXiv:1509.05028](https://arxiv.org/abs/1509.05028) (read) | numerical tree-level bounce of $V=\lambda(h)h^4/4$ with running MS-bar $\lambda$ at the SM point; comparison with $8\pi^2/(3\vert\lambda(1/R)\vert)$ (eqs. 39–42); heuristic lower-bound remark | prior numerical solution of the same problem at a single point |
| 5 | A. Andreassen, W. Frost, M. D. Schwartz, *Scale invariant instantons and the complete lifetime of the Standard Model*, PRD 97 (2018) 056006, [arXiv:1707.08124](https://arxiv.org/abs/1707.08124) (read in part) | consistent treatment of the dilatation collective coordinate; running enters at one loop through $\ln(\mu R)$ terms (eq. 3.33); "the leading prediction must start at one loop" | SMVac's tree-level RG-improved bounce is a different resummation from the consistent expansion (§5) |
| 6 | D. Buttazzo et al., *Investigating the near-criticality of the Higgs boson*, JHEP 12 (2013) 089, [arXiv:1307.3536](https://arxiv.org/abs/1307.3536) (read) | §6.3, eqs. (80)–(83): $S=8\pi^2/(3\vert\lambda(\Lambda_B)\vert)$ "replaced by the one-loop action" of ref. 7; $\wp_0=0.15\Lambda_B^4e^{-S}/H_0^4$; Planck prescription; eqs. (64), (66) stability bounds | the comparison target (§5) |
| 7 | G. Isidori, G. Ridolfi, A. Strumia, *On the metastability of the Standard Model vacuum*, NPB 609 (2001) 387, [hep-ph/0104016](https://arxiv.org/abs/hep-ph/0104016) (abstract) | full one-loop decay rate around the Fubini bounce | the one-loop action used by Buttazzo et al.; not computed in SMVac |
| 8 | S. Chigusa, T. Moroi, Y. Shoji, PRL 119 (2017) 211801, [arXiv:1707.09301](https://arxiv.org/abs/1707.09301); PRD 97 (2018) 116012, [arXiv:1803.03902](https://arxiv.org/abs/1803.03902) (abstracts) | state-of-the-art one-loop decay rate in the SM | beyond SMVac's tree-level scope |
| 9 | P. Baratella et al., *Revising the full one-loop gauge prefactor in electroweak vacuum stability*, [arXiv:2406.05180](https://arxiv.org/abs/2406.05180) (abstract) | updated full one-loop rate; $\mu=1/R$; Planck cut-off on the $R$ integral | current state of the art |
| 10 | I. Affleck, *On constrained instantons*, NPB 191 (1981) 429, [doi](https://doi.org/10.1016/0550-3213(81)90307-2) | origin of the scaling (Affleck) constraint for broken scale invariance | SMVac's dilatation identity |
| 11 | K. Lee, E. J. Weinberg, *Tunneling without barriers*, NPB 267 (1986), [doi](https://doi.org/10.1016/0550-3213(86)90150-1) | decay in scale-invariant / barrierless potentials | context for Fubini–Lipatov decay |
| 12 | J. R. Espinosa, *Tunneling without bounce*, PRD 100 (2019) 105002, [arXiv:1908.01730](https://arxiv.org/abs/1908.01730) | decay without a bounce; action minimum pushed to infinite field ("pseudo-bounces") | the same mechanism as SMVac's Planck-plateau points |
| 13 | J. R. Espinosa, *A fresh look at the calculation of tunneling actions*, [arXiv:1805.03680](https://arxiv.org/abs/1805.03680); J. R. Espinosa, T. Konstandin, JCAP 01 (2019) 051, [arXiv:1811.09185](https://arxiv.org/abs/1811.09185) | tunnelling-potential formulation in which the action is a *minimum*; trial functions give upper bounds | an alternative route to rigorous upper estimates |
| 14 | S. Coleman, V. Glaser, A. Martin, CMP 58 (1978) 211 | least-action bounce characterised variationally; O(4) symmetry | basis of $S_{\rm CGM}$ and its attainment caveat |
| 15 | H. Brezis, L. Nirenberg, CPAM 36 (1983) 437, [doi](https://doi.org/10.1002/cpa.3160360405) | existence for critical Sobolev growth when the infimum lies below the bubble level | the rigorous framework for existence and non-existence on the plateau |
| 16 | T. Aubin, J. Diff. Geom. 11 (1976) 573; G. Talenti, Ann. Mat. Pura Appl. 110 (1976) 353 (bibliographic details as cited in the repository; not re-verified here) | sharp Sobolev constant and extremals | ingredient of the lower bound |
| 17 | G. Isidori, V. S. Rychkov, A. Strumia, N. Tetradis, PRD 77 (2008) 025034, [arXiv:0712.0242](https://arxiv.org/abs/0712.0242); A. Salvio, A. Strumia, N. Tetradis, A. Urbano, JHEP 09 (2016) 054, [arXiv:1608.02555](https://arxiv.org/abs/1608.02555) | gravitational corrections to SM vacuum decay | not included in SMVac; relevant for Planck-sensitive points |
| 18 | V. Branchina, E. Messina, PRL 111 (2013) 241801, [arXiv:1307.5193](https://arxiv.org/abs/1307.5193) | Planck-suppressed operators can change the lifetime | context for the legacy $\phi^6$ term |
| 19 | A. Aravind, D. Lorshbough, S. Paban, PRD 89 (2014) 103535, [arXiv:1401.1230](https://arxiv.org/abs/1401.1230) | lower bounds for multi-field quartic bounce actions | related bounds literature |
| 20 | B. Garbrecht, P. Millington, PRD 98 (2018) 016001, [arXiv:1804.04944](https://arxiv.org/abs/1804.04944) | radiative corrections around the Fubini–Lipatov instanton | related |
| 21 | Numerical bounce solvers: C. L. Wainwright, [arXiv:1109.4189](https://arxiv.org/abs/1109.4189) (CosmoTransitions); A. Masoumi, K. D. Olum, B. Shlaer, [arXiv:1610.06594](https://arxiv.org/abs/1610.06594); P. Athron et al., CPC 244 (2019) 448, [arXiv:1901.03714](https://arxiv.org/abs/1901.03714); V. Guada, M. Nemevšek, M. Pintar, [arXiv:2002.00881](https://arxiv.org/abs/2002.00881); S. Chigusa, T. Moroi, Y. Shoji, [arXiv:1906.10829](https://arxiv.org/abs/1906.10829); R. Sato, PRD 101 (2020) 016012, [arXiv:1907.02417](https://arxiv.org/abs/1907.02417) | shooting, polygonal, gradient-flow and multi-field methods | SMVac's spectral-collocation approach is standard numerical analysis; not a methodological novelty |

**On the 13/72 and 3/32 results [E/I].**

* $13/72$ exists in equivalent form: Espinosa 2020, eq. 3.21.
* I found no occurrence of the general finiteness criterion $q>\frac12{\rm Var}(\ln g)$, of the closed form
  $q_c(n)$, or of $3/32$. Searches combined the variational bounce estimate with trial-profile finiteness,
  thresholds and log-quadratic couplings.
* Because these are elementary, the absence of a citation is weak evidence of novelty. A specialist check would
  still be needed before any claim.

---

## 5. Comparison with Buttazzo et al. (2013)

**What Buttazzo et al. compute [E].**

* NNLO matching at $M_t$ (eqs. 55–60).
* Three-loop RGEs (App. B).
* Two-loop $\lambda_{\rm eff}$ with the anomalous-dimension factor, for the stability bound.
* For the lifetime, $S(\Lambda_B)=8\pi^2/(3|\lambda(\Lambda_B)|)$ with the tree-level action "replaced by the
  one-loop action" of Isidori et al. (2001).
* $\Lambda_B$ chosen to maximise $\Lambda_B^4e^{-S}$, which "roughly amounts to minimising $\lambda(\Lambda_B)$".
* $\wp_0=0.15\Lambda_B^4e^{-S}/H_0^4$ with $H_0=67.4$ km/s/Mpc (eqs. 80–83).
* For $\Lambda_B>M_{\rm Pl}$, only a lower bound on the tunnelling probability, obtained with
  $\lambda(\Lambda_B)=\lambda(M_{\rm Pl})$.

**What SMVac tests [I].** SMVac compares $S_{\rm FL}=8\pi^2/(3|\min\lambda_P|)$ with the exact *tree-level* bounce of
the potential $V=\lambda_P(\phi)\phi^4/4$ for the same $\lambda_P$, at $\mu=\phi$. Two prescriptions are used:

* `tree`: the MS-bar $\lambda$. Here $S_{\rm FL}$ coincides with the tree-level part of Buttazzo's formula,
  because $\min_\mu\lambda(\mu)=\lambda(\Lambda_B)$ at $\beta_\lambda=0$.
* `cw1`: the one-loop CW $\lambda_{\rm eff}$, without the factor $e^{4\Gamma}$ and without two-loop terms.

Ingredient by ingredient [E]:

* *Inputs.* Matching uses the same fits. RG loop orders are the same. Action normalisation is the same, with
  $V=\lambda h^4/4$ for the real radial field. The metastability criterion is the same formula.
* *Not the same.*
  1. The effective potential: no two-loop $\lambda_{\rm eff}$ and no $e^{4\Gamma}$.
  2. The quantum corrections to the action: the one-loop action of Isidori et al. (2001) is absent.
  3. The scale identification for the prefactor: SMVac uses $\Lambda_B=1/R$ from the numerical bounce, and also
     reports $\arg\min\lambda$.

**Quantitative cross-check of inputs [E].**

* At $M_t=173.34$ GeV, SMVac's absolute-stability line lies at $M_h^{\rm crit}=129.74$ GeV (`cw1`) and 129.84
  GeV (`tree`). These are interpolated from the stored scans. Buttazzo's eq. (64) gives $129.6\pm0.3$ GeV, so
  the agreement is within their quoted theoretical error.
* At $M_h=125.15$ GeV, SMVac gives $M_t^{\rm crit}=171.12$ GeV (`cw1`). Buttazzo's eq. (66) quotes
  $171.53\pm0.42$ GeV, while linear inversion of their eq. (64) gives 171.1 GeV. Their two equations differ by
  about 0.4 GeV at this point. I have not resolved this, so the comparison is inconclusive at the 0.4 GeV level.
* $\lambda$ at $M_{\rm Pl}$ reproduces eq. (61).

**Which conclusions concern the analytical approximation [I].**

* The measured $S_{\rm num}-S_{\rm FL}$ is a statement about one approximation: replacing the tree-level bounce
  of the RG-improved potential by the Fubini action at $\lambda_{\min}$. That statement is genuine and
  well-founded.
* In the consistent expansion (Andreassen et al. 2017; Isidori et al. 2001), the running of $\lambda$ across the
  bounce enters through the one-loop $\ln(\mu R)$ terms and the collective-coordinate integral. SMVac's tree-level
  "running correction" $\frac{13}{72q}S_{\rm FL}$ is $O(\beta'/\lambda^2)$. That is the same parametric order as
  pieces of the one-loop correction that Buttazzo et al. include and SMVac does not.
* Therefore SMVac has **not** shown that Buttazzo's lifetime is in error by $7\times10^{-4}$ in the exponent. It
  has shown that, within the RG-improved tree-level scheme, the Fubini formula is accurate to that level.
* The finding that the FL error (0.0025 GeV shift of the $\wp=1$ boundary) is much smaller than the
  prescription and scale ambiguities (0.09–0.12 GeV) concerns the tree-level scheme's internal ambiguities. It is
  consistent with the literature's view that a consistent one-loop treatment is needed. It is not a correction to
  the published result.
* SMVac does **not** improve the Buttazzo calculation. It quantifies one component of its approximation chain.

---

## 6. Novelty assessment

| # | result | class | basis |
|---|---|---|---|
| 1 | Validated numerical reference for the tree-level RG-improved bounce across the $(M_h,M_t)$ plane | 1 | §3; second-order analytic agreement (§2.2); reproducibility (§1.3). Single-point prior work: Di Luzio et al. 2016 |
| 2 | $S\ge S_{\rm FL}$ (sharp, rigorous for finite-action bounces) | 2 | Wipf 1985 method; weaker Sato–Takimoto 2018; heuristic remark in Di Luzio et al. 2016 |
| 3 | $\frac{13}{72q}$ thick-wall law and the Fubini-weighted average | 2 | Espinosa 2020, eqs. 3.21, 3.25–3.29 |
| 4 | Second-order coefficient | 2 (numerically confirmed by SMVac, class 1 as a check) | Espinosa 2020, eq. 6.5; §2.2 |
| 5 | $S_{\rm curv}$, $S_{\rm CGM}$ | 2 | resummation and variational form of #3; Wipf 1985 trial-function bounds |
| 6 | $S_{\rm CGM2}$; finite iff $q>3/32$; closed-form $q_c(n)$ | 3 | proved here (§2.3); ad hoc family; elementary |
| 7 | General finiteness criterion $q>\frac12{\rm Var}(\ln g)$ | 3 | elementary; not found in the literature reviewed |
| 8 | Bounce identity $S=\frac{\pi^2}{2}\lambda''M(q-\frac12{\rm Var}_b)$ | 3 (elementary corollary) | Derrick + Affleck constraint (Espinosa eq. 3.5) |
| 9 | Universality of $S/S_{\rm FL}-1$ in $q$ to 2–3% over $0.17\le q\le600$ | 4 | stored scans; expected from the structure of the expansion |
| 10 | Thresholds (1%/10%), crossovers $q^*=0.39$ ($S_{\rm CGM}$) and $0.16$ ($S_{\rm CGM2}$); SM $q/\delta\approx121$–$141$ GeV$^{-1}$ | 4 | stored data |
| 11 | Unresummed first-order formula accurate to $\le1.5$% over the scans (reviewer finding) | 4 | §2.4; unexplained |
| 12 | "Only numerics reliable for $0.05\lesssim q\lesssim1$" | 5 | contradicted by #11 |
| 13 | "Metastable by the bound" for plateau points, as worded in VALIDATION §5.2 | 5 as worded | corrected argument in §2.1 |
| 14 | No-finite-bounce region; Planck-sensitive bounces | 6 | frozen couplings; numerical classification only |
| 15 | Legacy construction is a first-order expansion that is unbounded below | 1 (internal) | reproduced analytically and numerically |

**Answers to the questions asked.**

* **What is genuinely new?**
  * A systematic, reproducible numerical map of the error of the Fubini formula for the RG-improved tree-level SM
    bounce, from thick to thin walls and across the controlled $(M_h,M_t)$ window, with explicit near-critical
    and Planck-plateau classification.
  * The explicit finiteness criterion and the closed-form two-parameter threshold. Both are minor and
    elementary.
  * The empirical observation that the unresummed first-order formula stays accurate deep into the near-critical
    regime. This came from this review, using the project's data.
* **What is useful numerical confirmation of existing theory?**
  * Espinosa's first- and second-order thick-wall results, confirmed to about $10^{-3}$ in the second-order
    coefficient.
  * The lower bound.
  * The Fubini-average construction.
  * The qualitative near-critical breakdown of the Fubini picture, which Espinosa noted in his footnote 3.
* **What would make a defensible paper?** A benchmark study, for example in PRD, JCAP or EPJC, titled roughly
  *"How accurate are analytical approximations to the electroweak bounce near criticality?"*. It would:
  1. confirm the known expansion to second order;
  2. show quantitatively why resummed and fixed-shape variational estimates fail, through the finiteness
     criterion;
  3. document and attempt to explain the robustness of the unresummed first-order formula;
  4. characterise the thin-wall crossover;
  5. state the RG-scheme context honestly.

  Without items 3 and 4, or a new uniformly valid approximation, the contribution would be thin.
* **Is the finiteness criterion publishable alone as a short note?** In my judgement, no. It is a few lines of
  algebra. For the FL family it locates the obvious pole of the resummed first-order formula, and the 3/32 result
  concerns an ad hoc family. It works well as a lemma inside the benchmark paper above.
* **What would make the work more substantial?**
  1. A uniformly valid analytical approximation through the crossover, with a stated error. Candidates are
     matched asymptotics between the perturbed-Fubini and thin-wall regimes (allowing for the power-law tail of a
     massless false vacuum), or a trial family that decouples the tail amplitude from the core.
  2. An explanation, or bounds, for the accuracy of the unresummed first-order formula.
  3. A rigorous existence criterion for the frozen-coupling potential (Brezis–Nirenberg type), replacing the
     numerical no-bounce classification.
  4. Comparison with a consistent one-loop rate (e.g. ELVAS, Chigusa et al., Baratella et al.), to measure the
     scheme dependence that SMVac currently only bounds by the `cw1`/`tree` difference.
  5. Gravitational corrections near the Planck-sensitive region (Espinosa 2020; Isidori et al. 2008).
  6. Current inputs ($M_t$, $\alpha_s$) and propagation of parametric uncertainties.

---

## 7. Scientific limitations

1. **Tree level only.** There is no one-loop determinant, so the decay exponent is not computed to the precision
   of the literature, and "errors of the literature formula" are scheme-internal (§5).
2. **Simplified effective potential.** `cw1` has no two-loop terms and no $e^{4\Gamma}$, so the stability line
   is not at Buttazzo's precision. Agreement is 0.14 GeV with eq. (64); the comparison with eq. (66) is
   inconclusive.
3. **False vacuum at $\phi=0$.** The electroweak mass term is dropped. This is standard and harmless at these
   scales, but it is part of the model.
4. **Planck scale.** Frozen couplings decide 459 controlled `cw1` points (no finite bounce) and 78
   Planck-sensitive bounces. The near-critical analysis along $M_h=128$ GeV is impossible within the model.
5. **Deep thin wall.** $q\lesssim0.05$ is unresolved numerically. The thin-wall asymptotics ($S\propto q^{-3}$)
   are not reached, with $S/S_{\rm TW}=4.2$ at $q=0.05$. The leading-order $S_{\rm TW}$ used for this comparison
   is the reviewer-checked expression $0.0812/(\lambda''q^3)$.
6. **Matching domain.** Controlled statements apply only near the fit point. The exploratory scan is
   uncontrolled.
7. **Analytic results are formal.** The perturbative statements in §2.2 have no rigorous remainder bounds.
8. **Provenance relies on records outside the repository** for the exact scan build (§1.3).

## 8. Recommended next research directions

In priority order [I]:

1. **Correct the record.**
   * Attribute the 13/72 law, the Fubini average and the second-order term to Espinosa (2020).
   * Attribute the bound method to Wipf (1985), and cite Sato & Takimoto (2018) and Di Luzio et al. (2016).
   * Revise "only numerics are reliable" and the "every finite-action configuration" wording.
   * Label $S_{\rm curv}$ and $S_{\rm CGM}$ as resummed or variational forms of known results.
   * Add the unresummed first- and second-order formulas to the comparison.
2. **Study the empirical robustness** of the unresummed first-order formula over the whole $(M_h,M_t)$ plane and
   on other model families, e.g. cubic terms in $\ln\phi$. Seek an explanation.
3. **Develop and test a crossover approximation**, such as matched asymptotics or a tail-decoupled trial family,
   with error estimates. This is the most promising route to an original analytical result.
4. **Rigorous plateau criterion** via concentration–compactness.
5. **Connect to a one-loop computation**, to measure what fraction of SMVac's tree-level running correction is
   double-counted or modified in a consistent expansion.

## 9. Claims that are safe to make

* SMVac provides a validated numerical reference for the tree-level O(4) bounce of $V=\lambda_P(\phi)\phi^4/4$ in
  the SM.
  * Its numerical uncertainty is $\sim10^{-12}$ or better in the standard tier, supported by checks independent
    of solver agreement.
  * It reproduces the published second-order analytic expansion (Espinosa 2020) to about $10^{-3}$ in the
    coefficient.
* $S_{\rm FL}\le S$ for every finite-action bounce of this potential, strictly for running $\lambda$, by the
  Pohozaev identity and the sharp Sobolev inequality. This is an application of known methods (Wipf 1985).
* In the controlled region, $S_{\rm FL}$ underestimates the tree-level action by $3\times10^{-4}$ to
  $1.7\times10^{-3}$ for $S<5000$, in agreement with the known first-order law $\frac{13}{72q}$ (Espinosa 2020).
  It underestimates by up to 34% (`cw1`) and 50% (`tree`) near the stability boundary.
* Within the RG-improved tree-level scheme, the FL approximation shifts the metastability boundary by about
  0.0025 GeV. That is much less than the scheme ambiguities (0.09–0.12 GeV).
* Fixed-shape variational estimates for log-quadratic $\lambda$ are finite iff $q>\frac12{\rm Var}(\ln g)$. This
  gives the poles at $13/72$ (Fubini) and $3/32$ (the two-parameter family considered). It explains the observed
  failure of these estimates near criticality.
* The unresummed first-order formula $S_{\rm FL}(1+\frac{13}{72q})$ is empirically accurate to about 1.5% over the
  controlled scans. Present this as an observation, without a theoretical explanation.
* For `cw1` and $M_h\gtrsim127$ GeV near criticality, the frozen-coupling potential has no finite bounce. State
  this as numerical evidence about that model.

## 10. Claims that should not yet be made

* That the $13/72$ law, the Fubini average, the lower bound or the variational estimate are **new**.
* That SMVac **improves** or **corrects** Buttazzo et al. (2013), or that it shows an error in their lifetime.
* That closed-form approximations fail and **only numerics are reliable** for $0.05\lesssim q\lesssim1$.
* That $S\ge S_{\rm FL}$ holds for **every finite-action configuration**. It holds for critical points and for
  the mountain-pass barrier.
* That the no-finite-bounce classification is **proved**, or that Planck-plateau and Planck-sensitive results are
  **physical predictions**.
* That the finiteness criterion is a **significant new theorem**. At most it is a correct, possibly new,
  elementary lemma.
* That the complete scan was generated from **exactly** the tree of `ad63002`. Use the provenance wording in
  §1.3.
* Any lifetime statement at better than tree-level precision.

## 11. Complete references

1. D. Buttazzo, G. Degrassi, P. P. Giardino, G. F. Giudice, F. Sala, A. Salvio, A. Strumia, *Investigating the near-criticality of the Higgs boson*, JHEP 12 (2013) 089, [arXiv:1307.3536](https://arxiv.org/abs/1307.3536).
2. J. R. Espinosa, *Vacuum decay in the Standard Model: analytical results with running and gravity*, JCAP 06 (2020) 052, [doi:10.1088/1475-7516/2020/06/052](https://doi.org/10.1088/1475-7516/2020/06/052), [arXiv:2003.06219](https://arxiv.org/abs/2003.06219).
3. A. Wipf, *Upper and lower bounds for the bounce action*, J. Phys. A: Math. Gen. 18 (1985) 2521–2529, [author PDF](https://www.tpi.uni-jena.de/qfphysics/homepage/wipf/publications/papers/bounceaction.pdf).
4. R. Sato, M. Takimoto, *Absolute lower bound on the bounce action*, Phys. Rev. Lett. 120 (2018) 091802, [arXiv:1707.01099](https://arxiv.org/abs/1707.01099).
5. L. Di Luzio, G. Isidori, G. Ridolfi, *Stability of the electroweak ground state in the Standard Model and its extensions*, Phys. Lett. B 753 (2016) 150, [arXiv:1509.05028](https://arxiv.org/abs/1509.05028).
6. A. Andreassen, W. Frost, M. D. Schwartz, *Scale invariant instantons and the complete lifetime of the Standard Model*, Phys. Rev. D 97 (2018) 056006, [arXiv:1707.08124](https://arxiv.org/abs/1707.08124).
7. G. Isidori, G. Ridolfi, A. Strumia, *On the metastability of the Standard Model vacuum*, Nucl. Phys. B 609 (2001) 387, [hep-ph/0104016](https://arxiv.org/abs/hep-ph/0104016).
8. S. Chigusa, T. Moroi, Y. Shoji, *State-of-the-art calculation of the decay rate of electroweak vacuum in the Standard Model*, Phys. Rev. Lett. 119 (2017) 211801, [arXiv:1707.09301](https://arxiv.org/abs/1707.09301); *Decay rate of electroweak vacuum in the Standard Model and beyond*, Phys. Rev. D 97 (2018) 116012, [arXiv:1803.03902](https://arxiv.org/abs/1803.03902).
9. P. Baratella, M. Nemevšek, Y. Shoji, K. Trailović, L. Ubaldi, *Revising the full one-loop gauge prefactor in electroweak vacuum stability*, [arXiv:2406.05180](https://arxiv.org/abs/2406.05180).
10. G. Degrassi, S. Di Vita, J. Elias-Miró, J. R. Espinosa, G. F. Giudice, G. Isidori, A. Strumia, *Higgs mass and vacuum stability in the Standard Model at NNLO*, JHEP 08 (2012) 098, [arXiv:1205.6497](https://arxiv.org/abs/1205.6497).
11. I. Affleck, *On constrained instantons*, Nucl. Phys. B 191 (1981) 429–444, [doi:10.1016/0550-3213(81)90307-2](https://doi.org/10.1016/0550-3213(81)90307-2).
12. K. Lee, E. J. Weinberg, *Tunneling without barriers*, Nucl. Phys. B 267 (1986), [doi:10.1016/0550-3213(86)90150-1](https://doi.org/10.1016/0550-3213(86)90150-1).
13. J. R. Espinosa, *Tunneling without bounce*, Phys. Rev. D 100 (2019) 105002, [arXiv:1908.01730](https://arxiv.org/abs/1908.01730).
14. J. R. Espinosa, *A fresh look at the calculation of tunneling actions*, [arXiv:1805.03680](https://arxiv.org/abs/1805.03680).
15. J. R. Espinosa, T. Konstandin, *A fresh look at the calculation of tunneling actions in multi-field potentials*, JCAP 01 (2019) 051, [arXiv:1811.09185](https://arxiv.org/abs/1811.09185).
16. S. Coleman, V. Glaser, A. Martin, *Action minima among solutions to a class of Euclidean scalar field equations*, Commun. Math. Phys. 58 (1978) 211.
17. H. Brezis, L. Nirenberg, *Positive solutions of nonlinear elliptic equations involving critical Sobolev exponents*, Comm. Pure Appl. Math. 36 (1983) 437–477, [doi:10.1002/cpa.3160360405](https://doi.org/10.1002/cpa.3160360405).
18. T. Aubin, J. Differential Geom. 11 (1976) 573; G. Talenti, Ann. Mat. Pura Appl. 110 (1976) 353 (as cited in the repository; not re-verified here).
19. S. Fubini, Nuovo Cim. A 34 (1976) 521; L. N. Lipatov, Sov. Phys. JETP 45 (1977) 216 (as cited by refs. 2 and 5).
20. S. Coleman, Phys. Rev. D 15 (1977) 2929; C. G. Callan, S. Coleman, Phys. Rev. D 16 (1977) 1762.
21. P. Arnold, Phys. Rev. D 40 (1989) 613; P. Arnold, S. Vokos, Phys. Rev. D 44 (1991) 3620 (as cited by ref. 2).
22. G. Isidori, V. S. Rychkov, A. Strumia, N. Tetradis, *Gravitational corrections to Standard Model vacuum decay*, Phys. Rev. D 77 (2008) 025034, [arXiv:0712.0242](https://arxiv.org/abs/0712.0242).
23. A. Salvio, A. Strumia, N. Tetradis, A. Urbano, *On gravitational and thermal corrections to vacuum decay*, JHEP 09 (2016) 054, [arXiv:1608.02555](https://arxiv.org/abs/1608.02555).
24. V. Branchina, E. Messina, *Stability, Higgs boson mass and new physics*, Phys. Rev. Lett. 111 (2013) 241801, [arXiv:1307.5193](https://arxiv.org/abs/1307.5193).
25. A. Aravind, D. Lorshbough, S. Paban, *Lower bound for the multi-field bounce action*, Phys. Rev. D 89 (2014) 103535, [arXiv:1401.1230](https://arxiv.org/abs/1401.1230).
26. B. Garbrecht, P. Millington, *Fluctuations about the Fubini-Lipatov instanton for false vacuum decay in classically scale invariant models*, Phys. Rev. D 98 (2018) 016001, [arXiv:1804.04944](https://arxiv.org/abs/1804.04944).
27. C. L. Wainwright, *CosmoTransitions*, [arXiv:1109.4189](https://arxiv.org/abs/1109.4189); A. Masoumi, K. D. Olum, B. Shlaer, [arXiv:1610.06594](https://arxiv.org/abs/1610.06594); P. Athron et al., *BubbleProfiler*, Comput. Phys. Commun. 244 (2019) 448, [arXiv:1901.03714](https://arxiv.org/abs/1901.03714); V. Guada, M. Nemevšek, M. Pintar, *FindBounce*, [arXiv:2002.00881](https://arxiv.org/abs/2002.00881); S. Chigusa, T. Moroi, Y. Shoji, *Bounce configuration from gradient flow*, [arXiv:1906.10829](https://arxiv.org/abs/1906.10829); R. Sato, *Simple gradient flow equation for the bounce solution*, Phys. Rev. D 101 (2020) 016012, [arXiv:1907.02417](https://arxiv.org/abs/1907.02417).
28. L. N. Trefethen, *Spectral Methods in MATLAB*, SIAM (2000).

**Not verified in this review:**

* the full texts of refs. 7–9, which were checked at abstract level only;
* the bibliographic details of refs. 18–21, which are given as cited elsewhere;
* the attainment conditions of the Coleman–Glaser–Martin theorem for the frozen-coupling potential;
* whether the finiteness criterion appears in the mathematical literature on variational methods for critical
  elliptic problems.
