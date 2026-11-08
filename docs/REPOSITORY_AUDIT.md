# Repository Audit (Stage A)

**Project:** SMVac — Naman Goyal (namangoyal31), Indian Institute of Science (IISc), Bangalore
**Audited revision:** `main` @ `a7c070ae501b0669b37e6a583a0a18da575e6fe1` ("Initial commit: High-precision
electroweak vacuum stability analysis"). This is the only commit; there is no earlier history.
**Working branch:** `formalization`. `main` is untouched; nothing has been pushed.
**Audit date:** 2026-11-08.
**Environment used:** Windows 11; MinGW.org GCC 6.3.0 (32-bit); Python 3.14, numpy 2.5,
pandas 3.0, scipy 1.18, matplotlib 3.11; CMake 4.4.

Scope: every source file, header, script, build/CI file, dataset, reference file, figure and
document in the repository was inspected. The execution flow was traced from the code. The
physics and numerics are reviewed in `docs/SCIENTIFIC_METHODOLOGY.md` and
`docs/NUMERICAL_METHODS_REVIEW.md`. All quantitative statements here can be reproduced with the
scripts in `audit/` (see `audit/README.md`).

---

## 1. Repository structure (95 tracked files)

```
apps/solver_numerical.cpp        488-line monolithic solver: matching, RGEs, lambda_eff, "S_exact", grid driver
scripts/                         run_numerical.py, run_analytical.py (parallel chunk drivers),
                                 plot_*.py (4 plotting scripts), plot_traces.py (orphan)
data/numerical_data.csv          49.6 MB, 1,737,982 rows (see section 5)
data/analytical_data.csv         18.4 MB, 1,002,001 rows
figures/*.png                    8 figures
SMVacuumDecay/                   modular C++ re-packaging of the same code
  include/SMVacuumDecay/*.hpp    RGE, EffectivePotential, CanonicalBounce, FubiniLipatov, numerics/{RK4,Simpson,GoldenSectionSearch}
  src/*.cpp                      the same physics as apps/solver_numerical.cpp (verified identical results)
  apps/                          generate_phase_diagram.cpp (adds --analytical/--dense modes), benchmark_scan.cpp
  tests/                         5 executables (regression tests, section 6)
  reference/v1.0/                frozen outputs: JSON benchmarks, rge/potential/phase_scan/bounce CSVs
  python/                        near-duplicates of scripts/ with results/ paths, plus run_dense.py
  CMakeLists.txt, LICENSE        EMPTY files (0 bytes)
  CITATION.cff                   placeholder URL github.com/YourUsername/SMVacuumDecay
docs/                            29 files, ~18,000 lines (section 7)
LICENSES/                        Bertini and GMP/MPFR licences (software not used by this repository)
README.md, CITATION.cff, CHANGELOG.md, CONTRIBUTING.md, requirements.txt, .github/workflows/ci.yml
```

**Duplication.** The physics exists twice: `apps/solver_numerical.cpp` and `SMVacuumDecay/src`.
The plotting and run scripts also exist twice, differing only in paths and titles.

## 2. Execution flow (as implemented)

1. `scripts/run_numerical.py` splits 1,002,001 grid points ($M_t,M_h\in[0,250]$ GeV, step 0.25)
   into 12 chunks. It launches `apps/generate_phase_diagram.exe --numerical s e` per chunk, which
   writes `data/data_numerical_chunk_<s>.csv`, then globs *all* chunk files, concatenates them into
   `data/numerical_data.csv`, and deletes the chunks.
2. For each point, `classify_stability(Mh, Mt)`:
   * matching at $\mu=M_t$ (linear fits);
   * RK4 RG evolution to $1.22\times10^{19}$ GeV in steps Δln μ² = 0.1, tracking
     $\lambda_{\rm eff}$ on nodes above $\max(M_t,v)$;
   * status 4 if $|\lambda|$ or $y_t$ exceeds $4\pi$; status 1 if $\lambda_{\rm eff}\ge0$ everywhere;
   * otherwise `S_approx` $=8\pi^2/(3|\lambda_{\rm eff,min}|)$ and `S_exact` = golden-section minimum
     of `evaluate_action_at_R`;
   * status 2 if `S_exact` $>4\ln(T\mu_{\rm inst})$, else 3.
3. `scripts/plot_numerical_contours.py` reads the CSV, replaces `S_exact<0` by $10^8$ and
   `S_exact>10^{10}` by 0, and draws `tricontour` lines at
   S = 10, 50, 100, 200, 450, 1000, 3000, 500000 → `figures/phase_diagram_contours.png`.
4. **`--analytical` in `apps/solver_numerical.cpp` changes only the output file name.** The
   computation is identical. `data/analytical_data.csv` (4 columns) was produced by a different
   program: `classify_buttazzo` in the `SMVacuumDecay` version, or its predecessor, which uses
   `S_approx` for the classification. Its `S_approx` column equals the numerical file's
   `S_approx` exactly at all 253,412 common unstable/metastable points.
5. The README smoke test says the analytical chunk goes to `results/`, but the code writes to
   `data/`.

## 3. Scientific objective and framework (summary)

The code classifies $(M_h,M_t)$ as stable, metastable, unstable or non-perturbative. It uses
three-loop SM running from $M_t$ to $M_{\rm Pl}$, a one-loop large-field effective quartic, and
a tree-level Fubini–Lipatov estimate of the tunnelling action. The headline claim is that an
"RG-improved Fubini–Lipatov action" (`S_exact`) deviates from the conformal estimate
(`S_approx`) near phase boundaries. The full reconstruction, with equation-by-equation
verification, is in `docs/SCIENTIFIC_METHODOLOGY.md`.

## 4. Build, dependencies, CI

* **No build system.** `SMVacuumDecay/CMakeLists.txt` is empty. The README compile line uses
  `-fopenmp`, which fails to link with MinGW.org GCC 6.3 (`cannot find -lpthread`). There are no
  OpenMP pragmas anywhere, so the flag is inert. Without it, both trees compile:
  `apps/solver_numerical.cpp` gives 13 `-Wextra` warnings, and the library tests compile cleanly
  apart from 2 warnings in `test_FL.cpp`.
* **Python:** `requirements.txt` lists pandas, numpy and matplotlib. `plot_overlay.py` also needs
  **scipy**, which is missing. No versions are pinned.
* **CI** (`.github/workflows/ci.yml`) compiles the monolith and runs two 10-point smoke tests at
  $M_t=0$, which are all status 4. It never builds or runs the `SMVacuumDecay` tests. Per the
  owner, GitHub Actions is currently unavailable.
* External libraries: none. Licences for Bertini and GMP/MPFR are present, but neither is used.

## 5. Existing results and their provenance

### 5.1 Datasets

| File | Facts (verified with `audit/` scripts) | Reproducible? |
|---|---|---|
| `data/numerical_data.csv` | 1,737,982 data rows, but only **1,002,001 unique grid points**. **735,980 rows are exact duplicates**, an artefact of `run_numerical.py` globbing chunk files from several runs. One **corrupted row** at line 1,320,389 (`2.0,472.97,472.932,,`, a truncated copy of `166.75,85.0,2.0,472.97,472.932`). The `Stability` column is float because of that row. | **Yes**, for the unique points. 390 stratified random points were recomputed with the current `apps/solver_numerical.cpp`. All 390 classifications match. The actions match to the 6 stored significant digits, except one pathological point (S ≈ −3.98) that differs at $1.5\times10^{-4}$ because of floating-point sensitivity. |
| `data/analytical_data.csv` | 1,002,001 unique rows, no duplicates | Its generator (`classify_buttazzo`) exists only in `SMVacuumDecay/`, and no script in the repo writes this file name from the top-level tree. The values are consistent with that code. |
| `data/*dense_closeup.csv` | **missing** | Required by `plot_numerical.py`, `plot_analytical.py` and `plot_overlay.py`. |

Content issues in the numerical data (unique points):

* 8,740 points have **negative** `S_exact`, as low as $-8.1\times10^8$; all are classified unstable.
* 217 points carry the sentinel $10^{100}$; **all are classified metastable**.
* 617 points have `S_exact` $>5\times10^5$, from the kinetic-only early exit.
* Inside $M_h\in[110,140]$, $M_t\in[165,185]$, `S_exact/S_approx − 1` ranges from −84% to
  $+5\times10^6$.
* Numerical and analytical classifications differ at 7,648 points.

README validation numbers: −0.01% at (115, 180) is actually −0.0151%, and 32.1% at
(134.75, 176.5) is reproduced. Both are explained in METHODOLOGY §6.6.

### 5.2 Figures

| Figure | Generator in repo | Reproducible | Notes |
|---|---|---|---|
| `phase_diagram_contours.png` (README "Primary Result") | `plot_numerical_contours.py` (its title text differs: "FL …" vs "Numerical …" in the PNG) | partially | contour at 500000 sits on a code artefact; negative actions → $10^8$, sentinel → 0 before contouring; overlapping labels ("S=1S0=100"); duplicated points plotted; linear-fit extrapolation to $M_t\to0$ |
| `overlay.png` | none (script writes `overlay_plot.png`) | no (name mismatch) | blue "disagreement" band at $M_h<60$, $M_t\approx98$–115 is the negative-`S_exact` artefact region |
| `overlay_closeup_plot.png` | `plot_overlay.py` | **no** (dense data missing) | |
| `FL_stability_plot.png` | `plot_numerical.py` | yes | titled "FL" although it shows the `S_exact` classification |
| `FL_stability_closeup.png`, `analytical_stability_closeup.png` | `plot_numerical.py`, `plot_analytical.py` | **no** (dense data missing) | |
| `analytical_stability_plot.png` | `plot_analytical.py` | yes (from the committed CSV) | |
| `convergence.png` | **none** | **no** | README calls it an "optimization convergence" diagnostic. It actually shows a *shooting* bounce for $(M_h,M_t)=(5,105)$ from a solver absent from this repo. It does show that a finite bounce exists where `S_exact` is negative. |

### 5.3 Reference files (`SMVacuumDecay/reference/v1.0`)

These are outputs of the legacy code, frozen to 16 digits, including the negative action
`Numerical_action = -1364436.824475072` at (5, 105). `bounce_reference.csv` (a profile with
$\phi(0)=23725.67$ GeV) and `phase_scan_reference.csv` are used by no code.

## 6. Testing status

The five `SMVacuumDecay/tests` executables were built manually (no CMake) and run from the
repository root. **All five pass.** All are regression tests against frozen outputs of the same
code, with tolerances of $10^{-12}$ relative or $10^{-9}$ absolute:

| Test | What it asserts | Problem |
|---|---|---|
| `test_rge` | couplings vs `rge_reference.csv` | the reference inherits the RGE bugs (METHODOLOGY §3) |
| `test_potential` | λ_eff, V vs reference at (5, 105) | inherits the Z-coupling bug; this point is outside the validity domain |
| `test_bounce` | `S_exact(5,105) = −1364436.824475072` | **enshrines a physically impossible negative action** |
| `test_FL` | FL action and status at two points | regression only |
| `test_scientific` | statuses and actions at (125.1, 173.1) | regression only; no independent reference |

None of the tests checks an analytic limit, a literature value, a convergence property, or a
failure mode. `S_exact` is platform-sensitive at the $10^{-7}$ level (x87 vs SSE2), so the
$10^{-12}$ tolerances pass only on the reference platform and points.

## 7. Documentation status

* **`README.md`** contains claims contradicted by the code:
  * OpenMP parallelisation;
  * "NNLO matching" (actually linear fits, with a wrong $g_2$);
  * the validation numbers (method artefacts);
  * the convergence figure description;
  * "diagnostic scripts under `src/`" (no such directory);
  * no mention of the $\phi^6$ operator.
* **`SMVacuumDecay/README.md`** claims an adaptive RK4 (unused), "exact evaluation of the
  Euclidean action" (false), and tests that "match legacy literature" (they match legacy
  *outputs*). It also says CMake is pending.
* **`SMVacuumDecay/CODE_PROVENANCE.md`** claims "100% migrated … byte-for-byte". This is
  consistent with our reproduction, but it refers to a predecessor repository ("Threshold",
  commit 7c9eeff) that is not part of this repository.
* **`docs/`** (29 files) is planning and review material from the predecessor project:
  * paper planning and figure planning;
  * mock referee audits and a rebuttal;
  * error-analysis claims;
  * three prompts written for an AI assistant (`prompt_for_claude.md`,
    `critique_prompt_for_claude.md`, `reformulate_prompt_for_claude.md`) and an agent
    configuration (`AGENTS.md`);
  * a 12,319-line file inventory of the predecessor repository;
  * `docs/README`, which is the README of the unrelated *Bertini* package.
* Several docs already diagnose parts of the problem: `SOLVER_VERIFICATION.md` states that no
  bounce ODE is solved, and `SANDWICH_ANALYSIS.md` states that the negative actions are an ansatz
  breakdown. These conclusions were never propagated to the code, tests, README or figures.
* **Licence:** the repository has **no licence**; `SMVacuumDecay/LICENSE` is empty. No licence
  should be added without the owner's decision.
* **Citation metadata:** the top-level `CITATION.cff` is reasonable.
  `SMVacuumDecay/CITATION.cff` has a placeholder URL.

## 8. Related repositories

`SMVac-Research` is **not present** in the parent directory, so no comparison was made. The
predecessor "Threshold" repository named in `CODE_PROVENANCE.md` exists locally but was **not
inspected** in this audit. All conclusions here derive from this repository alone.

## 9. Known bugs and weaknesses (ranked)

See METHODOLOGY §0 and NUMERICAL REVIEW §9 for details and evidence.

1. `S_exact` is not a bounce action. It is unsound, with negative values and sentinels classified
   as metastable. **This is fundamental to the headline result.**
2. A $\phi^6/(6M_{\rm Pl}^2)$ operator is added in the `S_exact` integrand only. It is not
   documented in the README or code; the predecessor project's notes (`docs/AGENTS.md`) list it
   as a deliberate paper topic (the "c₆ operator"). That file also gives
   $T_U=4.3\times10^{17}$ s, whereas the code uses $3.15\times10^{17}$ s.
3. The metastability threshold uses the wrong scale and time (45–65 units too low).
4. $g_2(M_t)$ is wrong in the matching (0.65355 vs 0.64779); the $g_1$ and $y_t$ slopes are off.
5. RGE transcription errors: the 4-loop sign, and two 3-loop terms each missing a factor $y_t^2$.
6. Z-boson coupling normalisation in $\lambda_{\rm eff}$.
7. Linear matching fits are extrapolated to $M_t,M_h\in[0,250]$ GeV ($y_t<0$, $\ln 0$).
8. Data integrity: duplicates, a corrupted row, missing datasets, figures that cannot be regenerated.
9. Plot-level data manipulation; a contour placed on an artefact.
10. Tests are regression-only and enshrine wrong values; there is no build system; CI does not
    run the tests.
11. Duplicated code trees; magic numbers; `CLAUDE_M_PLANCK` naming; failure modes are silent.
12. Documentation over-claims; AI-prompt and third-party files in `docs/`; no licence.
