# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0-dev] — formalization branch (unreleased)

### Stage D/E — near-critical study, final validation and documentation
- **Near-critical study** (`validation/near_critical.py`, `validation/analyze_near_critical.py`):
  - parabolic model with $q$ from 1000 to 0.05;
  - SM lines at M_h = 120, 125.15 and 128 GeV (`cw1`) and 125.15 GeV (`tree`), with M_t - M_t^crit from 1e-4
    to 3 GeV;
  - outputs: error thresholds, crossovers, thick-wall coefficients, thin-wall comparison, q-delta map, universality.
- **Variational failure explained (derived and verified):**
  - A fixed-shape trial estimate is finite iff q > (1/2) Var(ln g). This gives the pole at 13/72 for the FL trial.
  - The two-parameter estimate S_CGM2 (`python/smvac/reference/analytic.py`) moves the pole to 3/32 but does not
    remove it.
  - The exact bounce obeys S = (pi^2/2) lambda'' M (q - (1/2) Var_b(ln phi)), with <ln phi>_b = x_min.
  - Checks 6–8 are in `validation/verify_derivations.py`.
- **Fixed:** overflow in the two-parameter trial profile for large n, which made S_CGM2 spuriously infinite.
  A grid pre-scan was added before the Nelder–Mead refinement.
- **Solver failures resolved** (`validation/fallback_shooting.py`). All 360 `cw1` collocation failures are
  Planck-sensitive:
  - 78 are metastable, solved by M1/M2 shooting (agreement 5e-10);
  - 282 have no finite bounce in the frozen-coupling potential (Planck plateau), classified by explicit criteria
    that distinguish them from solver failures;
  - the plateau points are metastable by the rigorous bound S >= S_FL > 559.
- **Validation report** completed: scan-level checks (§5) and a separate uncertainty budget (§8): numerical,
  approximation, parametric, theoretical, Planck-model.
- **Figures** regenerated from the final data, with captions in `figures/README.md`. Added `near_critical.png`
  and `variational_failure.png`. Fixed the legacy panel, which now shows the divergence of the legacy objective.
- **Reproducibility:**
  - `validation/reproduce_sample.py` recomputes a seeded, stratified sample of scan, fallback and near-critical
    model points;
  - `validation/check_links.py` checks the documentation references;
  - committed scan outputs are gzip-compressed.

### Stage C — new implementation and validation
- **Scientific objective restated:** quantify the accuracy of the analytical Fubini–Lipatov approximation against
  a validated numerical reference for the tree-level bounce of the RG-improved SM potential (docs/STAGE_C_PLAN.md).
- **Numerical reference:** Chebyshev spectral collocation (adaptive resolution, Derrick/dilatation identity checks,
  two-resolution uncertainty). Selected over shooting (two formulations), adaptive low-order collocation and the
  legacy construction by a systematic comparison (docs/NUMERICAL_METHOD_COMPARISON.md).
- **Analytical estimates:** S_FL (rigorous lower bound), S_curv (curvature-corrected, coefficient 13/72) and
  S_CGM (variational upper estimate).
- **Physics corrections** (each validated separately in `audit/exp_defect_impact.py`):
  - Matching follows Buttazzo et al. (2013) eqs. 55–60. The legacy value g2(M_t) = 0.65355 is replaced by 0.64779,
    and the g1 and y_t slopes are corrected.
  - RGEs are transcribed afresh from their App. B. This fixes the sign of the 4-loop QCD term and the y_t^2
    factors missing in the 3-loop beta_lambda and beta_yt terms.
  - The Z-boson term of lambda_eff uses g2^2 + gY^2.
  - The metastability criterion follows Buttazzo et al. eqs. 80–82, with Lambda_B = 1/R of the bounce and
    H0 = 67.4 km/s/Mpc. The legacy criterion used the instability scale and T = 10 Gyr.
  - The phi^6/M_Pl^2 operator is removed from the SM calculation and kept only as an opt-in, flagged deformation.
- **Domain of validity:**
  - Headline results are restricted to M_h in [120,130] and M_t in [166,182] GeV.
  - The wide [0,250]^2 scan is recomputed as a labelled exploratory extrapolation.
  - Points whose bounce probes phi > M_Pl are flagged Planck-sensitive.
- **Software:** C++14 library with CMake and CTest (7 test suites), plus an independent Python reference
  implementation with tests. Scans are config-driven and parallel, and their outputs carry metadata headers. Every
  figure is regenerated from data.
- **Repository:**
  - The original implementation is preserved unchanged in `legacy/`.
  - AI-assistant prompt files, the AI-agent configuration and the unrelated Bertini README are removed.
  - Predecessor documents are archived in `docs/archive/`.
  - The licence is unchanged; the owner's decision is pending.

### Stage A/B — review (no change to the scientific implementation)
- `docs/REPOSITORY_AUDIT.md`, `docs/SCIENTIFIC_METHODOLOGY.md`, `docs/NUMERICAL_METHODS_REVIEW.md`:
  independent audit of the original implementation (physics, numerics, data, figures, tests).
- `audit/`: reproducible audit tooling, comprising:
  - a validated Python port of the legacy solver;
  - a loop-order power-counting check;
  - an independent RG solver;
  - the lower bound and variational estimate;
  - a shooting bounce;
  - defect-impact and landscape experiments.

## [1.0.0] — original implementation (commit a7c070a)
- Monolithic C++ solver (`legacy/apps/solver_numerical.cpp`) and its modular copy (`legacy/SMVacuumDecay/`).
- Python run/plot scripts, 1,002,001-point wide scan, eight figures.
