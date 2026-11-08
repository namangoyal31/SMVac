# Legacy implementation (frozen)

This directory holds the original SMVac implementation exactly as it was at commit `a7c070a`. The
files were moved here with `git mv`; their content is unchanged, and their history is available
through `git log --follow`. They are kept for **traceability and comparison only**. They are not
the scientific reference of this repository.

The audit (`docs/REPOSITORY_AUDIT.md`, `docs/SCIENTIFIC_METHODOLOGY.md`) found that:

* the quantity called `S_exact` is a first-order perturbative estimate around the Fubini–Lipatov
  instanton, not a bounce action. It produces negative values, sentinels classified as
  metastable, and values below a rigorous lower bound;
* the `S_exact` integrand contains a $\phi^6/(6M_{\rm Pl}^2)$ operator that is not Standard Model physics;
* the metastability threshold uses the instability scale and $T\simeq10$ Gyr;
* $g_2(M_t)$ is 0.65355 instead of the published 0.64779;
* there are three RGE transcription errors (the 4-loop QCD sign and two 3-loop terms each missing
  $y_t^2$);
* the Z-boson coupling in $\lambda_{\rm eff}$ uses the GUT-normalised $g_1^2$.

## Contents

| Path | What it is |
|---|---|
| `apps/solver_numerical.cpp` | monolithic legacy solver (wide scan, `S_exact`, `S_approx`) |
| `SMVacuumDecay/` | modular re-packaging of the same code, with regression tests and frozen references (`reference/v1.0`); its `CMakeLists.txt` and `LICENSE` are empty files, as in the original |
| `scripts/` | legacy run and plotting scripts (they expect paths relative to the old layout) |
| `figures/` | the eight figures that were committed originally; provenance is in REPOSITORY_AUDIT §5.2 |
| `data/numerical_data.csv` | historical output of the legacy wide scan (see below) |
| `data/analytical_data.csv` | historical output of `classify_buttazzo` (FL classification) |

## Data dictionary (historical outputs; exploratory, not physical predictions)

Grid: $M_t, M_h\in[0,250]$ GeV, step 0.25 GeV (1,002,001 points). This range lies far outside
the validity of the matching formulas.

**`numerical_data.csv`** — columns:

* `Mt` [GeV];
* `Mh_calc` [GeV]: the input $M_h$; nothing is calculated;
* `Stability` (1 stable, 2 metastable, 3 unstable, 4 non-perturbative);
* `S_exact` (legacy construction; −1 when not computed);
* `S_approx` $=8\pi^2/(3|\lambda_{\rm eff,min}|)$, with the legacy physics.

Known integrity issues:

* 735,980 exact duplicate rows, from the chunk-merging script;
* one truncated row at line 1,320,389;
* `Stability` stored as float.

The unique content (1,002,001 points) is reproducible with `apps/solver_numerical.cpp`; 390
sampled points were verified. To load a clean copy:

```python
import pandas as pd
df = pd.read_csv("legacy/data/numerical_data.csv")
df = df[df.Mh_calc <= 250].drop_duplicates(["Mt", "Mh_calc"])
```

**`analytical_data.csv`** — columns `Mt`, `Mh_calc`, `Stability`, `S_approx`; 1,002,001 unique rows.

## Rebuilding the legacy solver

```bash
g++ -std=c++17 -O3 legacy/apps/solver_numerical.cpp -o legacy_solver   # -fopenmp is not needed (no OpenMP code)
```

The legacy binary writes `data/data_numerical_chunk_<start>.csv` relative to the current
directory. For point-wise evaluation, use `audit/legacy_driver.cpp` (see `audit/README.md`).
