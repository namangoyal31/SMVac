# Minimal SMVac Migration Plan

This plan aggressively distills the Threshold codebase into the **smallest possible** repository required to fully reproduce the scientific results and figures presented in the paper. 

Redundant scripts, intermediate debugging utilities, standalone test cases, and obsolete experiments have been completely removed. Shared logic is recommended for consolidation.

## 1. C++ Source Code Consolidation

Currently, the C++ code suffers from heavy copy-pasting (especially the 3-loop RGE running and potential evaluation). To minimize files and maintain DRY (Don't Repeat Yourself) principles:

### A. Core Solvers -> `src/`

| Current Files | Recommended New File | Explanation |
|---------------|----------------------|-------------|
| `solver_analytical.cpp`, `solver_analytical_dense.cpp`, `solver_numerical.cpp`, `solver_numerical_dense.cpp` | **`src/scan_solver.cpp`** | **MERGE**: These four solvers differ only by the integration method (conformal approx vs. exact RK4) and the grid density. Merge into a single executable that accepts `--mode analytical|numerical` and `--grid standard|dense`. |
| `stationary_finder.cpp`, `verify_potential.cpp`, `potential_analysis_rge.cpp`, `potential_scan.cpp` | **`src/potential_analyzer.cpp`** | **MERGE**: These all extract the effective potential and its derivatives, locating barriers and minima. Merge into a single diagnostic utility for 1D potential profiling and stationary point extraction. |
| `run_bisection.cpp`, `sweep_interval.cpp` | **`src/canonical_bounce.cpp`** | **MERGE**: Both handle the exact canonical ODE shooting. `sweep_interval` just scans an array, while `run_bisection` hones in. Combine into one tool that can either `--sweep` an interval or `--bisect` to exact convergence. |

### B. Shared Library -> `include/`

| Current Files | Recommended New File | Explanation |
|---------------|----------------------|-------------|
| *Embedded RGEs in all `*.cpp` files* | **`include/StandardModelRGE.hpp`** | **EXTRACT**: The 100+ lines of 3-loop $\beta$-functions and NNLO matching conditions are currently copy-pasted across 15+ files. Extract into a single header included by all solvers. |

## 2. Python Execution & Plotting Consolidation

### A. Execution Drivers -> `scripts/runners/`

| Current Files | Recommended New File | Explanation |
|---------------|----------------------|-------------|
| `run_analytical.py`, `run_numerical.py`, `run_dense.py`, `run_ode_grid.py` | **`scripts/run_grid.py`** | **MERGE**: These scripts all chunk a 500x500 grid over $M_h$ and $M_t$ and invoke a C++ solver in parallel. A single `run_grid.py` with arguments (e.g. `--solver numerical --dense`) is sufficient. |

### B. Plotting Scripts -> `scripts/plotting/`

| Current Files | Recommended New File | Explanation |
|---------------|----------------------|-------------|
| `plot_analytical.py`, `plot_numerical.py`, `plot_numerical_contours.py`, `plot_overlay.py` | **`scripts/plot_phase_diagram.py`** | **MERGE**: These all parse CSV grids and produce $M_h$ vs $M_t$ phase diagrams. A single script can generate the individual plots and the disagreement overlay. |
| `plot_bisection.py`, `plot_traces.py` | **`scripts/plot_bounce.py`** | **MERGE**: Both plot $\phi(r)$ vs $r$. Merge into one script handling bounce profiles and over/undershoot traces. |

## 3. Exclusions (Do Not Migrate)

The following groups of files are strictly excluded from the final repository:

| Category | Representative Files | Reason for Exclusion |
|----------|----------------------|----------------------|
| **Debugging / Intermediate Solvers** | `solver_canonical.cpp`, `solver_instrumented.cpp`, `solver_diagnostic.cpp`, `solver_maxsteps_diag.cpp`, `deriv_solver.cpp`, `validation_solver.cpp` | Obsolete. The successful $\epsilon$-scaled $r_0$ fix makes all prior diagnostic solvers irrelevant. |
| **Stress Tests** | `stress_solver.cpp`, `test_nan.cpp`, `test_v.cpp`, `test_dv.cpp` | One-off internal validation tests; not required for reproducing the paper. |
| **Archival Prototyping** | Entire `src/sm_stability/`, `archive/`, `src/cosmology/` | Legacy code superseded by the main solvers. |
| **Python Diagnostics** | `_sandwich_*.py`, `_audit_*.py`, `_check_chunk*.py`, `verify_continuity.py`, `analyze_maxsteps.py` | Disposable ad-hoc scripts used during the interactive sessions to verify data integrity. |
| **AI Work Logs** | `*.md` (except `README.md` and finalized docs) | The AI reasoning logs (e.g. `POTENTIAL_STRUCTURE_REPORT.md`, `STORY_ANALYSIS.md`) are conversational artifacts, not scientific documentation. |

## 4. Final Directory Structure

The resulting repository is extremely lean, modular, and directly maps to the publication results.

```
SMVac/
├── README.md                           # Compilation instructions & reproducibility workflow
├── LICENSE
├── include/
│   └── StandardModelRGE.hpp            # Extracted single-source-of-truth for RGEs
├── src/
│   ├── scan_solver.cpp                 # Phase diagram grid generator
│   ├── potential_analyzer.cpp          # Stationary point / sandwich structure
│   └── canonical_bounce.cpp            # Exact exact ODE bisection & traces
├── scripts/
│   ├── run_grid.py                     # Parallel cluster execution driver
│   ├── plot_phase_diagram.py           # Reproduces all phase boundary / overlay figures
│   └── plot_bounce.py                  # Reproduces the canonical bounce traces
├── data/                               # CSVs (tracked via git-lfs or generated via scripts)
├── paper/                              # LaTeX source
│   ├── main.tex
│   ├── references.bib
│   └── figures/                        # Generated PNGs imported by LaTeX
```

## Summary

This plan reduces the repository from **150+ source files** to exactly **1 header, 3 C++ files, and 3 Python scripts**—without losing a single piece of scientific functionality. Every paper figure can be regenerated by running `run_grid.py` followed by the two plotting scripts.
