# AGENTS.md — Threshold Repository

## Quick Start

```bash
# All commands run from repo root (C:\Users\LENOVO\code\Threshold)

# Compile and run 1M-point scan (analytical or numerical solver)
g++ -std=c++17 -O3 -fopenmp solver_analytical.cpp -o solver_analytical.exe
python run_analytical.py                    # parallel chunks -> results/analytical_data.csv

g++ -std=c++17 -O3 -fopenmp solver_numerical.cpp -o solver_numerical.exe
python run_numerical.py                     # parallel chunks -> results/numerical_data.csv

# Plotting (reads CSV -> outputs PNG)
python plot_analytical.py                   # results/analytical_stability_plot.png
python plot_numerical.py                    # results/numerical_stability_plot.png
python plot_numerical_contours.py           # results/numerical_stability_contours.png
python plot_overlay.py                      # results/overlay_plot.png + overlay_closeup_plot.png
```

## Architecture

- **C++ solvers:** Custom 3-loop SM RGE integrator + bounce action evaluator. No external dependencies (no SimpleBounce, no GSL). Two variants:
  - `solver_analytical.cpp` — uses conformal approximation S = 8π²/3|λ_min|
  - `solver_numerical.cpp` — exact numerical bounce with RK4, Simpson quadrature, GSS radius optimization
- **Python plotting:** `matplotlib` + `pandas` + `scipy`. Each script self-contained. 300 dpi PNG output.
- **Paper:** `paper/main.tex` (revtex4-2, PRD format), `paper/references.bib`. Build with `pdflatex` + `bibtex` or `tectonic`.
- **Results:** CSV files cached; generated PNGs under `results/` and `paper/figures/`.

## Key Constraints

- **Do NOT modify `paper/main.tex`** unless explicitly asked. The paper is in a fragile draft state.
- **Do NOT delete or overwrite CSV data** under `results/`. Solver runs take hours.
- **All C++ solvers share RGE code** via copy-paste. Edits to `solver_numerical.cpp` RGE section must be ported to `solver_analytical.cpp`, `solver_numerical_dense.cpp`, `solver_analytical_dense.cpp`, `archive/solver_c6_sweep.cpp`.
- **Windows-only toolchain:** Mingw-w64 g++ from `external/w64devkit.exe` or system g++. OpenMP threading via `-fopenmp`.
- **Python venv** at `venv/`. Dependencies: `pandas`, `matplotlib`, `numpy`, `scipy`, `tqdm`.

## File Map

```
Threshold/
├── solver_*.cpp               # C++ bounce action solvers (root-level, production)
├── run_*.py                   # Python driver: parallel chunk execution + CSV merge
├── plot_*.py                  # Python plotting scripts (root-level, production)
├── paper/
│   ├── main.tex               # LaTeX paper (DO NOT MODIFY)
│   ├── main.pdf               # Compiled paper
│   ├── references.bib         # 28 entries
│   └── figures/               # *.png copies used by paper
├── results/                   # CSV data + generated PNGs
├── src/
│   ├── sm_stability/          # Older solvers, error scaling tests, builder scripts
│   ├── plotting/              # Older plot scripts (plot_phase_diagram, plot_scatter, etc.)
│   ├── cosmology/             # Higgs inflation, landscape (separate project)
│   └── build_and_run.py       # Archive: old compile/run pipeline (SimpleBounce-based)
├── archive/                   # Old compiled exes, C++ variants, tectonic
└── prompt_for_claude.md       # Analysis notes for LLM discussion of negative action physics
```

## Scientific Context

- **Paper claims:** (1) Conformal approx fails near metastability boundary (20–30% error). (2) c₆ φ⁶/M_Pl² operator rescues ~110 points (0.3–0.5 GeV in M_t). (3) Effect universal for c₆ > 0. (4) Fubini-Lipatov ansatz breaks down at low M_h → "sandwiched" unstable region via virial ratio S_kin/|S_pot| ≠ 2.
- **Phase classification:** 1=Stable, 2=Metastable, 3=Unstable, 4=Non-perturbative.
- **Adopted convention:** S_threshold(μ_inst) = 4 ln(T_U μ_inst), T_U = 4.3×10¹⁷ s, M_Pl = 1.22×10¹⁹ GeV.
- **Parameter ranges:** M_t ∈ [155, 185] GeV, M_h ∈ [115, 135] GeV, 500×500 grid (250k pts). Dense runs: M_t ∈ [155, 185] GeV, M_h ∈ [110, 140] GeV, finer spacing.
- **NNLO matching:** `get_nnlo_matching(Mh, Mt)` in solver files; includes threshold corrections δλ = -0.00313, δy_t = -0.0587.

## Accuracy Targets

- RK4 tolerance: 10⁻¹³, GSS tolerance: 10⁻¹³, Simpson N: 2048
- Total numerical error < 10⁻¹¹ relative action
- Physical uncertainty dominated by δM_t ~ 0.7 GeV (1σ)

## Plotting Conventions

- Stable=green, Metastable=yellow, Unstable=red, Non-perturbative=black
- Experimental point: star at (M_h=125.1, M_t=173.3) in blue with black edge
- Phase classification: 1=Stable, 2=Metastable, 3=Unstable, 4=Non-perturbative
- Standard DPI: 300, standard size: 10×8 or 11×9 inches
- Scripts load CSV from `results/`, save PNG to `results/`

## Script Execution Order

1. `run_analytical.py` → C++ analytical solver → `results/analytical_data.csv`
2. `run_numerical.py` → C++ numerical solver → `results/numerical_data.csv`
3. (Optional) `run_dense.py` → fine-grid solver → `results/*_dense_closeup.csv`
4. `plot_analytical.py` → reads analytical CSV → phase diagram PNGs
5. `plot_numerical.py` → reads numerical CSV → phase diagram PNGs
6. `plot_numerical_contours.py` → reads numerical CSV → contour overlay PNG
7. `plot_overlay.py` → reads both CSVs → disagreement overlay PNGs

## Avoiding Pitfalls

- **CSV column naming:** Analytical solver writes `Mh_calc`; numerical solver writes `Mh_calc`. Overlay script rounds to 6 decimals and uses `cKDTree` for alignment — mismatches cause empty disagreement plots.
- **Point count confusion:** `TOTAL_POINTS = 1002001` in run scripts is NOT 250k. The 250k figure in the paper corresponds to the (500×500) grid; the 1M-point runs are supersets used for smooth contouring.
- **c₆ sweep:** `archive/solver_c6_sweep.cpp` is a separate solver not wired into the main pipeline. It outputs four-column CSV with c₆ as a parameter.
- **Negative actions:** Points with S ≤ 0 are classified as unstable (the Fubini ansatz has broken down, not the physics). Do not treat negative S as numerical error unless |S| > 10⁶.
- **Non-perturbative region:** Points where Yukawa or gauge couplings blow up before M_Pl are flagged (Stability=4). These appear as black dots at large M_t / M_h.
