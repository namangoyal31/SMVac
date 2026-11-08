# Audit tooling

Scripts that reproduce every quantitative statement in `docs/REPOSITORY_AUDIT.md`,
`docs/SCIENTIFIC_METHODOLOGY.md` and `docs/NUMERICAL_METHODS_REVIEW.md`. None of them modifies
the original implementation. Run them from the repository root with Python ≥ 3.10, numpy, scipy
and matplotlib.

| Script | Purpose |
|---|---|
| `legacy_port.py` | Faithful Python port of `apps/solver_numerical.cpp`. The beta functions are extracted mechanically into `_legacy_betas_extracted.py`. Options switch individual defects on or off. |
| `validate_port.py` | Compares the port with the compiled legacy binary (`legacy_driver.cpp`). It agrees exactly with an SSE2 build; x87 builds differ at ~1e-6 in `S_exact` only. |
| `legacy_driver.cpp` | Evaluates the legacy `classify_stability` at points read from stdin. |
| `check_power_counting.py` | Literature-independent loop-order consistency test of the beta functions. |
| `independent.py` | Independent reference numerics: DOP853 RG solution, rigorous lower bound, CGM variational upper estimate, numerical O(4) shooting bounce. |
| `test_independent.py` | Analytic-limit tests of `independent.py`. |
| `exp_action_landscape.py`, `exp_action_landscape_zoom.py` | Landscape and roughness of the legacy `S_exact(R)` objective. |
| `exp_bounce_comparison.py` | Legacy actions vs bounds vs numerical bounce at benchmark points. |
| `exp_defect_impact.py` | Effect of each physics defect; comparison with Buttazzo et al. 2013 eq. (61); RK4 convergence. |

Building the legacy driver (64-bit IEEE arithmetic, matching the port):

```bash
g++ -std=c++17 -O3 -msse2 -mfpmath=sse -Dmain=legacy_main -c apps/solver_numerical.cpp -o solver.o
g++ -std=c++17 -O3 -msse2 -mfpmath=sse -c audit/legacy_driver.cpp -o driver.o
g++ solver.o driver.o -o legacy_driver
printf "125.1 173.1\n5 105\n" | ./legacy_driver > out.csv
python -m audit.validate_port out.csv
```

Outputs are written to `audit/results/`.
