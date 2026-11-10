# Contributing to SMVac

## Development setup

* C++14 compiler, CMake ≥ 3.10 (Ninja or Make).
* Python ≥ 3.11 with the packages in `requirements.txt`.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
ctest --test-dir build --output-on-failure
python -m unittest discover -s python/tests
```

## Ground rules for scientific changes

* **Separate physics from numerics.**
  * Physics inputs live in `src/matching.cpp`, `src/beta_functions.cpp` and `src/potential.cpp`.
  * Algorithms live in `src/bounce.cpp`, `src/rge.cpp` and `include/smvac/numerics/`.
  * Change one at a time, and document the reason in `CHANGELOG.md`.
* **No silent physics.** Anything beyond the Standard Model (e.g. the $\phi^6$ deformation) must be opt-in,
  recorded in the output metadata, and excluded from SM figures.
* **Tests need independent references.** Use analytic limits, literature values (with equation numbers),
  convergence orders, or values from the independent Python implementation. Never freeze an output of the code
  under test as its own reference.
* **Report uncertainties honestly.** A numerical result is a numerical solution with an estimated uncertainty, not
  an exact value. Call a bound *rigorous* only when its assumptions are proven.
* **Do not loosen tolerances to make a test pass.** Investigate the failure and document the cause.
* **Regenerate, do not edit, results.** Figures and CSVs in `results/` and `figures/` come from the scripts in
  `validation/` and `python/smvac/`. Large intermediate outputs are not committed.

## Pull requests

Describe the scientific motivation, the validation performed (tests, comparisons, convergence), and the effect on
the published results (`results/summary.json`).
