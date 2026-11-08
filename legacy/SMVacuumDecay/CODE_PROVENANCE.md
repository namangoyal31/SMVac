# Code Provenance & Inventory

This document traces every function from the canonical legacy implementations (`solver_numerical.cpp` and `solver_analytical.cpp`) to their destination in the new `SMVacuumDecay` architecture. It ensures no physics logic is lost or accidentally modified during the migration.

Line numbers refer to the frozen Threshold repository prior to migration (commit `7c9eeff` or local equivalent).

## 1. `solver_numerical.cpp`

| Function / Class | Legacy Lines | Classification | New File Destination | Validation Test |
| :--- | :--- | :--- | :--- | :--- |
| `RGEHelper` | 39-62 | Numerical algorithm (Interpolation) | `include/SMVacuumDecay/RGE.hpp` | `test_rge.cpp` |
| `V1_loop` | 65-78 | Core physics | `src/EffectivePotential.cpp` | `test_potential.cpp` |
| `get_Mh_calc` | 80-91 | Core physics | `src/RGE.cpp` (Matching) | `test_rge.cpp` |
| `get_nnlo_matching` | 96-108 | Core physics | `src/RGE.cpp` (Matching) | `test_rge.cpp` |
| `get_lambda_eff` | 111-126 | Core physics | `src/EffectivePotential.cpp` | `test_potential.cpp` |
| `betaG1sq` | 128-143 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaG2sq` | 144-161 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaG3sq` | 162-180 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaLambda` | 181-227 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaYt2` | 228-238 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaYb2` | 239-247 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `betaYtau2` | 248-256 | Core physics | `src/RGE.cpp` | `test_rge.cpp` |
| `rk4_single_step` | 258-280 | Numerical algorithm | `include/SMVacuumDecay/numerics/RK4.hpp` | `test_rge.cpp` |
| `rk4_adaptive_step` | 282-304 | Numerical algorithm | `include/SMVacuumDecay/numerics/RK4.hpp` | `test_rge.cpp` |
| `get_pure_sm_lambda` | 312-320 | Core physics | `src/CanonicalBounce.cpp` | `test_bounce.cpp` |
| `integrand_u` | 322-342 | Core physics | `src/CanonicalBounce.cpp` | `test_bounce.cpp` |
| `evaluate_action_at_R` | 344-370 | Numerical algorithm | `src/CanonicalBounce.cpp` | `test_bounce.cpp` |
| `find_minimum_action` | 372-415 | Numerical algorithm | `include/SMVacuumDecay/numerics/GoldenSectionSearch.hpp` | `test_bounce.cpp` |
| `classify_stability` | 418-483 | Core physics / Application | `src/CanonicalBounce.cpp` | `test_scientific.cpp` |
| `main` | 485-528 | CSV/output | `apps/generate_phase_diagram.cpp` | `test_scientific.cpp` |

---

## 2. `solver_analytical.cpp`

| Function / Class | Legacy Lines | Classification | New File Destination | Validation Test |
| :--- | :--- | :--- | :--- | :--- |
| `betaG1sq` to `betaYtau2` | 32-192 | Core physics | *Duplicated in numerical* | `test_rge.cpp` |
| `rk4_single_step` | 195-217 | Numerical algorithm | *Duplicated in numerical* | `test_rge.cpp` |
| `classify_buttazzo` | 219-268 | Core physics / Application | `src/FubiniLipatov.cpp` | `test_FL.cpp` / `test_scientific.cpp` |
| `main` | 270-312 | CSV/output | `apps/generate_phase_diagram.cpp` | `test_scientific.cpp` |

*(Note: While `solver_numerical.cpp` is considered the canonical numerical source, `solver_analytical.cpp` provides the canonical implementation for `classify_buttazzo`, which maps to the Fubini-Lipatov analytical logic).*

## Audit Conclusion
100% of the research logic from the canonical files `solver_numerical.cpp` and `solver_analytical.cpp` has been successfully migrated to the `SMVacuumDecay` architecture.
No code was lost, no magic numbers were accidentally modified, and standard benchmark results reproduce byte-for-byte exact data as confirmed by the continuous validation suites.
