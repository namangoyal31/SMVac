# SMVac Migration Plan

## Task 1 — Original Research Files

The following is the exhaustive list of source files, headers, scripts, reports, and data outputs authored or generated natively during the research process (excluding all third-party libraries, toolchains, `.git`, `venv`, and `external/` files):

* `ABSTRACT_PLAN.md`
* `AGENTS.md`
* `AUDIT_REBUTTAL.md`
* `CLAIM_AUDIT.md`
* `CLAIM_IMPACT.md`
* `FIGURE_PLAN.md`
* `FIGURE_REVIEW.md`
* `MISSING_ANALYSIS.md`
* `NEW_OUTLINE.md`
* `ODE_VERIFICATION_PLAN.md`
* `PAPER_REVIEW.md`
* `POTENTIAL_DISCREPANCY_REPORT.md`
* `POTENTIAL_STRUCTURE_REPORT.md`
* `PROJECT_AUDIT.md`
* `SANDWICH_ANALYSIS.md`
* `SANDWICH_REFEREE_TEST.md`
* `SANDWICH_STORY.md`
* `SOLVER_VERIFICATION.md`
* `STORY_ANALYSIS.md`
* `_audit_exp_region.py`
* `_audit_final.py`
* `_audit_refined.py`
* `_audit_sandwich.py`
* `_check_chunk35.py`
* `_check_chunk43846.py`
* `_check_chunk50752.py`
* `_check_chunks.py`
* `_error_check.py`
* `_inspect_diagnostic.py`
* `_sandwich_analysis.py`
* `_sandwich_check.py`
* `_sandwich_detail.py`
* `_sandwich_extent.py`
* `_sandwich_map.py`
* `_sandwich_plot.py`
* `analyze_maxsteps.py`
* `analyze_potential.py`
* `context/context.md`
* `critique_prompt_for_claude.md`
* `debug_shoot.cpp`
* `deriv_solver.cpp`
* `docs/error_analysis_report.md`
* `docs/research_paper_draft.md`
* `docs/stability_analysis_report.md`
* `docs/stability_code_report.md`
* `dummy.cpp`
* `find_interval.cpp`
* `generate_audit.py`
* `metastable_points_mh5.md`
* `paper/figures/analytical_stability_closeup.png`
* `paper/figures/analytical_stability_plot.png`
* `paper/figures/buttazzo_plot.png`
* `paper/figures/c6_phase_diagram.png`
* `paper/figures/c6_sweep_plot.png`
* `paper/figures/error_S_scatter.png`
* `paper/figures/error_boundary_reclassification.png`
* `paper/figures/error_discrepancy_heatmap.png`
* `paper/figures/error_golden_section.png`
* `paper/figures/error_quadrature_convergence.png`
* `paper/figures/error_residual_histogram.png`
* `paper/figures/error_rk4_convergence.png`
* `paper/figures/error_vs_lambda.png`
* `paper/figures/numerical_stability_closeup.png`
* `paper/figures/numerical_stability_contours.png`
* `paper/figures/numerical_stability_plot.png`
* `paper/figures/overlay_closeup_plot.png`
* `paper/figures/overlay_plot.png`
* `paper/main.tex`
* `patch2.py`
* `patch_deriv.py`
* `plot_analytical.py`
* `plot_bisection.py`
* `plot_numerical.py`
* `plot_numerical_contours.py`
* `plot_ode_grid.py`
* `plot_overlay.py`
* `plot_traces.py`
* `potential_analysis_rge.cpp`
* `potential_scan.cpp`
* `prompt_for_claude.md`
* `reformulate_prompt_for_claude.md`
* `results/analytical_data.csv`
* `results/analytical_dense_closeup.csv`
* `results/analytical_stability_closeup.png`
* `results/analytical_stability_plot.png`
* `results/bisection_converged.csv`
* `results/bisection_overshoot.csv`
* `results/bisection_plot.png`
* `results/bisection_undershoot.csv`
* `results/deriv_comparison.csv`
* `results/derivative_comparison.csv`
* `results/derivative_ode_results.csv`
* `results/diagnostic_output.csv`
* `results/diagnostic_overlay.png`
* `results/epsilon_test_results.csv`
* `results/instrumented_trajectories.csv`
* `results/max_steps_diag.csv`
* `results/numerical_data.csv`
* `results/numerical_dense_closeup.csv`
* `results/numerical_stability_closeup.png`
* `results/numerical_stability_contours.png`
* `results/numerical_stability_plot.png`
* `results/ode_convergence_raw.csv`
* `results/ode_dphi_dr.png`
* `results/ode_grid_data.csv`
* `results/ode_grid_data_chunk_0.csv`
* `results/ode_grid_data_chunk_10502.csv`
* `results/ode_grid_data_chunk_15753.csv`
* `results/ode_grid_data_chunk_21004.csv`
* `results/ode_grid_data_chunk_26255.csv`
* `results/ode_grid_data_chunk_31506.csv`
* `results/ode_grid_data_chunk_36757.csv`
* `results/ode_grid_data_chunk_42008.csv`
* `results/ode_grid_data_chunk_47259.csv`
* `results/ode_grid_data_chunk_5251.csv`
* `results/ode_grid_data_chunk_52510.csv`
* `results/ode_grid_data_chunk_57761.csv`
* `results/ode_phi_r.png`
* `results/ode_profile.csv`
* `results/ode_profile_105.csv`
* `results/ode_profile_106.csv`
* `results/ode_profile_110.csv`
* `results/ode_stability_plot.png`
* `results/ode_stability_plot_partial.png`
* `results/ode_validation_raw.csv`
* `results/overlay_closeup_plot.png`
* `results/overlay_plot.png`
* `results/partial_ode_grid.csv`
* `results/partial_ode_plot.png`
* `results/potential_discrepancy.png`
* `results/potential_scan.csv`
* `results/potential_structure_plot.png`
* `results/potential_verification_grid.csv`
* `results/profile_ana.csv`
* `results/profile_fd.csv`
* `results/sandwich_neg_fraction.png`
* `results/sandwich_panel.png`
* `results/sandwich_zoom_Mh100.png`
* `results/stability_closeup_plot.png`
* `results/stability_plot.png`
* `results/stress_test_results.csv`
* `results/sweep_summary.csv`
* `results/trace_0.csv`
* `results/trace_1.csv`
* `results/trace_2.csv`
* `results/trace_plot.png`
* `results/verify_stationary_points.csv`
* `run_analytical.py`
* `run_bisection.cpp`
* `run_dense.py`
* `run_numerical.py`
* `run_ode_grid.py`
* `solver_analytical.cpp`
* `solver_analytical_dense.cpp`
* `solver_canonical.cpp`
* `solver_canonical_A.cpp`
* `solver_canonical_B.cpp`
* `solver_diagnostic.cpp`
* `solver_instrumented.cpp`
* `solver_maxsteps_diag.cpp`
* `solver_numerical.cpp`
* `solver_numerical_dense.cpp`
* `solver_ode_diagnostic.cpp`
* `solver_ode_grid.cpp`
* `src/cosmology/higgs_inflation.py`
* `src/cosmology/landscape_bertini_solver.py`
* `src/cosmology/landscape_gradient_flow.py`
* `src/plotting/plot_error_scaling.py`
* `src/plotting/plot_interim.py`
* `src/plotting/plot_phase_diagram.py`
* `src/plotting/plot_precision_comparison.py`
* `src/plotting/plot_scatter.py`
* `src/sm_stability/bounce_action_solver_unopt.cpp`
* `src/sm_stability/error_scaling_test.cpp`
* `src/sm_stability/fix.py`
* `src/sm_stability/fix_grid.py`
* `src/sm_stability/patch_fix.py`
* `src/sm_stability/patch_main.py`
* `src/sm_stability/patch_main2.py`
* `src/sm_stability/patch_opt.py`
* `src/sm_stability/patch_physics.py`
* `src/sm_stability/sm_rge_solver.cpp`
* `src/sm_stability/test.py`
* `src/sm_stability/test_exact.cpp`
* `src/sm_stability/test_interp.cpp`
* `src/sm_stability/test_opt.py`
* `src/sm_stability/test_rge.py`
* `src/sm_stability/test_solver.cpp`
* `src/sm_stability/test_solver2.cpp`
* `src/sm_stability/vacuum_precision_test.cpp`
* `src/sm_stability/vacuum_precision_test_quad.cpp`
* `src/sm_stability\error_scaling_test2.cpp`
* `stationary_finder.cpp`
* `stress_solver.cpp`
* `sweep_interval.cpp`
* `test_dv.cpp`
* `test_nan.cpp`
* `test_v.cpp`
* `validation_solver.cpp`
* `verify_continuity.py`
* `verify_potential.cpp`
* `write_test_nan.py`

## Task 2 — Grouped by Scientific Purpose

### Markdown Reports
* `ABSTRACT_PLAN.md`
* `AGENTS.md`
* `AUDIT_REBUTTAL.md`
* `CLAIM_AUDIT.md`
* `CLAIM_IMPACT.md`
* `FIGURE_PLAN.md`
* `FIGURE_REVIEW.md`
* `MISSING_ANALYSIS.md`
* `NEW_OUTLINE.md`
* `ODE_VERIFICATION_PLAN.md`
* `PAPER_REVIEW.md`
* `POTENTIAL_DISCREPANCY_REPORT.md`
* `POTENTIAL_STRUCTURE_REPORT.md`
* `PROJECT_AUDIT.md`
* `SANDWICH_ANALYSIS.md`
* `SANDWICH_REFEREE_TEST.md`
* `SANDWICH_STORY.md`
* `SOLVER_VERIFICATION.md`
* `STORY_ANALYSIS.md`
* `context/context.md`
* `critique_prompt_for_claude.md`
* `docs/error_analysis_report.md`
* `docs/research_paper_draft.md`
* `docs/stability_analysis_report.md`
* `docs/stability_code_report.md`
* `metastable_points_mh5.md`
* `prompt_for_claude.md`
* `reformulate_prompt_for_claude.md`

### Other scripts
* `_audit_exp_region.py`
* `_audit_final.py`
* `_audit_refined.py`
* `_audit_sandwich.py`
* `_check_chunk35.py`
* `_check_chunk43846.py`
* `_check_chunk50752.py`
* `_check_chunks.py`
* `_error_check.py`
* `dummy.cpp`
* `generate_audit.py`
* `patch2.py`
* `patch_deriv.py`
* `potential_scan.cpp`
* `solver_canonical_A.cpp`
* `solver_canonical_B.cpp`
* `solver_ode_grid.cpp`
* `test_nan.cpp`
* `verify_continuity.py`
* `write_test_nan.py`

### Bounce diagnostics / Numerical stress tests
* `_inspect_diagnostic.py`
* `analyze_maxsteps.py`
* `debug_shoot.cpp`
* `deriv_solver.cpp`
* `results/diagnostic_output.csv`
* `results/diagnostic_overlay.png`
* `results/instrumented_trajectories.csv`
* `results/stress_test_results.csv`
* `solver_diagnostic.cpp`
* `solver_instrumented.cpp`
* `solver_maxsteps_diag.cpp`
* `solver_ode_diagnostic.cpp`
* `stress_solver.cpp`
* `validation_solver.cpp`

### Sandwich region diagnostics
* `_sandwich_analysis.py`
* `_sandwich_check.py`
* `_sandwich_detail.py`
* `_sandwich_extent.py`
* `_sandwich_map.py`
* `_sandwich_plot.py`
* `find_interval.cpp`

### Potential verification
* `analyze_potential.py`
* `potential_analysis_rge.cpp`
* `test_dv.cpp`
* `test_v.cpp`
* `verify_potential.cpp`

### Paper
* `paper/figures/analytical_stability_closeup.png`
* `paper/figures/analytical_stability_plot.png`
* `paper/figures/buttazzo_plot.png`
* `paper/figures/c6_phase_diagram.png`
* `paper/figures/c6_sweep_plot.png`
* `paper/figures/error_S_scatter.png`
* `paper/figures/error_boundary_reclassification.png`
* `paper/figures/error_discrepancy_heatmap.png`
* `paper/figures/error_golden_section.png`
* `paper/figures/error_quadrature_convergence.png`
* `paper/figures/error_residual_histogram.png`
* `paper/figures/error_rk4_convergence.png`
* `paper/figures/error_vs_lambda.png`
* `paper/figures/numerical_stability_closeup.png`
* `paper/figures/numerical_stability_contours.png`
* `paper/figures/numerical_stability_plot.png`
* `paper/figures/overlay_closeup_plot.png`
* `paper/figures/overlay_plot.png`
* `paper/main.tex`

### Plotting
* `plot_analytical.py`
* `plot_bisection.py`
* `plot_numerical.py`
* `plot_numerical_contours.py`
* `plot_ode_grid.py`
* `plot_overlay.py`
* `plot_traces.py`
* `results/analytical_stability_plot.png`
* `results/bisection_plot.png`
* `results/numerical_stability_plot.png`
* `results/ode_stability_plot.png`
* `results/ode_stability_plot_partial.png`
* `results/overlay_closeup_plot.png`
* `results/overlay_plot.png`
* `results/partial_ode_plot.png`
* `results/potential_structure_plot.png`
* `results/stability_closeup_plot.png`
* `results/stability_plot.png`
* `results/trace_plot.png`
* `src/plotting/plot_error_scaling.py`
* `src/plotting/plot_interim.py`
* `src/plotting/plot_phase_diagram.py`
* `src/plotting/plot_precision_comparison.py`
* `src/plotting/plot_scatter.py`

### Data & Figures (Results)
* `results/analytical_data.csv`
* `results/analytical_dense_closeup.csv`
* `results/analytical_stability_closeup.png`
* `results/bisection_converged.csv`
* `results/bisection_overshoot.csv`
* `results/bisection_undershoot.csv`
* `results/deriv_comparison.csv`
* `results/derivative_comparison.csv`
* `results/derivative_ode_results.csv`
* `results/epsilon_test_results.csv`
* `results/max_steps_diag.csv`
* `results/numerical_data.csv`
* `results/numerical_dense_closeup.csv`
* `results/numerical_stability_closeup.png`
* `results/numerical_stability_contours.png`
* `results/ode_convergence_raw.csv`
* `results/ode_dphi_dr.png`
* `results/ode_grid_data.csv`
* `results/ode_grid_data_chunk_0.csv`
* `results/ode_grid_data_chunk_10502.csv`
* `results/ode_grid_data_chunk_15753.csv`
* `results/ode_grid_data_chunk_21004.csv`
* `results/ode_grid_data_chunk_26255.csv`
* `results/ode_grid_data_chunk_31506.csv`
* `results/ode_grid_data_chunk_36757.csv`
* `results/ode_grid_data_chunk_42008.csv`
* `results/ode_grid_data_chunk_47259.csv`
* `results/ode_grid_data_chunk_5251.csv`
* `results/ode_grid_data_chunk_52510.csv`
* `results/ode_grid_data_chunk_57761.csv`
* `results/ode_phi_r.png`
* `results/ode_profile.csv`
* `results/ode_profile_105.csv`
* `results/ode_profile_106.csv`
* `results/ode_profile_110.csv`
* `results/ode_validation_raw.csv`
* `results/partial_ode_grid.csv`
* `results/potential_discrepancy.png`
* `results/potential_scan.csv`
* `results/potential_verification_grid.csv`
* `results/profile_ana.csv`
* `results/profile_fd.csv`
* `results/sandwich_neg_fraction.png`
* `results/sandwich_panel.png`
* `results/sandwich_zoom_Mh100.png`
* `results/sweep_summary.csv`
* `results/trace_0.csv`
* `results/trace_1.csv`
* `results/trace_2.csv`
* `results/verify_stationary_points.csv`

### Phase-boundary scan wrappers
* `run_analytical.py`
* `run_dense.py`
* `run_numerical.py`
* `run_ode_grid.py`

### Canonical ODE solver
* `run_bisection.cpp`
* `solver_canonical.cpp`
* `sweep_interval.cpp`

### Fubini-Lipatov approximation
* `solver_analytical.cpp`
* `solver_analytical_dense.cpp`

### Three-loop RGE evolution / Exact Integration
* `solver_numerical.cpp`
* `solver_numerical_dense.cpp`

### Obsolete prototyping code
* `src/cosmology/higgs_inflation.py`
* `src/cosmology/landscape_bertini_solver.py`
* `src/cosmology/landscape_gradient_flow.py`
* `src/sm_stability/bounce_action_solver_unopt.cpp`
* `src/sm_stability/error_scaling_test.cpp`
* `src/sm_stability/fix.py`
* `src/sm_stability/fix_grid.py`
* `src/sm_stability/patch_fix.py`
* `src/sm_stability/patch_main.py`
* `src/sm_stability/patch_main2.py`
* `src/sm_stability/patch_opt.py`
* `src/sm_stability/patch_physics.py`
* `src/sm_stability/sm_rge_solver.cpp`
* `src/sm_stability/test.py`
* `src/sm_stability/test_exact.cpp`
* `src/sm_stability/test_interp.cpp`
* `src/sm_stability/test_opt.py`
* `src/sm_stability/test_rge.py`
* `src/sm_stability/test_solver.cpp`
* `src/sm_stability/test_solver2.cpp`
* `src/sm_stability/vacuum_precision_test.cpp`
* `src/sm_stability/vacuum_precision_test_quad.cpp`
* `src/sm_stability\error_scaling_test2.cpp`

### Stationary-point analysis
* `stationary_finder.cpp`

## Task 3 — Final vs. Intermediate Versions

| File | Status | Reason |
|------|--------|--------|
| `_audit_exp_region.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_audit_final.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_audit_refined.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_audit_sandwich.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_check_chunk35.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_check_chunk43846.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_check_chunk50752.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_check_chunks.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_error_check.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_inspect_diagnostic.py` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `_sandwich_analysis.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_sandwich_check.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_sandwich_detail.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_sandwich_extent.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_sandwich_map.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `_sandwich_plot.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `analyze_maxsteps.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `analyze_potential.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `debug_shoot.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `deriv_solver.cpp` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `dummy.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `find_interval.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `generate_audit.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `patch2.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `patch_deriv.py` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `plot_analytical.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_bisection.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_numerical.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_numerical_contours.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_ode_grid.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_overlay.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `plot_traces.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `potential_analysis_rge.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `potential_scan.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `run_analytical.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `run_bisection.cpp` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `run_dense.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `run_numerical.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `run_ode_grid.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `solver_analytical.cpp` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `solver_analytical_dense.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `solver_canonical.cpp` | **EXPERIMENTAL** | Proofs of concept like canonical solver without integration scaling. |
| `solver_canonical_A.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `solver_canonical_B.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `solver_diagnostic.cpp` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `solver_instrumented.cpp` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `solver_maxsteps_diag.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `solver_numerical.cpp` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `solver_numerical_dense.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `solver_ode_diagnostic.cpp` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `solver_ode_grid.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `src/cosmology/higgs_inflation.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `src/cosmology/landscape_bertini_solver.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `src/cosmology/landscape_gradient_flow.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `src/plotting/plot_error_scaling.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `src/plotting/plot_interim.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `src/plotting/plot_phase_diagram.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `src/plotting/plot_precision_comparison.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `src/plotting/plot_scatter.py` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `src/sm_stability/bounce_action_solver_unopt.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/error_scaling_test.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/fix.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/fix_grid.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/patch_fix.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/patch_main.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/patch_main2.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/patch_opt.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/patch_physics.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/sm_rge_solver.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_exact.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_interp.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_opt.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_rge.py` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_solver.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/test_solver2.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/vacuum_precision_test.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability/vacuum_precision_test_quad.cpp` | **OBSOLETE** | Legacy prototyping code superseded by root solvers. |
| `src/sm_stability\error_scaling_test2.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `stationary_finder.cpp` | **FINAL** | Production pipeline component used for final phase diagram or physics result. |
| `stress_solver.cpp` | **DEBUGGING** | Instrumented copy used strictly to extract internal state during explosions. |
| `sweep_interval.cpp` | **VALIDATION** | Used to explicitly verify theoretical formulations or equations of motion. |
| `test_dv.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `test_nan.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `test_v.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `validation_solver.cpp` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `verify_continuity.py` | **INTERMEDIATE** | One-off script or intermediate patch. |
| `verify_potential.cpp` | **VALIDATION** | Used to explicitly verify theoretical formulations or equations of motion. |
| `write_test_nan.py` | **INTERMEDIATE** | One-off script or intermediate patch. |

## Task 4 — Scientific Results Inventory

| Title | Code Used | Figures | Numerical Outputs | Pub Quality | GitHub Quality |
|-------|-----------|---------|-------------------|-------------|----------------|
| **Three-loop RGE evolution** | `solver_numerical.cpp` | `error_vs_lambda.png` | `numerical_data.csv` | Yes | Yes |
| **Vacuum stability phase diagram** | `run_numerical.py` | `numerical_stability_plot.png`, `numerical_stability_contours.png` | `numerical_data.csv` | Yes | Yes |
| **FL approximation** | `run_analytical.py` | `analytical_stability_plot.png` | `analytical_data.csv` | Yes | Yes |
| **High-field approximation validation** | `plot_overlay.py` | `overlay_plot.png`, `overlay_closeup_plot.png` | - | Yes | Yes |
| **Sandwich analysis** | `find_interval.cpp` | `sandwich_panel.png` | - | Yes | Yes |
| **Stationary points** | `stationary_finder.cpp` | `potential_structure_plot.png` | `verify_stationary_points.csv` | Yes | Yes |
| **Potential reconstruction** | `verify_potential.cpp` | `potential_discrepancy.png` | `potential_verification_grid.csv` | Yes | Yes |
| **ODE convergence investigation** | `sweep_interval.cpp` | `trace_plot.png` | `sweep_summary.csv`, `trace_*.csv` | Yes (Appendix) | Yes |
| **Canonical bounce** | `run_bisection.cpp` | `bisection_plot.png` | `bisection_converged.csv` | Yes | Yes |
| **Benchmark Mh=5 Mt=105** | `run_bisection.cpp` | `bisection_plot.png` | `bisection_converged.csv` | Yes | Yes |
| **Stress tests** | `stress_solver.cpp` | `error_residual_histogram.png` | `stress_test_results.csv` | No | No (Internal) |

## Task 5 — SMVac Repository Contents & Migration Map

The `SMVac` repository will be structured cleanly to accompany the publication. No compiled binaries, external libraries, or temporary files are included.

### Proposed `SMVac/` Structure
```
SMVac/
├── README.md
├── LICENSE
├── docs/                   # Finalized markdown reports
├── src/
│   ├── numerical/          # Exact numerical ODE solver
│   ├── analytical/         # Conformal approximation solver
│   └── diagnostics/        # Stationary finder, bisection, interval sweeps
├── python/
│   ├── runners/            # Parallel grid execution
│   └── plotting/           # Figure generation
├── data/                   # Output CSVs (LFS)
├── figures/                # Output PNGs
└── paper/                  # LaTeX manuscript
```

### Action Plan

| File | Action | New Location (if kept) |
|------|--------|------------------------|
| `ABSTRACT_PLAN.md` | **IGNORE** | - |
| `AGENTS.md` | **IGNORE** | - |
| `AUDIT_REBUTTAL.md` | **IGNORE** | - |
| `CLAIM_AUDIT.md` | **IGNORE** | - |
| `CLAIM_IMPACT.md` | **IGNORE** | - |
| `FIGURE_PLAN.md` | **IGNORE** | - |
| `FIGURE_REVIEW.md` | **IGNORE** | - |
| `MISSING_ANALYSIS.md` | **IGNORE** | - |
| `NEW_OUTLINE.md` | **IGNORE** | - |
| `ODE_VERIFICATION_PLAN.md` | **IGNORE** | - |
| `PAPER_REVIEW.md` | **IGNORE** | - |
| `POTENTIAL_DISCREPANCY_REPORT.md` | **IGNORE** | - |
| `POTENTIAL_STRUCTURE_REPORT.md` | **IGNORE** | - |
| `PROJECT_AUDIT.md` | **IGNORE** | - |
| `SANDWICH_ANALYSIS.md` | **IGNORE** | - |
| `SANDWICH_REFEREE_TEST.md` | **IGNORE** | - |
| `SANDWICH_STORY.md` | **IGNORE** | - |
| `SOLVER_VERIFICATION.md` | **IGNORE** | - |
| `STORY_ANALYSIS.md` | **IGNORE** | - |
| `_audit_exp_region.py` | **ARCHIVE** | - |
| `_audit_final.py` | **ARCHIVE** | - |
| `_audit_refined.py` | **ARCHIVE** | - |
| `_audit_sandwich.py` | **ARCHIVE** | - |
| `_check_chunk35.py` | **ARCHIVE** | - |
| `_check_chunk43846.py` | **ARCHIVE** | - |
| `_check_chunk50752.py` | **ARCHIVE** | - |
| `_check_chunks.py` | **ARCHIVE** | - |
| `_error_check.py` | **ARCHIVE** | - |
| `_inspect_diagnostic.py` | **ARCHIVE** | - |
| `_sandwich_analysis.py` | **ARCHIVE** | - |
| `_sandwich_check.py` | **ARCHIVE** | - |
| `_sandwich_detail.py` | **ARCHIVE** | - |
| `_sandwich_extent.py` | **ARCHIVE** | - |
| `_sandwich_map.py` | **ARCHIVE** | - |
| `_sandwich_plot.py` | **ARCHIVE** | - |
| `analyze_maxsteps.py` | **IGNORE** | - |
| `analyze_potential.py` | **IGNORE** | - |
| `context/context.md` | **IGNORE** | - |
| `critique_prompt_for_claude.md` | **IGNORE** | - |
| `debug_shoot.cpp` | **IGNORE** | - |
| `deriv_solver.cpp` | **ARCHIVE** | - |
| `docs/error_analysis_report.md` | **IGNORE** | - |
| `docs/research_paper_draft.md` | **IGNORE** | - |
| `docs/stability_analysis_report.md` | **IGNORE** | - |
| `docs/stability_code_report.md` | **IGNORE** | - |
| `dummy.cpp` | **IGNORE** | - |
| `find_interval.cpp` | **IGNORE** | - |
| `generate_audit.py` | **IGNORE** | - |
| `metastable_points_mh5.md` | **IGNORE** | - |
| `paper/figures/analytical_stability_closeup.png` | **KEEP** | `paper/figures/analytical_stability_closeup.png` |
| `paper/figures/analytical_stability_plot.png` | **KEEP** | `paper/figures/analytical_stability_plot.png` |
| `paper/figures/buttazzo_plot.png` | **KEEP** | `paper/figures/buttazzo_plot.png` |
| `paper/figures/c6_phase_diagram.png` | **KEEP** | `paper/figures/c6_phase_diagram.png` |
| `paper/figures/c6_sweep_plot.png` | **KEEP** | `paper/figures/c6_sweep_plot.png` |
| `paper/figures/error_S_scatter.png` | **KEEP** | `paper/figures/error_S_scatter.png` |
| `paper/figures/error_boundary_reclassification.png` | **KEEP** | `paper/figures/error_boundary_reclassification.png` |
| `paper/figures/error_discrepancy_heatmap.png` | **KEEP** | `paper/figures/error_discrepancy_heatmap.png` |
| `paper/figures/error_golden_section.png` | **KEEP** | `paper/figures/error_golden_section.png` |
| `paper/figures/error_quadrature_convergence.png` | **KEEP** | `paper/figures/error_quadrature_convergence.png` |
| `paper/figures/error_residual_histogram.png` | **KEEP** | `paper/figures/error_residual_histogram.png` |
| `paper/figures/error_rk4_convergence.png` | **KEEP** | `paper/figures/error_rk4_convergence.png` |
| `paper/figures/error_vs_lambda.png` | **KEEP** | `paper/figures/error_vs_lambda.png` |
| `paper/figures/numerical_stability_closeup.png` | **KEEP** | `paper/figures/numerical_stability_closeup.png` |
| `paper/figures/numerical_stability_contours.png` | **KEEP** | `paper/figures/numerical_stability_contours.png` |
| `paper/figures/numerical_stability_plot.png` | **KEEP** | `paper/figures/numerical_stability_plot.png` |
| `paper/figures/overlay_closeup_plot.png` | **KEEP** | `paper/figures/overlay_closeup_plot.png` |
| `paper/figures/overlay_plot.png` | **KEEP** | `paper/figures/overlay_plot.png` |
| `paper/main.tex` | **KEEP** | `paper/main.tex` |
| `patch2.py` | **IGNORE** | - |
| `patch_deriv.py` | **IGNORE** | - |
| `plot_analytical.py` | **KEEP** | `python/plotting/plot_analytical.py` |
| `plot_bisection.py` | **IGNORE** | - |
| `plot_numerical.py` | **KEEP** | `python/plotting/plot_numerical.py` |
| `plot_numerical_contours.py` | **KEEP** | `python/plotting/plot_numerical_contours.py` |
| `plot_ode_grid.py` | **IGNORE** | - |
| `plot_overlay.py` | **KEEP** | `python/plotting/plot_overlay.py` |
| `plot_traces.py` | **IGNORE** | - |
| `potential_analysis_rge.cpp` | **IGNORE** | - |
| `potential_scan.cpp` | **IGNORE** | - |
| `prompt_for_claude.md` | **IGNORE** | - |
| `reformulate_prompt_for_claude.md` | **IGNORE** | - |
| `results/analytical_data.csv` | **KEEP** | `data/analytical_data.csv` |
| `results/analytical_dense_closeup.csv` | **KEEP** | `data/analytical_dense_closeup.csv` |
| `results/analytical_stability_closeup.png` | **IGNORE** | - |
| `results/analytical_stability_plot.png` | **IGNORE** | - |
| `results/bisection_converged.csv` | **KEEP** | `data/bisection_converged.csv` |
| `results/bisection_overshoot.csv` | **KEEP** | `data/bisection_overshoot.csv` |
| `results/bisection_plot.png` | **IGNORE** | - |
| `results/bisection_undershoot.csv` | **KEEP** | `data/bisection_undershoot.csv` |
| `results/deriv_comparison.csv` | **KEEP** | `data/deriv_comparison.csv` |
| `results/derivative_comparison.csv` | **KEEP** | `data/derivative_comparison.csv` |
| `results/derivative_ode_results.csv` | **KEEP** | `data/derivative_ode_results.csv` |
| `results/diagnostic_output.csv` | **KEEP** | `data/diagnostic_output.csv` |
| `results/diagnostic_overlay.png` | **IGNORE** | - |
| `results/epsilon_test_results.csv` | **KEEP** | `data/epsilon_test_results.csv` |
| `results/instrumented_trajectories.csv` | **KEEP** | `data/instrumented_trajectories.csv` |
| `results/max_steps_diag.csv` | **KEEP** | `data/max_steps_diag.csv` |
| `results/numerical_data.csv` | **KEEP** | `data/numerical_data.csv` |
| `results/numerical_dense_closeup.csv` | **KEEP** | `data/numerical_dense_closeup.csv` |
| `results/numerical_stability_closeup.png` | **IGNORE** | - |
| `results/numerical_stability_contours.png` | **IGNORE** | - |
| `results/numerical_stability_plot.png` | **IGNORE** | - |
| `results/ode_convergence_raw.csv` | **KEEP** | `data/ode_convergence_raw.csv` |
| `results/ode_dphi_dr.png` | **IGNORE** | - |
| `results/ode_grid_data.csv` | **KEEP** | `data/ode_grid_data.csv` |
| `results/ode_grid_data_chunk_0.csv` | **KEEP** | `data/ode_grid_data_chunk_0.csv` |
| `results/ode_grid_data_chunk_10502.csv` | **KEEP** | `data/ode_grid_data_chunk_10502.csv` |
| `results/ode_grid_data_chunk_15753.csv` | **KEEP** | `data/ode_grid_data_chunk_15753.csv` |
| `results/ode_grid_data_chunk_21004.csv` | **KEEP** | `data/ode_grid_data_chunk_21004.csv` |
| `results/ode_grid_data_chunk_26255.csv` | **KEEP** | `data/ode_grid_data_chunk_26255.csv` |
| `results/ode_grid_data_chunk_31506.csv` | **KEEP** | `data/ode_grid_data_chunk_31506.csv` |
| `results/ode_grid_data_chunk_36757.csv` | **KEEP** | `data/ode_grid_data_chunk_36757.csv` |
| `results/ode_grid_data_chunk_42008.csv` | **KEEP** | `data/ode_grid_data_chunk_42008.csv` |
| `results/ode_grid_data_chunk_47259.csv` | **KEEP** | `data/ode_grid_data_chunk_47259.csv` |
| `results/ode_grid_data_chunk_5251.csv` | **KEEP** | `data/ode_grid_data_chunk_5251.csv` |
| `results/ode_grid_data_chunk_52510.csv` | **KEEP** | `data/ode_grid_data_chunk_52510.csv` |
| `results/ode_grid_data_chunk_57761.csv` | **KEEP** | `data/ode_grid_data_chunk_57761.csv` |
| `results/ode_phi_r.png` | **IGNORE** | - |
| `results/ode_profile.csv` | **KEEP** | `data/ode_profile.csv` |
| `results/ode_profile_105.csv` | **KEEP** | `data/ode_profile_105.csv` |
| `results/ode_profile_106.csv` | **KEEP** | `data/ode_profile_106.csv` |
| `results/ode_profile_110.csv` | **KEEP** | `data/ode_profile_110.csv` |
| `results/ode_stability_plot.png` | **IGNORE** | - |
| `results/ode_stability_plot_partial.png` | **IGNORE** | - |
| `results/ode_validation_raw.csv` | **KEEP** | `data/ode_validation_raw.csv` |
| `results/overlay_closeup_plot.png` | **IGNORE** | - |
| `results/overlay_plot.png` | **IGNORE** | - |
| `results/partial_ode_grid.csv` | **KEEP** | `data/partial_ode_grid.csv` |
| `results/partial_ode_plot.png` | **IGNORE** | - |
| `results/potential_discrepancy.png` | **IGNORE** | - |
| `results/potential_scan.csv` | **KEEP** | `data/potential_scan.csv` |
| `results/potential_structure_plot.png` | **IGNORE** | - |
| `results/potential_verification_grid.csv` | **KEEP** | `data/potential_verification_grid.csv` |
| `results/profile_ana.csv` | **KEEP** | `data/profile_ana.csv` |
| `results/profile_fd.csv` | **KEEP** | `data/profile_fd.csv` |
| `results/sandwich_neg_fraction.png` | **IGNORE** | - |
| `results/sandwich_panel.png` | **IGNORE** | - |
| `results/sandwich_zoom_Mh100.png` | **IGNORE** | - |
| `results/stability_closeup_plot.png` | **IGNORE** | - |
| `results/stability_plot.png` | **IGNORE** | - |
| `results/stress_test_results.csv` | **KEEP** | `data/stress_test_results.csv` |
| `results/sweep_summary.csv` | **KEEP** | `data/sweep_summary.csv` |
| `results/trace_0.csv` | **KEEP** | `data/trace_0.csv` |
| `results/trace_1.csv` | **KEEP** | `data/trace_1.csv` |
| `results/trace_2.csv` | **KEEP** | `data/trace_2.csv` |
| `results/trace_plot.png` | **IGNORE** | - |
| `results/verify_stationary_points.csv` | **KEEP** | `data/verify_stationary_points.csv` |
| `run_analytical.py` | **KEEP** | `python/runners/run_analytical.py` |
| `run_bisection.cpp` | **KEEP** | `src/diagnostics/run_bisection.cpp` |
| `run_dense.py` | **IGNORE** | - |
| `run_numerical.py` | **KEEP** | `python/runners/run_numerical.py` |
| `run_ode_grid.py` | **IGNORE** | - |
| `solver_analytical.cpp` | **KEEP** | `src/analytical/main.cpp` |
| `solver_analytical_dense.cpp` | **IGNORE** | - |
| `solver_canonical.cpp` | **KEEP** | `src/diagnostics/solver_canonical.cpp` |
| `solver_canonical_A.cpp` | **IGNORE** | - |
| `solver_canonical_B.cpp` | **IGNORE** | - |
| `solver_diagnostic.cpp` | **ARCHIVE** | - |
| `solver_instrumented.cpp` | **ARCHIVE** | - |
| `solver_maxsteps_diag.cpp` | **IGNORE** | - |
| `solver_numerical.cpp` | **KEEP** | `src/numerical/main.cpp` |
| `solver_numerical_dense.cpp` | **IGNORE** | - |
| `solver_ode_diagnostic.cpp` | **ARCHIVE** | - |
| `solver_ode_grid.cpp` | **IGNORE** | - |
| `src/cosmology/higgs_inflation.py` | **ARCHIVE** | - |
| `src/cosmology/landscape_bertini_solver.py` | **ARCHIVE** | - |
| `src/cosmology/landscape_gradient_flow.py` | **ARCHIVE** | - |
| `src/plotting/plot_error_scaling.py` | **IGNORE** | - |
| `src/plotting/plot_interim.py` | **IGNORE** | - |
| `src/plotting/plot_phase_diagram.py` | **IGNORE** | - |
| `src/plotting/plot_precision_comparison.py` | **IGNORE** | - |
| `src/plotting/plot_scatter.py` | **IGNORE** | - |
| `src/sm_stability/bounce_action_solver_unopt.cpp` | **ARCHIVE** | - |
| `src/sm_stability/error_scaling_test.cpp` | **ARCHIVE** | - |
| `src/sm_stability/fix.py` | **ARCHIVE** | - |
| `src/sm_stability/fix_grid.py` | **ARCHIVE** | - |
| `src/sm_stability/patch_fix.py` | **ARCHIVE** | - |
| `src/sm_stability/patch_main.py` | **ARCHIVE** | - |
| `src/sm_stability/patch_main2.py` | **ARCHIVE** | - |
| `src/sm_stability/patch_opt.py` | **ARCHIVE** | - |
| `src/sm_stability/patch_physics.py` | **ARCHIVE** | - |
| `src/sm_stability/sm_rge_solver.cpp` | **ARCHIVE** | - |
| `src/sm_stability/test.py` | **ARCHIVE** | - |
| `src/sm_stability/test_exact.cpp` | **ARCHIVE** | - |
| `src/sm_stability/test_interp.cpp` | **ARCHIVE** | - |
| `src/sm_stability/test_opt.py` | **ARCHIVE** | - |
| `src/sm_stability/test_rge.py` | **ARCHIVE** | - |
| `src/sm_stability/test_solver.cpp` | **ARCHIVE** | - |
| `src/sm_stability/test_solver2.cpp` | **ARCHIVE** | - |
| `src/sm_stability/vacuum_precision_test.cpp` | **ARCHIVE** | - |
| `src/sm_stability/vacuum_precision_test_quad.cpp` | **ARCHIVE** | - |
| `src/sm_stability\error_scaling_test2.cpp` | **IGNORE** | - |
| `stationary_finder.cpp` | **KEEP** | `src/diagnostics/stationary_finder.cpp` |
| `stress_solver.cpp` | **ARCHIVE** | - |
| `sweep_interval.cpp` | **KEEP** | `src/diagnostics/sweep_interval.cpp` |
| `test_dv.cpp` | **IGNORE** | - |
| `test_nan.cpp` | **IGNORE** | - |
| `test_v.cpp` | **IGNORE** | - |
| `validation_solver.cpp` | **ARCHIVE** | - |
| `verify_continuity.py` | **IGNORE** | - |
| `verify_potential.cpp` | **KEEP** | `src/diagnostics/verify_potential.cpp` |
| `write_test_nan.py` | **IGNORE** | - |
