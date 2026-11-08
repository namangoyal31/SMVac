# Archived documents (historical; superseded)

These documents were written during the predecessor project ("Threshold") and the first version of this
repository. They are kept unchanged for the record, and **they are not part of the current scientific argument**.
Many of their statements were later found to be incorrect or unsupported:

* "exact" bounce actions;
* error budgets of $10^{-11}$;
* a 20–30% failure of the conformal approximation near the metastability boundary;
* the "sandwiched unstable region";
* the $\phi^6$ "rescue".

See `docs/REPOSITORY_AUDIT.md` and `docs/SCIENTIFIC_METHODOLOGY.md` for the audit, and `docs/VALIDATION_REPORT.md`
for the current results.

| File | Original purpose |
|---|---|
| `PROJECT_AUDIT.md`, `SMVAC_MIGRATION_PLAN.md`, `MINIMAL_MIGRATION_PLAN.md` | inventory of the predecessor repository and migration plans |
| `research_paper_draft.md`, `ABSTRACT_PLAN.md`, `NEW_OUTLINE.md`, `STORY_ANALYSIS.md`, `PAPER_REVIEW.md`, `FIGURE_PLAN.md`, `FIGURE_REVIEW.md`, `MISSING_ANALYSIS.md` | manuscript planning |
| `SOLVER_VERIFICATION.md`, `CLAIM_AUDIT.md`, `CLAIM_IMPACT.md`, `AUDIT_REBUTTAL.md` | mock referee review of the legacy solver and a rebuttal. They correctly noted that no bounce equation was solved. |
| `SANDWICH_ANALYSIS.md`, `SANDWICH_REFEREE_TEST.md`, `SANDWICH_STORY.md`, `metastable_points_mh5.md` | analysis of the negative-action region, an artefact of the legacy construction (METHODOLOGY §6.3) |
| `ODE_VERIFICATION_PLAN.md`, `POTENTIAL_STRUCTURE_REPORT.md`, `POTENTIAL_DISCREPANCY_REPORT.md` | plans and notes for an ODE-based check at (5, 105) GeV |
| `error_analysis_report.md`, `stability_analysis_report.md`, `stability_code_report.md` | numerical-error claims and code reports for the legacy solver |

Removed from the repository, with their history preserved in git:

* three AI-assistant prompt files and an AI-agent configuration file (`AGENTS.md`);
* `docs/README`, which was the README of the unrelated Bertini package. No file in the repository uses Bertini.
