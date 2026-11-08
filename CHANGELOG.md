# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added (formalization branch — Stage A/B review, no change to the scientific implementation)
- `docs/REPOSITORY_AUDIT.md`, `docs/SCIENTIFIC_METHODOLOGY.md`, `docs/NUMERICAL_METHODS_REVIEW.md`:
  independent audit of the original implementation (physics, numerics, data, figures, tests).
- `audit/`: reproducible audit tooling (validated Python port of the legacy solver, loop-order
  power-counting check, independent RG solver, rigorous bounce-action bounds, exact shooting
  bounce, defect-impact and landscape experiments).

### Added (earlier)
- Standard research directory structure (`apps/`, `scripts/`, `data/`, `figures/`).
- Separation of python orchestration from core C++ numerical solvers.
- Integration of `CHANGELOG.md` and `CONTRIBUTING.md`.

### Changed
- Converted analytical vs numerical benchmarks to use the precise `data/numerical_data.csv` backend.
- Renamed and reorganized output figures (e.g., `phase_diagram.png`).
- Fully removed legacy `paper/` directory and build artifacts to maintain focus on the research software itself.
