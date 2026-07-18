# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Standard research directory structure (`apps/`, `scripts/`, `data/`, `figures/`).
- Separation of python orchestration from core C++ numerical solvers.
- Integration of `CHANGELOG.md` and `CONTRIBUTING.md`.

### Changed
- Converted analytical vs numerical benchmarks to use the precise `data/numerical_data.csv` backend.
- Renamed and reorganized output figures (e.g., `phase_diagram.png`).
- Fully removed legacy `paper/` directory and build artifacts to maintain focus on the research software itself.
