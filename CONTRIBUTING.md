# Contributing to SMVac

Thank you for your interest in contributing!

## Development Setup
- Clone the repository.
- Ensure you have a C++17 compatible compiler with OpenMP support.
- Python 3.8+ with `pandas`, `numpy`, `matplotlib`, and `scipy`.

## Guidelines
- The core integration solvers are contained in `apps/`.
- Make sure to test changes against both analytical and numerical bounds.
- Update `scripts/` when adding new plotting parameters.

## Pull Requests
- Please provide a clear description of the issue resolved or the feature added.
- All patches should avoid disrupting existing CSV outputs unless there is a physical motivation.
