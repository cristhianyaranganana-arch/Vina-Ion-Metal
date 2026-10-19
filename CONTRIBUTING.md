# Contributing to VinaMetal

Thank you for your interest in improving VinaMetal. This document describes how
to build the project, run the tests, and submit changes.

## Code of Conduct

Be respectful and constructive. Technical disagreement is welcome; personal
attacks, harassment and discriminatory language are not.

## Reporting Bugs

Open an issue and include:

- VinaMetal version (`vina --version`) and commit hash.
- Operating system, compiler version, CMake version and Boost version.
- The exact command line or configuration file used.
- The complete log, including the `VinaMetal:` startup lines.
- Whether the problem also occurs with `--metals off`. If it does, it is likely
  inherited from AutoDock Vina 1.2.7.

Do not attach proprietary receptor structures. A minimal synthetic reproduction
is always preferred.

## Building

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

Requirements: CMake ≥ 3.12, a C++11 compiler and Boost ≥ 1.65 with the
`system`, `thread`, `serialization`, `filesystem` and `program_options`
libraries. See [docs/INSTALL.md](docs/INSTALL.md) for platform details.

## Tests

All contributions must keep the following green:

| Test | What it checks |
|---|---|
| `metal_gradient` | Analytic metal forces match finite differences. |
| `cn_detection` | The binary reports the expected CN targets for GlfT2 (Mn²⁺) and 1s63 (Zn²⁺). |

Additionally, changes that touch scoring, parsing or initialization must not
break the no-metal regression:

    tests/test_regression_no_metal.sh <path-to-vinametal>/vina <path-to-vina-1.2.7>

The script compares `--local_only` output byte for byte. Run it whenever your
change can affect the standard Vina terms.

If you add a new metal or change metal parameters, add or update a validation
log under `docs/validation/` with the exact command line, seed, box and
resulting RMSD. Do not edit existing logs to match new numbers; add a new log
and explain the change in the pull request.

## Coding Style

- C++11. Do not introduce a newer standard without discussion.
- Match the style of the surrounding files: 4-space indentation, `fl` for
  floating point, `sz`/`szv` for indices and index vectors, `VINA_FOR` /
  `VINA_FOR_IN` for loops.
- Prefer analytic derivatives over numerical ones. If you add an energy term,
  implement its gradient in the same function and extend
  `tests/metal_gradient_test.cpp`.
- Never silently change the behavior of the standard Vina scoring terms when
  metal scoring is disabled.
- Keep functions focused. Comments should explain non-obvious physics or
  numerical choices, not restate the code.

## Pull Requests

1. Fork the repository and create a topic branch from `main`.
2. Keep commits focused and write imperative commit messages
   (`Add Mg2+ cutoff`, not `added stuff`).
3. Describe *why* the change is needed, not only what it does.
4. State which tests you ran and paste the relevant output.
5. Update `CHANGELOG.md` under `[Unreleased]` for user-visible changes.
6. If the change affects input formats, options or validation numbers, update
   the corresponding file in `docs/`.

By submitting a pull request you agree that your contribution is licensed under
the Apache License 2.0.

## Binaries

Do not commit compiled binaries, object files, build directories or docking
outputs. Platform binaries are produced by the CI workflow and attached to
GitHub Releases. The repository must remain source-only.

## Scientific Claims

Phrase performance claims precisely. In particular:

- Metal parameters that have not been calibrated must be described as
  placeholders, both in code and documentation.
- Report the exact protocol (box, seed, exhaustiveness, receptor provenance) for
  any validation number you add.
- Do not present preliminary results as a benchmark. A benchmark requires a
  defined, documented set of complexes and a reproducible protocol.
