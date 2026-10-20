# Changelog

All notable changes to VinaMetal are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

No changes yet.

## [1.0.0] - 2026-10-19

First public release. VinaMetal is a fork of AutoDock Vina 1.2.7 that adds
native transition-metal coordination scoring to the Vina scoring function.

### Added

- Metal-specific 12-6 van der Waals pair potential between each parametrized
  ion and ligand donor atoms (N/O/S), using AD4 parameters and Harding
  coordination distances.
- Smooth, differentiable coordination-number (CN) penalty driven by a sigmoid
  occupancy function, with analytic gradients for the BFGS minimizer.
- Bounded coordination-geometry bonus for tetrahedral, octahedral and
  square-planar sites, aggregated over donor pairs.
- Automatic receptor occupancy detection: before scoring, each metal counts the
  rigid receptor donors inside its coordination cutoff and reduces the ligand CN
  target to `max(0, CN_ideal - receptor-bound donors)`. The adjustment is
  reported at startup (`VinaMetal: <ion> CN target = ideal ...`).
- Built-in parameters for Zn²⁺, Mn²⁺, Mg²⁺, Ca²⁺, Fe, Cu²⁺ and Ni²⁺, with
  per-ion coordination cutoffs and ideal geometries.
- `--metals auto|off` switch to enable or disable metal coordination scoring.
- `--metal_params` JSON override for per-site ion parameters, including a
  site-specific `ideal_cn` for metals with unusual coordination.
- `--metal_oxidation` to declare oxidation states of multi-valent metals
  (for example `Fe=3`).
- `--weight_metal_vdw`, `--weight_metal_cn` and `--weight_metal_geom` to scale
  the three metal energy terms independently.
- `--metal_cn_strict` to penalize any deviation from the ideal coordination
  number, including incomplete shells.
- Stouten desolvation fields (`vol`, `solpar`) ported and exposed in the metal
  parameter set. The metal-specific desolvation term is not yet applied.
- CMake build (`vina`, `vina_split`, `vinametal_core` static library) with
  CTest integration.
- C++ gradient test comparing analytic metal forces against finite differences.
- End-to-end CN detection test on the two bundled examples.
- No-metal regression script comparing `--local_only` output against unmodified
  AutoDock Vina 1.2.7.
- Example systems: GlfT2 (PDB 4FIY chain A, Mn²⁺) and PDB 1s63 (Zn²⁺).
- Helper wrappers for Meeko ligand and receptor preparation.
- Documentation: installation, usage, limitations and validation logs.

### Validation

- **GlfT2 / PDB 4FIY chain A (Mn²⁺):** 4 receptor donors detected, ligand CN
  target 2; best pose RMSD 0.95 Å (mode 3 of 9).
- **PDB 1s63 (Zn²⁺):** 3 receptor donors detected, ligand CN target 1; top-1
  pose RMSD 0.33 Å.
- **No-metal regression:** byte-identical `--local_only` output against
  AutoDock Vina 1.2.7 on the 1iep test system.
- RMSD values are heavy-atom, symmetry-aware, computed with RDKit
  `GetBestRMS` against the crystal ligand.

### Known Limitations

See [docs/LIMITATIONS.md](docs/LIMITATIONS.md). In short: Fe oxidation states
are manual, Cu²⁺ and Ni²⁺ parameters are placeholders, the receptor is rigid,
occupancy counting uses Vina donor typing, and the metal weights are not yet
benchmark-calibrated.

### Release assets

- Linux x86-64 binary archive and examples archive are built and published by
  GitHub Actions on each `v*` tag.
- Prebuilt Windows binary is **deferred** for this release; Windows users can
  build from source.

### Licensing

- Source code: Apache License 2.0. The original Scripps Research Institute
  copyright notices are retained.
- `data/ad4_reference/`: GNU General Public License v2.0 (see
  `data/ad4_reference/COPYING`).

[Unreleased]: https://github.com/cristhianyaranganana-arch/VinaMetal/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/cristhianyaranganana-arch/VinaMetal/releases/tag/v1.0.0
