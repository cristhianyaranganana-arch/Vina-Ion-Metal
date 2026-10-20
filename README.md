# VinaMetal

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![CI](https://github.com/cristhianyaranganana-arch/VinaMetal/actions/workflows/build.yml/badge.svg)](https://github.com/cristhianyaranganana-arch/VinaMetal/actions/workflows/build.yml)
[![Release](https://img.shields.io/github/v/release/cristhianyaranganana-arch/VinaMetal)](https://github.com/cristhianyaranganana-arch/VinaMetal/releases)
![Platform: Linux x64](https://img.shields.io/badge/platform-Linux%20x64-blue.svg)

**AutoDock Vina fork with native, physics-based transition metal coordination
and automatic binding site occupancy detection.**

VinaMetal 1.0.0 is a fork of AutoDock Vina 1.2.7 developed for reproducible
docking into metalloprotein active sites.

> **Platform support:** the prebuilt release assets are **Linux x86-64 only**.
> Windows users can still build from source; a prebuilt Windows binary is
> planned for a future release.

## The Problem

Generic empirical scoring functions treat metal ions as ordinary heavy atoms.
As a result, they cannot reward the correct metal–ligand coordination number,
penalize over-coordination, or distinguish tetrahedral from octahedral sites.
Docking campaigns against metalloenzymes therefore depend on manual
intervention and are difficult to reproduce.

## The Solution

VinaMetal integrates AutoDock 4 metal physics directly into Vina's native
scoring function, with analytic gradients consumed by the unmodified ILM+BFGS
engine:

- Metal-specific 12-6 van der Waals parameters and Harding coordination
  distances, taken from the AD4 parameter tables.
- A smooth, differentiable coordination-number (CN) penalty.
- A bounded coordination-geometry bonus (tetrahedral, octahedral,
  square-planar).
- **Automatic Receptor Occupancy Detection**: before scoring, the rigid
  receptor donors already bound to each metal are counted within the ion's
  coordination cutoff, and the ligand CN target becomes
  `max(0, CN_ideal − receptor-bound donors)`. No manual per-site override is
  required.

Metal–donor pairs are excluded from the standard Vina Gaussian term to avoid
double counting. The Stouten desolvation parameters of AD4 are ported but
currently inactive; see [Limitations](#limitations).

## Validation

| System | Ion | Receptor donors detected | Ligand CN target | Best-pose RMSD |
|---|---|---|---|---|
| GlfT2, PDB 4FIY chain A | Mn²⁺ | 4 | 2 | **0.95 Å** (mode 3 of 9) |
| PDB 1s63 | Zn²⁺ | 3 | 1 | **0.33 Å** (top-1) |

RMSD is heavy-atom, symmetry-aware (RDKit `GetBestRMS`) against the crystal
ligand. The automatic adjustment was confirmed in the run logs, for example:

    VinaMetal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

For receptors without parametrized metal ions, the scoring and gradient path is
byte-identical to unmodified AutoDock Vina 1.2.7 (`--local_only` comparison);
see `docs/validation/`.

## Quick Start

Requirements: CMake ≥ 3.12, a C++11 compiler, and Boost ≥ 1.65 with the
`system`, `thread`, `serialization`, `filesystem` and `program_options`
libraries.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

Run the GlfT2 example (the Mn²⁺ is detected automatically):

    cd examples/glft2
    ../../build/vina --config config.txt --metals auto --out glft2_out.pdbqt

The log prints the occupancy-adjusted target:

    VinaMetal: parametrized metal ions in receptor: Mn(x1) - coordination scoring is ACTIVE
    VinaMetal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

Use `--metals off` to disable metal coordination scoring and recover plain Vina
behavior. See [docs/USAGE.md](docs/USAGE.md) for all options,
[docs/INSTALL.md](docs/INSTALL.md) for platform-specific builds, and
[docs/LIMITATIONS.md](docs/LIMITATIONS.md) before production use.

## Limitations

- Fe oxidation states must be declared manually (`--metal_oxidation Fe=3`).
- Cu²⁺ and Ni²⁺ parameters are placeholders pending calibration.
- The receptor is rigid; coordinating side chains do not reorganize.
- Receptor occupancy uses a hard integer cutoff and Vina's donor typing; flex
  residues are not counted as occupancy.
- The metal scoring weights are not yet benchmark-calibrated.

See [docs/LIMITATIONS.md](docs/LIMITATIONS.md) for the complete list.

## Citation

    @software{vinametal2026,
      title        = {VinaMetal},
      author       = {Christian Yaranga Nana},
      year         = {2026},
      version      = {1.0.0},
      doi          = {Pending (Will be updated via Zenodo)},
      url          = {https://github.com/cristhianyaranganana-arch/VinaMetal},
      note         = {Fork of AutoDock Vina 1.2.7 with native metal coordination}
    }

Please also cite the foundations of this work:

- Trott, O.; Olson, A. J. *J. Comput. Chem.* **2010**, 31, 455–461.
- Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S.
  *J. Chem. Inf. Model.* **2021**, 61, 3891–3898.
- Santos-Martins, D.; Forli, S.; Ribeiro, A. J. M.; Ramos, M. J.; Olson, A. J.
  *J. Chem. Inf. Model.* **2014**, 54, 1931–1941.
- Harding, M. M. *Acta Crystallogr. D* **2001**, 57, 401–411.

## Credits & License

VinaMetal is built on AutoDock Vina and AutoDock 4, developed at The Scripps
Research Institute. The original copyright notices are retained in the source
files. The example inputs are derived from PDB entries **4FIY** and **1s63**
(RCSB Protein Data Bank).

The source code is distributed under the Apache License, Version 2.0; see
[LICENSE](LICENSE). AutoDock 4 parameter files in `data/ad4_reference/` are
distributed under the GNU General Public License v2.0; see
`data/ad4_reference/COPYING`.
