# Vina-Ion-Metal

**English** | [Español](README.es.md)

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![CI](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/actions/workflows/build.yml/badge.svg)](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/actions/workflows/build.yml)
[![Release](https://img.shields.io/github/v/release/cristhianyaranganana-arch/Vina-Ion-Metal)](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/releases)
![Platform: Linux x64](https://img.shields.io/badge/platform-Linux%20x64-blue.svg)

**AutoDock Vina fork with native, physics-based transition-metal coordination
scoring, automatic binding-site occupancy detection, and a robust receptor
preparation pipeline.**

Vina-Ion-Metal 1.0.0 is the corrected and completed successor of
[VinaMetal 1.0.0](https://github.com/cristhianyaranganana-arch/VinaMetal). It
is built on AutoDock Vina 1.2.7 and incorporates the fixes for the 16 problems
identified by the 2026 technical audit of VinaMetal, plus additional defects
found during this revision (including a corrupted CI trigger that silently
disabled continuous integration).

> **What it is for:** docking into metalloenzyme active sites where ligand-to-metal
> coordination is part of the binding mechanism (e.g. GlfT2/Mn²⁺ in
> *M. tuberculosis*, metallo-β-lactamases, carbonic anhydrases, HDACs).
> Vina-Ion-Metal is a **mechanistic filter that guides poses toward correct
> metal coordination** — it is not an affinity predictor, and its scoring
> weights are not yet benchmark-calibrated (see
> [Limitations](#limitations)).

## The Problem

Generic empirical scoring functions treat metal ions as ordinary heavy atoms.
They cannot reward the correct coordination number, penalize over-coordination,
or distinguish tetrahedral from octahedral sites. Stock Vina places the
generic metal type at an optimal contact distance of ~2.9 Å, far from real
coordination distances (~2.0–2.4 Å), so docking campaigns against
metalloenzymes systematically miss the catalytic geometry.

## The Solution

Vina-Ion-Metal integrates AutoDock 4 metal physics directly into Vina's native
scoring function, with analytic gradients consumed by the unmodified ILM+BFGS
engine:

- **Metal-specific 12-6 van der Waals** with Harding coordination distances
  (per ion and donor element N/O/S), taken from the AD4 parameter tables.
- **Coordination-number (CN) penalty** with two regimes:
  - *lenient* (default): only over-coordination is penalized (smooth C-infinity
    hinge) — an incomplete ligand shell is not punished;
  - *strict* (`--metal_cn_strict`): any deviation from the CN target is
    penalized with the pure quadratic.
- **Bounded coordination-geometry bonus** (tetrahedral, octahedral,
  square-planar).
- **Automatic receptor occupancy detection**: rigid-receptor donors already
  bound to each metal are counted, and the ligand CN target becomes
  `max(0, CN_ideal − receptor-bound donors)` — no per-site configuration.
  A typing-independent geometric donor count is also reported as a
  preparation-quality diagnostic when the two counts disagree.
- **Robust receptor preparation** (`scripts/prepare_receptor.py`): altLoc
  resolution, chain-integrity checks, metal typing **including Cu²⁺/Ni²⁺
  (which Meeko cannot emit)**, and a coordination pre-flight report.

Metal–donor pairs involving parametrized ions are excluded from the standard
Vina potentials to avoid double counting. Receptors containing
**unparametrized** metals (Na⁺, K⁺, Co, Cd, Hg, …) keep the exact stock Vina
1.2.7 treatment, preserving the byte-identical promise.

## Validation

| System | Ion | Receptor donors detected | Ligand CN target | Result (recorded validation runs, seed 42, exhaustiveness 16) |
|---|---|---|---|---|
| GlfT2, PDB 4FIY chain A | Mn²⁺ | 4 | 2 | **v1.0.2 crystal redocking**: top-ranked pose −8.86 kcal/mol (metal-aware), 1.50 Å symmetry-aware RMSD against the deposited UDP pose, coordinates the Mn²⁺ with three oxygens at 2.12/2.18/2.51 Å (Harding Mn–O = 2.15 Å); all 9 modes coordinate the ion |
| PDB 1s63 | Zn²⁺ | 3 (geometric 4, diagnostic printed) | 1 | top-1 pose coordinates the Zn²⁺ at 2.26 Å (poses unchanged since the v1.0.0-era run) |
| PDB 1iep (no metal) | — | — | — | `--local_only` output **byte-identical** to AutoDock Vina 1.2.7 (CI-enforced) |

The automatic adjustment is printed before docking, for example:

    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

Pose-level metal distances are measured from the output PDBQT; RMSD should be
computed symmetry-aware (RDKit `GetBestRMS`). With respect to VinaMetal 1.0.0
the CN regime (lenient default) and default weights changed (audits #2 and
#9), so rankings differ; run `examples/*/config.txt` to reproduce. Since
v1.0.2 the GlfT2 example bundles the deposited crystal ligand, so its
redocking is measured by crystal-pose RMSD (see `docs/validation/`), and the
reported affinities include the metal term (fix F6).

## Quick Start

Requirements: CMake ≥ 3.12, a C++11 compiler, and Boost ≥ 1.65 with the
`system`, `thread`, `serialization`, `filesystem` and `program_options`
libraries. Static Boost is preferred by default (portable binaries); pass
`-DBoost_USE_STATIC_LIBS=OFF` if your system only ships dynamic libraries.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

Prepare a receptor (metals typed automatically, Cu/Ni supported):

    scripts/prepare_receptor.py --pdb receptor.pdb --out receptor.pdbqt

Run the GlfT2 example (the Mn²⁺ is detected automatically):

    cd examples/glft2
    ../../build/vina --config config.txt --metals auto --out glft2_out.pdbqt

The log prints the occupancy-adjusted target:

    Vina-Ion-Metal: parametrized metal ions in receptor: Mn(x1) - coordination scoring is ACTIVE
    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

Use `--metals off` to disable metal coordination scoring and recover plain
Vina behavior. See [docs/USAGE.md](docs/USAGE.md) for all options (including
the recommended multi-seed protocol), [docs/INSTALL.md](docs/INSTALL.md) for
platform-specific builds, and [docs/LIMITATIONS.md](docs/LIMITATIONS.md)
before production use.

## What was fixed relative to VinaMetal 1.0.0

The complete audit trail lives in [CHANGELOG.md](CHANGELOG.md). Summary of the
16 audited problems:

| # | Problem | Status in Vina-Ion-Metal 1.0.0 |
|---|---|---|
| 1 | Grid over-exclusion for any `Met`-element atom | **Fixed** — only parametrized ions are excluded; regression test added |
| 2 | Dead option `--metal_cn_strict` | **Implemented** — lenient (default) vs strict CN penalty, C-infinity smooth |
| 3 | Dead variable `dE_dcos` | **Removed** |
| 4 | Branding pointed to upstream Vina | **Fixed** — version string, issue tracker and citations are Vina-Ion-Metal |
| 5 | Cu²⁺/Ni²⁺ unreachable via Meeko | **Fixed** — `prepare_receptor.py` types them natively |
| 6 | `--metal_oxidation` silent failures | **Fixed** — strict validation (element, integer, range) |
| 7 | Silent oxidation-state fallback | **Fixed** — warning with the parametrized states at receptor load |
| 8 | Orphan `data/AD4Zn.dat` | **Removed** (it was a mislabeled copy of `AD4.1_bound.dat`) |
| 9 | Uncalibrated weights (1.0/1.0/1.0) | **Heuristic defaults** (1.2/0.5/1.0) with runtime notice; calibration benchmark still pending |
| 10 | Fragile occupancy detection | **Hardened** — geometric donor count + divergent-count warnings |
| 11 | Fragile receptor preparation | **Replaced** — own robust pipeline (altLoc, integrity checks, pre-flight report) |
| 12 | Non-portable dynamic-Boost binaries | **Fixed** — static Boost by default (with graceful fallback) |
| 13 | Regression test not runnable out of the box | **Fixed** — registered in CTest; CI downloads the Vina 1.2.7 reference |
| 14 | Linux-only CI | **Kept** (intentional scope); CI trigger bug found and fixed |
| 15 | RMSD by atom order (auditor's benchmark) | **Documented** — recommended protocol uses symmetry-aware RMSD |
| 16 | Single seed / low exhaustiveness | **Documented** — recommended protocol (≥3 seeds, exhaustiveness ≥16) |

Additional fixes not in the audit: corrupted `branches: ain]` YAML in the CI
workflow (CI never triggered), and the AutoDock4Zn citation in the README
carried wrong page numbers (see below).

## Limitations

- Metal scoring weights are **provisional heuristics**, not calibrated against
  a multi-protein benchmark; treat Vina-Ion-Metal as a mechanistic pose filter,
  not an affinity predictor.
- Fe²⁺/Fe³⁺ share the generic AD4 iron parameters; oxidation-state-specific
  refinement is pending.
- Cu²⁺/Ni²⁺ parameters remain AD4-style estimates (not calibrated in AD4);
  override with `--metal_params` for serious work.
- The receptor is rigid; coordinating side chains do not reorganize.
- Pregenerated maps must come from Vina-Ion-Metal itself (metal ions are
  excluded from donor-type maps).

See [docs/LIMITATIONS.md](docs/LIMITATIONS.md) for the complete list.

## Recommended docking protocol

1. Prepare the receptor with `scripts/prepare_receptor.py`; read the
   coordination pre-flight report and resolve any warnings.
2. Dock with `--exhaustiveness 16` (or higher) and **at least 3 different
   seeds**; per-seed rankings with low exhaustiveness are noise-dominated.
3. Evaluate poses by symmetry-aware heavy-atom RMSD (e.g. RDKit
   `GetBestRMS`), never by atom order.
4. Treat metal-coordinating hits as **hypothesis generators**; re-rank with
   MM/GBSA or short MD, and apply PAINS/aggregation counter-screens before any
   biochemical assay.

## Citation

    @software{vina_ion_metal2026,
      title        = {Vina-Ion-Metal: AutoDock Vina fork with native
                      transition-metal coordination scoring},
      author       = {Yaranga Nana, Christian},
      year         = {2026},
      version      = {1.0.0},
      url          = {https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal},
      note         = {Corrected successor of VinaMetal 1.0.0; based on
                      AutoDock Vina 1.2.7}
    }

This work builds on VinaMetal 1.0.0 (Yaranga Nana, 2026,
DOI 10.5281/zenodo.23499212) — please cite it when comparing results across
versions. Please also cite the foundations:

- Trott, O.; Olson, A. J. *J. Comput. Chem.* **2010**, 31, 455–461.
- Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S.
  *J. Chem. Inf. Model.* **2021**, 61, 3891–3898.
- Santos-Martins, D.; Forli, S.; Ramos, M. J.; Olson, A. J. AutoDock4Zn: an
  improved AutoDock force field for small-molecule docking to zinc
  metalloproteins. *J. Chem. Inf. Model.* **2014**, 54(8), 2371–2379.
- Harding, M. M. *Acta Crystallogr. D* **2001**, 57, 401–411.

## Credits & License

Vina-Ion-Metal is developed by **Christian Yaranga Nana**
(Universidad Norbert Wiener, Lima, Peru) as the corrected successor of
VinaMetal 1.0.0. It is built on AutoDock Vina and AutoDock 4, developed at
The Scripps Research Institute; their original copyright notices are retained
in the source files (see [NOTICE](NOTICE)). The example inputs are derived
from PDB entries **4FIY**, **1s63** and **1iep** (RCSB Protein Data Bank).

The source code is distributed under the Apache License, Version 2.0; see
[LICENSE](LICENSE). AutoDock 4 parameter files in `data/ad4_reference/` are
distributed under the GNU General Public License v2.0; see
`data/ad4_reference/COPYING`.
