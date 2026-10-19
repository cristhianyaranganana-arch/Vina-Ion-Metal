# Usage

VinaMetal accepts the AutoDock Vina command line. This page documents the metal
coordination options and the automatic CN adjustment.

## Basic run

    vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
         --center_x 142.54 --center_y 32.19 --center_z 48.55 \
         --size_x 22 --size_y 22 --size_z 22 \
         --seed 42 --exhaustiveness 16 --num_modes 9 \
         --out results.pdbqt

or with a configuration file (`receptor`, `ligand`, `center_*`, `size_*`,
`seed`, `exhaustiveness`, `num_modes`, `metals`):

    cd examples/glft2
    vina --config config.txt --out results.pdbqt

Paths inside a configuration file are resolved relative to the current working
directory, not the configuration file location.

Metal coordination scoring is enabled by default (`--metals auto`). It is
activated only when the receptor contains a metal ion with built-in parameters.

## Automatic receptor occupancy detection

Before scoring, VinaMetal inspects the rigid receptor (`grid_atoms`) and counts
the N/O/S atoms typed as coordination donors that lie within the ion's
coordination cutoff. The ligand CN target is then

    CN_target = max(0, CN_ideal - CN_receptor_bound)

The adjustment is computed once per receptor and is therefore thread-safe.
With `-v` (the default), the result is printed before docking:

    VinaMetal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

The penalty compares this target with the smooth count of movable donors in the
coordination shell. A receptor with four protein donors therefore expects two
additional ligand donors, without any per-site configuration.

Coordination cutoffs are ion-specific (for example 3.0 Å for Zn²⁺, 3.1 Å for
Mn²⁺, 3.3 Å for Ca²⁺; see `src/lib/metal_parameters.cpp`).

## Metal options

| Option | Default | Description |
|---|---|---|
| `--metals auto\|off` | `auto` | Enable or disable metal coordination scoring. |
| `--metal_params FILE` | — | JSON file overriding built-in ion parameters (`ideal_cn`, cutoffs, radii, geometry, weights). |
| `--metal_oxidation Fe=3` | — | Declare oxidation states for multi-valent metals. |
| `--weight_metal_vdw W` | 1.0 | Weight of the metal 12-6 pair term. |
| `--weight_metal_cn W` | 1.0 | Weight of the CN penalty. |
| `--weight_metal_geom W` | 1.0 | Weight of the geometry bonus. |
| `--metal_cn_strict` | off | Penalize any CN deviation even with an incomplete shell. |

The weights are placeholders pending benchmark calibration; see
`docs/LIMITATIONS.md`.

### Parameter override example

```json
{
  "Mn2+": {
    "ideal_cn": 6,
    "cn_cutoff": 3.1,
    "geometry": "octahedral"
  }
}
```

Run with `--metal_params override.json`. An overridden `ideal_cn` is still
reduced by the automatic receptor occupancy.

## Scoring-only and local-only runs

`--score_only` and `--local_only` use the same receptor preparation and print
the CN target line, so they can be used to check the detected occupancy:

    vina --receptor examples/glft2/receptor.pdbqt \
         --ligand examples/glft2/ligand.pdbqt \
         --center_x 142.54 --center_y 32.19 --center_z 48.55 \
         --size_x 22 --size_y 22 --size_z 22 --score_only

## Interpreting results

- Coordination energy is integrated into the total affinity. Library users can
  obtain per-term values (vdW, CN, geometry) from
  `metal_scoring::eval_breakdown`; the CLI does not print a breakdown.
- The log line `CN target = ideal X - receptor-bound Y = Z` documents the
  occupancy adjustment for each metal.
- With a receptor containing no parametrized metal, VinaMetal is equivalent to
  AutoDock Vina 1.2.7; see `docs/validation/regression_no_metal.log`.

## Reproducibility checklist

Record the binary version (`vina --version`), the SHA-256 of the binary, the
input PDBQT hashes, the seed, exhaustiveness, box, and the full log. Binaries
built against different Boost versions may explore different trajectories for
the same seed while producing identical energies at a given pose; use
`--local_only` for scoring comparisons.
