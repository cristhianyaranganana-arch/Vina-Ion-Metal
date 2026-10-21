# Usage

Vina-Ion-Metal accepts the AutoDock Vina command line. This page documents the
metal coordination options, the CN adjustment, the scoring weights, and the
recommended docking protocol.

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

## Receptor preparation

Prepare receptors with the bundled pipeline (audit #5/#11):

    scripts/prepare_receptor.py --pdb receptor.pdb --out receptor.pdbqt

It resolves alternate locations, checks chain integrity (numbering gaps,
missing backbone atoms), types **all** parametrized ions — including Cu²⁺/Ni²⁺,
which Meeko cannot emit — keeps unparametrized metals (Na⁺, K⁺, …) with the
standard Vina treatment, and prints a **coordination pre-flight report**:
the geometric N/O/S donor count around each metal, the estimated ligand CN
target, and warnings for empty coordination spheres.

Two backends: `meeko` (default when installed; explicit hydrogens and charges)
and `simple` (no dependencies; united-atom typing from amino-acid templates,
valid for the vina/vinardo scoring functions — not for AD4 map generation).
Useful flags: `--keep-waters`, `--remove-het`, `--altloc A`,
`--metal-charge "Cu=2,Fe=3"`, `--meeko-args "..."`.

## Automatic receptor occupancy detection

Before scoring, Vina-Ion-Metal inspects the rigid receptor (`grid_atoms`) and
counts the N/O/S atoms typed as coordination donors that lie within the ion's
coordination cutoff. The ligand CN target is then

    CN_target = max(0, CN_ideal - CN_receptor_bound)

The adjustment is computed once per receptor and is therefore thread-safe.
With `-v` (the default), the result is printed before docking:

    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

Coordination cutoffs are ion-specific (for example 3.0 Å for Zn²⁺, 3.1 Å for
Mn²⁺, 3.3 Å for Ca²⁺; see `src/lib/metal_parameters.cpp`).

**Preparation-quality diagnostics (audit #10).** A second, typing-independent
geometric count of N/O/S atoms is computed inside the same cutoff. When the
two counts disagree you will see, before docking:

- `WARNING: <ion> has N N/O/S receptor atom(s) within <cutoff> A but NONE is
  typed as a coordination donor` — the CN target was **not** reduced; a fully
  coordinated ligand may be over-penalized. Re-prepare the receptor or check
  for missing residues/atoms near the metal.
- `note: <ion> geometric donor count (N) exceeds the typed count (M)` —
  non-donor-typed N/O/S atoms were ignored by the target computation.

## The CN penalty: lenient (default) vs strict

The coordination-number penalty compares the smooth count of movable donors
in the coordination shell with the CN target (audit #2):

- **Lenient (default).** Only *over*-coordination is penalized. The hinge is a
  C-infinity softplus approximation of `max(0, cn - target)`, so an
  incomplete ligand shell (monodentate ligand in a bidentate site, solvent
  accessibility, partial occupancy) is not punished. Recommended for virtual
  screening.
- **Strict (`--metal_cn_strict`).** Any deviation from the target — including
  under-coordination — is penalized with the pure quadratic
  `w_cn * (cn - target)^2`. Use it when the biological hypothesis requires the
  full coordination (e.g. metal chelation is the pharmacophore).

Note for VinaMetal users: VinaMetal 1.0.0 always applied the strict quadratic
(the `--metal_cn_strict` flag was dead), so results may shift with the new
lenient default.

## Metal options

| Option | Default | Description |
|---|---|---|
| `--metals auto\|off` | `auto` | Enable or disable metal coordination scoring. |
| `--metal_params FILE` | — | JSON file overriding built-in ion parameters (`ideal_cn`, cutoffs, radii, geometry, weights). |
| `--metal_oxidation Fe=3` | — | Declare oxidation states for multi-valent metals. Strictly validated (audit #6): unknown elements, non-integer values and states without parameters are rejected or warned. |
| `--weight_metal_vdw W` | 1.2 | Weight of the metal 12-6 pair term. |
| `--weight_metal_cn W` | 0.5 | Weight of the CN penalty. |
| `--weight_metal_geom W` | 1.0 | Weight of the geometry bonus. |
| `--metal_cn_strict` | off | Penalize any CN deviation, including incomplete ligand shells (default: lenient, over-coordination only). |

Oxidation-state handling (audits #6/#7):

    vina --metal_oxidation Fe=3            # use the Fe3+ entry
    vina --metal_oxidation Fe=abc          # ERROR: not a valid integer (hard failure)
    vina --metal_oxidation Xx=2            # ERROR: unknown metal element
    vina --metal_oxidation Fe=4            # WARNING at receptor load: no Fe +4
                                          # parameters; falling back (states: +2, +3)

## Metal scoring weights (audit #9)

The defaults are **provisional heuristics**, pending the full calibration
benchmark, and a notice is printed whenever coordination scoring is active:

    Vina-Ion-Metal: NOTE: default metal weights are provisional heuristics
    (v1.0: vdw 1.2, cn 0.5, geom 1.0); override with
    --weight_metal_vdw/--weight_metal_cn/--weight_metal_geom for calibrated
    campaigns.

Rationale of the v1.0 heuristics:

- `w_vdw = 1.2` — the AD4-derived pair well depth for a typical M–O donor is
  only ~0.3–0.4 kcal/mol (e.g. Zn–O 0.33, Mn–O 0.42); AutoDock4Zn recovers
  roughly twice that per donor through its ZA/TZ directional pseudoatoms. A
  +20 % scaling nudges the single-centered 12-6 term toward that scale while
  staying conservative.
- `w_cn = 0.5` — the receptor occupancy is a hard integer count and therefore
  the dominant error source of the squared deviation term (a miscounted site
  distorts the score by `(Δcn)^2 × w_cn`). Halving the weight halves that
  worst-case damage while still suppressing over-coordination; combined with
  the lenient default (under-coordination unpenalized), the systematic risk
  drops substantially.
- `w_geom = 1.0` — the geometry bonus is already bounded (max 1.0 per donor
  pair), so it cannot explode; it is left unchanged.

For virtual-screening campaigns, calibrate `(w_vdw, w_cn, w_geom)` against
your target class (re-docking RMSD and enrichment) and pass the values
explicitly.

## Recommended docking protocol (audits #15/#16)

Two methodological pitfalls were documented for the VinaMetal experiments:
RMSD computed by atom order (which symmetry makes meaningless) and rankings
from a single seed with low exhaustiveness (stochastic variance dominates the
1–2 kcal/mol differences of interest). The recommended protocol:

1. **Seeds and exhaustiveness.** Run at least 3 different seeds with
   `--exhaustiveness 16` or higher (the bundled examples use 16). A single
   seed at exhaustiveness 8 does not constrain the ranking.
2. **RMSD.** Evaluate poses with symmetry-aware heavy-atom RMSD (RDKit
   `GetBestRMS`, `obrms -f`), never by atom order.
3. **Interpretation.** Treat Vina-Ion-Metal as a mechanistic pose filter:
   select poses that coordinate the metal when coordination is a requirement
   of activity, then re-rank with MM/GBSA or short MD.
4. **Hit triage.** Metal-chelating chemotypes (catechols, flavonoids,
   hydroxypyridinones) are frequent PAINS and aggregation offenders; apply
   counter-screens (free-metal competition, detergent controls, orthogonal
   readouts) before prioritizing any hit.
