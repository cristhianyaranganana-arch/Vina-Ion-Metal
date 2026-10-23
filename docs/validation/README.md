# Validation logs

These logs are the raw output of the runs summarized in the repository README.
Unless noted otherwise they were produced with the VinaMetal 1.0.0 binary
(predecessor of Vina-Ion-Metal) before the CMake build was added to the
repository; the scoring code is the same. Commands, seeds and box sizes are
recorded inside each log.

| File | System | Result |
|---|---|---|
| `glft2_redock.log` | GlfT2, PDB 4FIY chain A, Mn²⁺ — **v1.0.2 crystal redocking** (bundled ligand = deposited pose) | CN target 2 (4 receptor donors); top-ranked pose −8.861 kcal/mol (metal-aware) at 1.50 Å symmetry-aware RMSD from the crystal UDP, coordinating Mn²⁺ with three oxygens at 2.12/2.18/2.51 Å; all 9 modes coordinate Mn²⁺ (min O–Mn 2.10–2.62 Å). |
| `glft2_cnfix.log` | GlfT2, PDB 4FIY chain A, Mn²⁺ — v1.0.0-era run, **translated ligand input** (superseded by the row above) | CN target 2 (4 receptor donors); top modes coordinate Mn²⁺ at Harding-consistent O–Mn distances (~2.1–2.2 Å). Affinities predate the F6 fix (no metal term in the reported energy) and the input was not the crystal pose. |
| `1s63_cnfix.log` | PDB 1s63, Zn²⁺ | CN target 1 (3 receptor donors); top-1 pose has 0.33 Å RMSD. Poses are unchanged in v1.0.2 (F6 affects reported energies only). |
| `regression_no_metal.log` | 1iep, no metal | `--local_only` output byte-identical to AutoDock Vina 1.2.7. |

## Notes

- **Version caveat.** The metalloprotein numbers in the v1.0.0-era logs
  (`*_cnfix.log`) were obtained with VinaMetal 1.0.0 defaults (always-strict
  CN penalty, weights 1.0/1.0/1.0). Vina-Ion-Metal 1.0.0 changes the CN
  regime to lenient-by-default and the weights to 1.2/0.5/1.0, so re-running
  `examples/*/config.txt` may produce different rankings. The no-metal
  regression is unaffected (CI-enforced, byte-identical against stock Vina
  1.2.7).
- **Energy-reporting caveat (v1.0.2, fix F6).** Before v1.0.2 the metal
  coordination term was cancelled out of every reported affinity (it entered
  the total with opposite signs via the intra and unbound buckets). The
  affinities printed in the pre-1.0.2 logs therefore do **not** include the
  metal chemistry and are not comparable to v1.0.2 affinities. Poses are
  unaffected: F6 changes only reported energies, and since v1.0.2 any
  printed mode rescored with `--score_only` reproduces the reported
  affinity (mode 1 of the GlfT2 crystal redocking to 0.002 kcal/mol).
- RMSD values are heavy-atom, symmetry-aware, computed with RDKit `GetBestRMS`
  against the crystal ligand. They are not printed by Vina and are reported in
  the log header comment or in the example READMEs.
- The GlfT2 and 1s63 v1.0.0-era logs show the `AutoDock Vina v1.2.7` banner
  because the binary predates the `v1.0.0` version string. The `VinaMetal:`
  lines in those logs confirm metal scoring was active.
- The `glft2_cnfix.log` run used a translated ligand input (closest heavy
  atom ~5.6 Å from the Mn²⁺), so crystal-pose RMSD is not a reproducible
  metric for that log; the CN target is determined by the receptor and is
  independent of that translation. **Since v1.0.2 the bundled GlfT2 ligand is
  the deposited 4FIY crystal pose** (O1A 2.89 Å / O1B 2.79 Å coordination
  preserved), turning the example into a legitimate redocking documented in
  `glft2_redock.log`.
- Global docking trajectories are not expected to match binaries built against
  a different Boost version, because Boost's normal distribution defines the
  random stream. Scores and gradients for a fixed pose are unaffected; this is
  why the regression test uses `--local_only`.
