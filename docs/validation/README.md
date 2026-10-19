# Validation logs

These logs are the raw output of the runs summarized in the repository README.
They were produced with the VinaMetal 1.0.0 binary before the CMake build was
added to the repository; the scoring code is the same. Commands, seeds and box
sizes are recorded inside each log.

| File | System | Result |
|---|---|---|
| `glft2_cnfix.log` | GlfT2, PDB 4FIY chain A, Mn²⁺ | CN target 2 (4 receptor donors); mode 3 reproduces the crystal pose with 0.95 Å RMSD. |
| `1s63_cnfix.log` | PDB 1s63, Zn²⁺ | CN target 1 (3 receptor donors); top-1 pose has 0.33 Å RMSD. |
| `regression_no_metal.log` | 1iep, no metal | `--local_only` output byte-identical to AutoDock Vina 1.2.7. |

## Notes

- RMSD values are heavy-atom, symmetry-aware, computed with RDKit `GetBestRMS`
  against the crystal ligand. They are not printed by Vina and are reported in
  the log header comment or in the example READMEs.
- The GlfT2 and 1s63 logs show the `AutoDock Vina v1.2.7` banner because the
  binary predates the `v1.0.0` version string. The `VinaMetal:` lines in those
  logs confirm metal scoring was active.
- The bundled GlfT2 ligand was translated into the crystal frame so that it lies
  inside the docking box. The CN target is determined by the receptor and is
  independent of that translation; global-search mode ordering can differ with
  the Boost version (see the regression note below).
- Global docking trajectories are not expected to match binaries built against
  a different Boost version, because Boost's normal distribution defines the
  random stream. Scores and gradients for a fixed pose are unaffected; this is
  why the regression test uses `--local_only`.
