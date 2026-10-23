# GlfT2 (PDB 4FIY, chain A) — Mn2+ site

Crystal-pose redocking example for the galactofuranosyltransferase GlfT2 from
*Mycobacterium tuberculosis*.

## Inputs

| File | Description |
|---|---|
| `receptor.pdbqt` | Chain A of PDB entry **4FIY**, including the Mn2+ ion (chain A, residue 702). Prepared as a rigid receptor. |
| `ligand.pdbqt` | UDP (uridine-5'-diphosphate) **at the deposited 4FIY crystal pose** (chain A, residue 701), prepared with Meeko: the heavy-atom coordinates are preserved exactly, and the two oxygens that coordinate the Mn2+ in the crystal (O1B at 2.79 Å, O1A at 2.89 Å) are typed as the unprotonated (=O) sites. |
| `config.txt` | Search box covering the Mn2+/UDP site, 22 Å per side (center in the file), seed 42, exhaustiveness 16, 9 output modes. |

The second Mn2+ of 4FIY belongs to chain B and is not part of this receptor.
Before v1.0.2 this example shipped a translated input instead of the crystal
pose, so crystal-pose RMSD was not a reproducible metric (see EJ1 in the
v1.0.2 CHANGELOG).

## Expected behavior

Vina-Ion-Metal detects the Mn2+ automatically and counts the rigid receptor donors
already bound to it (Asp256, His396 and two Asp258 carboxylate oxygens; four
donors within the 3.1 Å coordination cutoff). The effective ligand CN target is
therefore `6 - 4 = 2`, and the log reports:

    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

In the recorded validation run (seed 42, exhaustiveness 16, Vina-Ion-Metal
v1.0.2 with static Boost 1.83; see `docs/validation/glft2_redock.log`) the
**top-ranked pose scores −8.86 kcal/mol** (metal-aware, see below), sits
**1.50 Å** (symmetry-aware heavy-atom RMSD) from the deposited UDP pose, and
coordinates the Mn2+ with three oxygens at 2.12/2.18/2.51 Å, consistent with
Harding Mn–O distances (~2.15 Å). All nine output modes coordinate the Mn2+
(minimum O–Mn 2.10–2.62 Å); the closest-to-crystal mode found (0.64 Å RMSD)
ranks 9th. Global-search trajectories depend on the Boost version, so a
different build may reorder the modes while keeping the metal distances.

Since v1.0.2 (fix F6) the reported affinities include the metal coordination
term; poses are unchanged relative to v1.0.1, and rescoring any printed mode
with `--score_only` reproduces the reported affinity (mode 1 above to
0.002 kcal/mol).

## Run

    vina --config config.txt --out glft2_out.pdbqt

## Source

- PDB entry: **4FIY** (RCSB Protein Data Bank).
- Ligand residue: UDP (chain A, residue 701).
- Metal: MN (chain A, residue 702).
