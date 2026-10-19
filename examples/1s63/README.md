# PDB 1s63 — Zn2+ site

Redocking example for a zinc metalloprotein.

## Inputs

| File | Description |
|---|---|
| `receptor.pdbqt` | Receptor from PDB entry **1s63**, including its Zn2+ ion. This is the `protein.pdbqt` file distributed with the AutoDock Vina example set. |
| `ligand.pdbqt` | Ligand extracted from PDB entry **1s63**, prepared with Meeko. |
| `config.txt` | Search box centered at (18, 134, 0), seed 42, exhaustiveness 8, 9 output modes. |

## Expected behavior

VinaMetal detects the Zn2+ automatically and counts the rigid receptor donors
already bound to it (three donors within the 3.0 Å coordination cutoff). The
effective ligand CN target is therefore `4 - 3 = 1`, and the log reports:

    VinaMetal: Zn2+ CN target = ideal 4 - receptor-bound 3 = 1

In the validation run, the top-ranked pose had a heavy-atom symmetry-aware RMSD
of 0.33 Å against the crystal ligand.

## Run

    vina --config config.txt --out 1s63_out.pdbqt

## Source

- PDB entry: **1s63** (RCSB Protein Data Bank).
- Receptor and ligand coordinates: PDB entry 1s63; receptor PDBQT from the
  AutoDock Vina example set.
