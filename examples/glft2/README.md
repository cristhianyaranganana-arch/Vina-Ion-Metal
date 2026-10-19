# GlfT2 (PDB 4FIY, chain A) — Mn2+ site

Redocking example for the galactofuranosyltransferase GlfT2 from
*Mycobacterium tuberculosis*.

## Inputs

| File | Description |
|---|---|
| `receptor.pdbqt` | Chain A of PDB entry **4FIY**, including the Mn2+ ion (chain A, residue 702). Prepared as a rigid receptor. |
| `ligand.pdbqt` | UDP (uridine-5'-diphosphate) extracted from the same PDB entry, prepared with Meeko. Its coordinates lie inside the docking box. |
| `config.txt` | Search box centered on the Mn2+ site, 22 Å per side, seed 42, exhaustiveness 16, 9 output modes. |

The second Mn2+ of 4FIY belongs to chain B and is not part of this receptor.

## Expected behavior

VinaMetal detects the Mn2+ automatically and counts the rigid receptor donors
already bound to it (Asp256, His396 and two Asp258 carboxylate oxygens; four
donors within the 3.1 Å coordination cutoff). The effective ligand CN target is
therefore `6 - 4 = 2`, and the log reports:

    VinaMetal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

In the validation run, the best pose (mode 3 of 9) reproduces the crystal
coordination with two phosphate oxygens at 2.11 and 2.14 Å from the Mn2+ and a
heavy-atom symmetry-aware RMSD of 0.95 Å against the crystal UDP.

## Run

    vina --config config.txt --out glft2_out.pdbqt

## Source

- PDB entry: **4FIY** (RCSB Protein Data Bank).
- Ligand residue: UDP (chain A, residue 701).
- Metal: MN (chain A, residue 702).
