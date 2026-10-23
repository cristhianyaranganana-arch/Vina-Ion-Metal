# GlfT2 (PDB 4FIY, chain A) — Mn2+ site

Redocking example for the galactofuranosyltransferase GlfT2 from
*Mycobacterium tuberculosis*.

## Inputs

| File | Description |
|---|---|
| `receptor.pdbqt` | Chain A of PDB entry **4FIY**, including the Mn2+ ion (chain A, residue 702). Prepared as a rigid receptor. |
| `ligand.pdbqt` | UDP (uridine-5'-diphosphate) derived from the same PDB entry, prepared with Meeko. **It is a translated input, not the crystal pose**: its closest heavy atom sits ~5.6 Å from the Mn2+, whereas the deposited 4FIY ligand coordinates the ion at 2.79/2.89 Å. |
| `config.txt` | Search box centered on the Mn2+ site, 22 Å per side, seed 42, exhaustiveness 16, 9 output modes. |

The second Mn2+ of 4FIY belongs to chain B and is not part of this receptor.

## Expected behavior

Vina-Ion-Metal detects the Mn2+ automatically and counts the rigid receptor donors
already bound to it (Asp256, His396 and two Asp258 carboxylate oxygens; four
donors within the 3.1 Å coordination cutoff). The effective ligand CN target is
therefore `6 - 4 = 2`, and the log reports:

    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

In a validation run (seed 42, exhaustiveness 16) the top modes coordinate the
Mn2+ with phosphate oxygens at ~2.1–2.2 Å, consistent with Harding Mn–O
distances (~2.15 Å). **No pose is expected to reproduce the crystal UDP pose**:
the bundled ligand is not the crystal pose (see Inputs), so this example
demonstrates metal-guided pose generation, not crystal-pose redocking;
crystal-pose RMSD is not a reproducible metric for this example.

## Run

    vina --config config.txt --out glft2_out.pdbqt

## Source

- PDB entry: **4FIY** (RCSB Protein Data Bank).
- Ligand residue: UDP (chain A, residue 701).
- Metal: MN (chain A, residue 702).
