# Examples

Each example directory contains a receptor, a ligand and a Vina configuration
file. Paths in `config.txt` are relative to the working directory, so run each
example from inside its own directory:

    cd glft2
    ../../vina --config config.txt --out glft2_out.pdbqt

| Directory | System | Ion | CN target | Validated RMSD (VinaMetal 1.0.0) |
|---|---|---|---|---|
| `glft2/` | GlfT2, PDB 4FIY chain A | Mn²⁺ | 2 | 0.95 Å (mode 3) |
| `1s63/` | PDB 1s63 | Zn²⁺ | 1 | 0.33 Å (top-1) |

Both configurations set `seed = 42` and `metals = auto`. The CN target printed
at startup is computed from the receptor and is not hard-coded. The RMSD
values above were obtained with VinaMetal 1.0.0 defaults; Vina-Ion-Metal 1.0.0
changes the CN regime (lenient default) and the default weights, so rankings
may shift — run the examples to obtain fresh numbers with the current
defaults.

The inputs are provided for testing and demonstration. See the repository
`README.md` for citation and licensing of the PDB-derived structures.
