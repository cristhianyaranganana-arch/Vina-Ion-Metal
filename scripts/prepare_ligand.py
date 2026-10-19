#!/usr/bin/env python3
"""Wrapper around Meeko's mk_prepare_ligand.py.

VinaMetal consumes PDBQT files; Meeko converts SDF/PDB ligands to PDBQT with
torsions and AutoDock atom types. This wrapper only forwards arguments, so the
upstream tool's options and behavior are unchanged.

Install Meeko first (it also provides mk_prepare_ligand.py), then run:

    scripts/prepare_ligand.py --mol ligand.sdf --out ligand.pdbqt

See https://github.com/forlilab/Meeko for the full option list.
"""

import shutil
import subprocess
import sys

TOOL = "mk_prepare_ligand.py"


def main():
    executable = shutil.which(TOOL)
    if executable is None:
        print(
            f"error: '{TOOL}' was not found on PATH. Install Meeko before "
            "preparing ligands.",
            file=sys.stderr,
        )
        return 127
    return subprocess.call([executable, *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
