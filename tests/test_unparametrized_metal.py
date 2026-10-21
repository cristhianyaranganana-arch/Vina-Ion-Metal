#!/usr/bin/env python3
"""Regression guard for the grid over-exclusion bug (audit problem #1).

VinaMetal 1.0.0 excluded from the donor-type affinity maps EVERY atom whose
element maps to the generic metal element (EL_TYPE_Met) - not only the
parametrized ions. A receptor containing an UNPARAMETRIZED metal (Na+, K+,
Co, Cd, Hg, ...) therefore lost the standard Vina interaction between that
metal and every ligand donor: with --metals auto the ligand saw a "ghost"
zero-interaction ion, breaking the byte-identical promise against stock Vina.

Vina-Ion-Metal only excludes PARAMETRIZED ions (metal_scoring::is_metal_ion),
so unparametrized metals keep the standard Vina map treatment.

This test appends a Na+ ion to the metal-free 1iep receptor and verifies that
--score_only energies are IDENTICAL between --metals auto and --metals off
(no parametrized ion present => the metal machinery must be a no-op).

Usage: tests/test_unparametrized_metal.py <vina-binary> <repo-root>
"""

import os
import subprocess
import sys
import tempfile


def read(path):
    with open(path) as f:
        return f.read()


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def score(vina, receptor, ligand, metals):
    cmd = [
        vina,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", "15.19", "--center_y", "53.903", "--center_z", "16.917",
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--metals", metals,
        "--score_only",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise AssertionError(
            f"vina exited with {result.returncode} (metals={metals})\n"
            + (result.stdout + result.stderr)[-2000:]
        )
    out = result.stdout + result.stderr
    # The "Estimated Free Energy of Binding" line carries the total energy.
    energy = None
    for line in out.splitlines():
        if "Estimated Free Energy of Binding" in line:
            energy = line.strip()
            break
    if energy is None:
        raise AssertionError(f"no energy line in output (metals={metals})\n{out[-2000:]}")
    return energy, out


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    vina, root = sys.argv[1], sys.argv[2]

    receptor = os.path.join(root, "tests", "data", "1iep_receptor.pdbqt")
    ligand = os.path.join(root, "tests", "data", "1iep_ligand.pdbqt")

    with tempfile.TemporaryDirectory() as tmp:
        # Append an unparametrized Na+ near (but not clashing with) the pocket.
        # Column layout mirrors the bundled receptor records (see parse_pdbqt):
        # x/y/z at 31-38/39-46/47-54, charge at 69-76, type at 78-79.
        base = read(receptor).rstrip("\n")
        na_line = (
            "ATOM   9999  NA  NA  A 900      15.500  53.900  17.500  "
            "1.00  0.00    +1.000 Na"
        )
        assert len(na_line) == 79, f"Na line width {len(na_line)} != 79"
        assert len(base.splitlines()) > 0
        na_receptor = os.path.join(tmp, "1iep_na.pdbqt")
        write(na_receptor, base + "\n" + na_line + "\n")

        e_auto, out_auto = score(vina, na_receptor, ligand, "auto")
        e_off, _ = score(vina, na_receptor, ligand, "off")

        print(f"metals=auto : {e_auto}")
        print(f"metals=off  : {e_off}")

        # 1) Energy must be IDENTICAL: no parametrized ion => no-op machinery.
        if e_auto != e_off:
            raise AssertionError(
                "UNPARAMETRIZED METAL CHANGED THE SCORE: "
                f"auto '{e_auto}' != off '{e_off}' - grid over-exclusion bug (#1) is back?"
            )

        # 2) The Na+ must be reported as unparametrized (warning), and metal
        #    coordination scoring must NOT be declared ACTIVE.
        if "coordination scoring is ACTIVE" in out_auto:
            raise AssertionError(
                "coordination scoring reported ACTIVE for an unparametrized metal (Na)"
            )
        if "standard Vina interactions will be used" not in out_auto:
            raise AssertionError("expected unparametrized-metal warning for Na in output:\n"
                                 + out_auto[-2000:])

    print("UNPARAMETRIZED METAL: PASS (auto == off, standard Vina treatment preserved)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
