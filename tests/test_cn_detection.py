#!/usr/bin/env python3
"""End-to-end check of the automatic receptor-occupancy CN adjustment.

Runs the actual VinaMetal binary with --score_only on the bundled examples and
checks the logged CN target line produced by metal_scoring::report_cn_targets().

Expected lines (recorded during validation):
  examples/glft2 : Mn2+ CN target = ideal 6 - receptor-bound 4 = 2
  examples/1s63  : Zn2+ CN target = ideal 4 - receptor-bound 3 = 1

Usage:
  tests/test_cn_detection.py <vina-binary> <repo-root>
"""

import os
import re
import subprocess
import sys

GRID = {
    "glft2": ("142.54", "32.19", "48.55", "22", "22", "22"),
    "1s63": ("18", "134", "0", "20", "20", "22"),
}

EXPECTED = {
    "glft2": r"Mn2\+ CN target = ideal 6 - receptor-bound 4 = 2",
    "1s63": r"Zn2\+ CN target = ideal 4 - receptor-bound 3 = 1",
}


def run_example(vina, root, name):
    example = os.path.join(root, "examples", name)
    cx, cy, cz, sx, sy, sz = GRID[name]
    cmd = [
        vina,
        "--receptor", os.path.join(example, "receptor.pdbqt"),
        "--ligand", os.path.join(example, "ligand.pdbqt"),
        "--center_x", cx, "--center_y", cy, "--center_z", cz,
        "--size_x", sx, "--size_y", sy, "--size_z", sz,
        "--metals", "auto",
        "--score_only",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise AssertionError(
            f"{name}: vina exited with {result.returncode}\n{output[-2000:]}"
        )
    match = re.search(r"VinaMetal: .*CN target.*", output)
    if not match:
        raise AssertionError(f"{name}: no CN target line in output\n{output[-2000:]}")
    line = match.group(0)
    print(f"{name}: {line}")
    if not re.search(EXPECTED[name], line):
        raise AssertionError(f"{name}: expected /{EXPECTED[name]}/, got '{line}'")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    vina, root = sys.argv[1], sys.argv[2]
    for name in GRID:
        run_example(vina, root, name)
    print("CN DETECTION: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
