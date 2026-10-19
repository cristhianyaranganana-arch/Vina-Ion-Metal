#!/usr/bin/env bash
# Byte-identical regression: receptor with no parametrized metal vs AutoDock Vina 1.2.7.
#
# Usage:
#   tests/test_regression_no_metal.sh /path/to/vina /path/to/vina_original_1.2.7
#
# Both binaries must accept the same CLI. The comparison uses --local_only,
# which exercises scoring and analytic gradients without Monte Carlo randomness.

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <vinametal-binary> <reference-vina-1.2.7-binary>" >&2
    exit 2
fi

VINA_BIN=$1
REFERENCE_BIN=$2
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DATA="$HERE/data"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

run() {
    "$1" --receptor "$DATA/1iep_receptor.pdbqt" \
         --ligand "$DATA/1iep_ligand.pdbqt" \
         --center_x 15.19 --center_y 53.903 --center_z 16.917 \
         --size_x 20 --size_y 20 --size_z 20 \
         --local_only --seed 42 --out "$2"
}

run "$VINA_BIN" "$TMP/vinametal.pdbqt" > "$TMP/vinametal.log" 2>&1
run "$REFERENCE_BIN" "$TMP/reference.pdbqt" > "$TMP/reference.log" 2>&1

echo "VinaMetal : $(grep -m1 'Estimated Free Energy' "$TMP/vinametal.log")"
echo "Reference : $(grep -m1 'Estimated Free Energy' "$TMP/reference.log")"

if cmp -s "$TMP/vinametal.pdbqt" "$TMP/reference.pdbqt"; then
    echo "RESULT: BYTE_IDENTICAL"
else
    echo "RESULT: DIFFERENT (the scoring/gradient path changed without metals)" >&2
    diff -u "$TMP/reference.pdbqt" "$TMP/vinametal.pdbqt" | head -40 >&2 || true
    exit 1
fi
