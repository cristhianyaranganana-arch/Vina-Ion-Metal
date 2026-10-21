#!/usr/bin/env bash
# Byte-identical regression: receptor with no parametrized metal vs AutoDock Vina 1.2.7.
#
# Vina-Ion-Metal (audit #13): this test is registered in CTest. It locates the
# reference Vina 1.2.7 binary in this order:
#   1. $VINA_REFERENCE_BIN environment variable (used by the CI, which downloads
#      the official release binary);
#   2. --reference argument (backwards compatible explicit path);
#   3. otherwise the test SKIPS with exit 0 and an explanatory message.
#
# Both binaries must accept the same CLI. The comparison uses --local_only,
# which exercises scoring and analytic gradients without Monte Carlo randomness.

set -euo pipefail

usage() {
    echo "usage: $0 <vina-ion-metal-binary> [reference-vina-1.2.7-binary]" >&2
    echo "       (or export VINA_REFERENCE_BIN=/path/to/vina_1.2.7)" >&2
}

VINA_BIN=""
REFERENCE_BIN=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --reference)
            REFERENCE_BIN="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            if [ -z "$VINA_BIN" ]; then
                VINA_BIN="$1"
            elif [ -z "$REFERENCE_BIN" ]; then
                REFERENCE_BIN="$1"
            else
                usage
                exit 2
            fi
            shift
            ;;
    esac
done

if [ -z "$VINA_BIN" ]; then
    usage
    exit 2
fi

# Locate the reference binary (audit #13)
if [ -z "$REFERENCE_BIN" ] && [ -n "${VINA_REFERENCE_BIN:-}" ]; then
    REFERENCE_BIN="$VINA_REFERENCE_BIN"
fi
if [ -z "$REFERENCE_BIN" ]; then
    echo "SKIP (regression_no_metal): no reference AutoDock Vina 1.2.7 binary available."
    echo "      To enable this test: export VINA_REFERENCE_BIN=/path/to/vina_1.2.7"
    echo "      (the CI downloads it from the official ccsb-scripps release assets)."
    exit 0
fi
if [ ! -x "$REFERENCE_BIN" ]; then
    echo "FAIL (regression_no_metal): reference binary '$REFERENCE_BIN' is not executable." >&2
    exit 1
fi

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

run "$VINA_BIN" "$TMP/vina_ion_metal.pdbqt" > "$TMP/vim.log" 2>&1
run "$REFERENCE_BIN" "$TMP/reference.pdbqt" > "$TMP/reference.log" 2>&1

echo "Vina-Ion-Metal : $(grep -m1 'Estimated Free Energy' "$TMP/vim.log")"
echo "Reference Vina : $(grep -m1 'Estimated Free Energy' "$TMP/reference.log")"

if cmp -s "$TMP/vina_ion_metal.pdbqt" "$TMP/reference.pdbqt"; then
    echo "RESULT: BYTE_IDENTICAL"
else
    echo "RESULT: DIFFERENT (the scoring/gradient path changed without parametrized metals)" >&2
    diff -u "$TMP/reference.pdbqt" "$TMP/vina_ion_metal.pdbqt" | head -40 >&2 || true
    exit 1
fi
