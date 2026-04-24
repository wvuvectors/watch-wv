#!/bin/bash
set -euo pipefail

cd $SCRATCH/watch-wv/seqr/

echo "[INFO] Starting parse_demix job on $(hostname)"

# Activate conda
source /shared/software/conda/etc/profile.d/conda.sh
conda activate seqr_freyja

# ---- Check argument ----
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <DEMIX_DIR>"
    exit 1
fi

DEMIX_DIR="$1"

# ---- Validate directory ----
if [[ ! -d "$DEMIX_DIR" ]]; then
    echo "[ERROR] DEMIX directory does not exist: $DEMIX_DIR"
    exit 1
fi

barcodes=()

# Loop through demix files
for file in "$DEMIX_DIR"/demix_*_barcode*; do
    [[ -f "$file" ]] || continue

    # Extract last two digits after "barcode"
    barcode=$(basename "$file" | grep -oP 'barcode\K[0-9]{2}')
    barcodes+=("$barcode")
done

echo "Found barcodes: ${barcodes[*]}"

# Pass barcodes to python script
python parse_demix.py "$DEMIX_DIR" "${barcodes[@]}"
