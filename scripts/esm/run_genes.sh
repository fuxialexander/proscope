#!/bin/bash
set -e 

GENES_FILE="${1}"
MODEL_DIR="${2}"
OUT_DIR="${3}"

for i in 1 2 3 4 5; do
    echo "Running ESM-1v seed v$i..."
    python get_gene.py --genes_file $GENES_FILE --model_location "${MODEL_DIR}/esm1v_t33_650M_UR90S_$i.pt" --out_dir $OUT_DIR
done
