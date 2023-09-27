#!/bin/bash
set -e 

GENES_FILE="${1}"
UNIPROT_FILE="${2}"
MODEL_DIR="${3}"
OUTPUT_DIR="${4}"

for i in 1 2 3 4 5; do
    echo "Running ESM-1v with seed $i..."
    python ./esm/run_esm1v.py --genes_file $GENES_FILE --uniprot_path $UNIPROT_FILE --model_location "${MODEL_DIR}/esm1v_t33_650M_UR90S_$i.pt" --output_dir $OUTPUT_DIR
done
