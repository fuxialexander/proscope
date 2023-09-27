# Scripts for running dependent packages

## ESM

Zero-shot variant effect predictions from ESM-1v ([Meier et al. 2021](https://www.biorxiv.org/content/10.1101/2021.07.09.450648v2)) are used to compute the ES score. Meta AI has released five comparable models for this task developed with different random seeds. First download the models and UniProt human reference proteome locally by running:

`sh ./esm/download_models.sh`

Once the models and mapping are downloaded, run predictions for each model with the following bash script and command-line arguments:

`sh ./esm/run_esm_preds.sh [genes_file] [mapping_file] [model_dir] [output_dir]`

1. `[genes_file]` points to a file with a list of gene names (UniProtKB accession numbers) separated by newlines.
2. `[mapping_file]` points to the path of the reference proteome downloaded in the previous step. 
3. `[model_dir]` points to the directory containing the five ESM models downloaded in the previous step.
4. Specify an `[output_dir]` where the results will be saved.

This script runs [run_esm1v.py](./esm/run_esm1v.py) for each model, which executes a saturated mutagenesis on each amino acid sequence corresponding to the provided genes. The provided ESM-1v model scores each generated sequence. The output is a feather file containing the scores for each possible mutation for all genes. 

## ColabFold

TODO

## MD

TODO
