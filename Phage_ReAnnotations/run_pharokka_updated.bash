#!/bin/bash
# run_pharokka.bash
# Usage: run_pharokka.bash <input_fasta> <output_dir>

set -e

# Grab command-line args
INPUT_FILE="$1"
PHAROKKA_OUT_DIR="$2"

# specify gene predictor and other params
GENE_PREDICTOR="phanotate"
PHAROKKA_PREFIX="pharokka"
LOCUS_TAG="Default"
FAST=false
META=true
META_HMM=false
FORCE=true

# Basic error-check
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: File ${INPUT_FILE} does not exist."
  exit 1
fi

# Construct the Pharokka command
COMMAND="pharokka.py -d pharokka_db -i ${INPUT_FILE} -t 4 -o ${PHAROKKA_OUT_DIR} \
                     -p ${PHAROKKA_PREFIX} -l ${LOCUS_TAG} -g ${GENE_PREDICTOR}"

if [ "$FORCE" = true ]; then
  COMMAND="${COMMAND} -f"
fi

if [ "$FAST" = true ]; then
  COMMAND="${COMMAND} --fast"
fi

if [ "$META" = true ]; then
  COMMAND="${COMMAND} -m"
fi

if [ "$META_HMM" = true ]; then
  COMMAND="${COMMAND} --meta_hmm"
fi

echo "Running pharokka command:"
echo "${COMMAND}"
eval ${COMMAND}

echo "pharokka completed successfully."
echo "Your output is in ${PHAROKKA_OUT_DIR}."
echo "Zipping the output directory so you can download it all in one go..."

ZIP_FILENAME="${PHAROKKA_OUT_DIR}.zip"
zip -r "${ZIP_FILENAME}" "${PHAROKKA_OUT_DIR}"

echo "Output directory has been zipped to ${ZIP_FILENAME}."

# Now call process_gbk_to_CDS.py on the pharokka .gbk
python process_gbk_to_CDS.py "${PHAROKKA_OUT_DIR}/pharokka.gbk"
