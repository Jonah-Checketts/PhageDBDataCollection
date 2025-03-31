#!/bin/bash
#Do this third (after using rename_and_combine to get desired combined_phages.fasta dataset)
# Or, do it once you have your desired combined .fasta file provided.
#Run this with the "submit_pharokka.sh" script
set -e

# Variables (adjust as needed)
#INPUT_FILE="combined_phages.fasta"
INPUT_FILE="cluster_BU_raw.fasta" #Make sure to modify this everytime
PHAROKKA_OUT_DIR="output_pharokka_BU_cluster" #Make sure to modify this everytime
GENE_PREDICTOR="phanotate"
PHAROKKA_PREFIX="pharokka"
LOCUS_TAG="Default"
FAST=false
META=true
META_HMM=false
FORCE=true

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: File ${INPUT_FILE} does not exist."
  echo "Please check the spelling and that you have uploaded it correctly."
  exit 1
else
  echo "Input file ${INPUT_FILE} exists."
fi

# Validate gene predictor
if [[ "$GENE_PREDICTOR" != "phanotate" && "$GENE_PREDICTOR" != "prodigal" && "$GENE_PREDICTOR" != "prodigal-gv" ]]; then
  echo "Invalid GENE_PREDICTOR. Please choose from: phanotate, prodigal, prodigal-gv."
  exit 1
fi

# Construct the Pharokka command
COMMAND="pharokka.py -d pharokka_db -i ${INPUT_FILE} -t 4 -o ${PHAROKKA_OUT_DIR} -p ${PHAROKKA_PREFIX} -l ${LOCUS_TAG} -g ${GENE_PREDICTOR}"

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
echo "Zipping the output directory so you can download it all in one go."

ZIP_FILENAME="${PHAROKKA_OUT_DIR}.zip"
# Zip the output directory (the -r flag recurses into directories)
zip -r "${ZIP_FILENAME}" "${PHAROKKA_OUT_DIR}"

echo "Output directory has been zipped to ${ZIP_FILENAME}."
