#!/bin/bash
# run_clusters.bash
# Usage: run_clusters.bash <cluster1> [<cluster2> ...]

if [ $# -lt 1 ]; then
  echo "Usage: $0 <cluster1> [<cluster2> ...]"
  exit 1
fi

# Loop over every cluster argument
for cluster in "$@"; do
  echo "===================================="
  echo "Processing cluster: $cluster"
  echo "===================================="
  
  # 1) Fetch the FASTA for this cluster
  #    "api_tool.py cluster <clusterName> <outfile>"
  #    e.g. BD -> BD_raw.fasta
  python ~/api_tool.py cluster "$cluster" "${cluster}_raw.fasta"
  
  # 2) Run pharokka on the newly created FASTA
  ./run_pharokka_updated.bash "${cluster}_raw.fasta" "output_pharokka_${cluster}_cluster"
  
  echo "Done processing cluster: $cluster"
  echo
done
