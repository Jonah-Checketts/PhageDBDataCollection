# PhageDBDataCollection
A tool to request phage-host relationship data from PhagesDB and NCBI online databases and analyze that data using several different tools.

## Re-Annotating the Phage Sequences
The Phage coding sequences require a new, reannotation as the NCBI annotations are incorrect.
We used Pharokka v1.7.5 to conduct this analysis. PHold should be used in conjunction with Pharokka to ensure the most accurate annotations, but the analysis took too long to do both.
All scripts for using Pharokka are in the **Phage_ReAnnotations** folder. Steps for running the annotations on a list of given phage clusters:
1. Get the names of each of your clusters (BD, BE, CZ, etc.)
2. Run **submit_pharokka.sh** with your clusters as arguments. Example: **submit_pharokka.sh BD BE P CZ**
3. This will run **run_clusters.bash** which runs **run_pharokka_updated** and **process_gbk_to_CDS** for each cluster. 
4. The final output will be a folder with the CDS fasta files for each phage in a given cluster. These can then be used in other analyses (like Codon Usage Bias).

## Running the Codon Usage Bias Analysis

