# PhageDBDataCollection
A tool to request phage-host relationship data from PhagesDB and NCBI online databases and analyze that data using several different tools.

## data folder
This stores a lot of the data we used in our analysis of phages and hosts. What specific data is used in what analyses are specified later.

## scripts folder
### api_tool.py
This is the main api tool to get data from phagesDB and NCBI. If you want information about what commands you can run type help. This will print out a list of commands, their usage and what data it returns. This tool was used to get the data for annotating the phage genomes we used, getting GC content of phages for our GC content analysis, and most of our other analyses.

### bacteria_tool.py
This tool looked at all of the files in the data/bacteria_fastas folder and manually calculated their GC content to analyze along with the phage GC content we got using api_tool.py

### gb_parser.py
This tool looks through a folder of genbank files and gets all of the translated proteins to be used for analysis in proteinOrtho. For our analysis we got proteins from all of the files in the bacteria_genbank and Streptomyces_phage_gb folders. The resulting fasta files can be found in the streptomyces_proteins.faa and streptomyces_phage_proteins.faa respectively.

### streptomyces_cluster.py
This generates the GC content cluster plot James and Jonah used on their poster using the strepotmyces_phages.csv file.

## Re-Annotating the Phage Sequences - Phage_ReAnnotations folder
The Phage coding sequences require a new, reannotation as the NCBI annotations are incorrect.
We used Pharokka v1.7.5 to conduct this analysis. PHold should be used in conjunction with Pharokka to ensure the most accurate annotations, but the analysis took too long to do both.
All scripts for using Pharokka are in the **Phage_ReAnnotations** folder. Steps for running the annotations on a list of given phage clusters:
1. Get the names of each of your clusters (BD, BE, CZ, etc.)
2. Run **submit_pharokka.sh** with your clusters as arguments. Example: **submit_pharokka.sh BD BE P CZ**
3. This will run **run_clusters.bash** which runs **run_pharokka_updated** and **process_gbk_to_CDS** for each cluster. 
4. The final output will be a folder with the CDS fasta files for each phage in a given cluster. These can then be used in other analyses (like Codon Usage Bias).

Running this process on one cluster usually takes 4 hours long, so the 10 example clusters used in the analysis are avaliable in the following drive:
[Coding Sequences Data](https://drive.google.com/drive/folders/1YSy4Ht9sVDazFc9X_gCZMCKSwuP4yqQR?usp=sharing)
You will use these sequences for the Codon Usage Bias Analysis below:

## Running the Codon Usage Bias Analysis - CodonUsageBiasAnalysis folder
To run the Codon Usage Bias Analysis and generate a .csv file that contains the RSCU values for each codon in each specified phage cluster. The .csv file will be named **RSCU_Difference_From_Bacteria_Multi_Cluster_df.csv**. A .fasta file with the bacterial CDS of the host is needed and needs to be specified. There is an example Streptomyces.fasta file for this analysis, but it can be replaced by downloading the coding sequences for the selected bacteria.
To run the analysis, use the **Calculating_RSCU_Values.py**. It requires command line arguments. The first argument should be the bacterial CDS fasta file, following by the bacterias GC% content. 