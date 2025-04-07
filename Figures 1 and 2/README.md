# Reproducing Figure 1 and Figure 2

Figure 1 (GC content scatterplot)

Make sure file names are consistent with this guide:

1. Use get_scatterplot_data.py on the command line "python get_scatterplot_data.py" when prompted for a command, put, "other_data BD phage_host_gc"
2. Run add_host_gc.py to take the data from "streptomyces_csv_data.csv" and add in the proper host gc data. you will get a new file: "phage_host_gc_full.csv". Note: I used the bacteria_tool.py in the "scripts" folder on github to generate the file "streptomyces_csv_data.csv". This file is already in the figure1 folder and the data folder for your convenience, but you can run bacteria_tool.py to get it yourself if you want. For three of the hosts, I had to manually add in the GC percents from NCBI. It is built into the add_host_gc.py file, but if you want to see where the data comes from here are the links:

Streptomyces azureus: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001270025.1/"
Streptomyces venezuelae: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008639165.1/"
Streptomyces coelicolor: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008931305.1/"

3. with the new file obtained in 3 and modified in 4, run "#Scatterplot.py"
4. Done, you can save the figure as a png if you want.


Figure 2 (trees):

Phage tree:
1. Use api_tool.py on the command line "python api_tool.py" when prompted for a command, put, "cluster_gb BD". This will take a minute to fully download the files.
2. With the resultant cluster_BD_gb_files folder in the same directory, run fig2_gb_parser.py file.
3. With the "streptomyces_phage_proteins.fasta" that results from that, in the terminal run "mafft --auto streptomyces_phage_proteins.fasta > concatenated_aligned.fasta" You would need to have MAFFT installed before this, I did it on my mac using "brew install mafft"
4. Install IQtree, I used "brew install iqtree". Then run the following command on the terminal using the "concatenated_aligned.fasta" file made from the MAFFT alignment in the previous step: "iqtree2 -s concatenated_aligned.fasta --alrt 1000 -B 1000 -T 4". Note: "iqtree2" will need to be replaced by the path to where you installed the iqtree2 Unix executable file. For me, it was found within the iqtree file that was installed with the following path "iqtree-2.3.6-macOS-arm/bin/iqtree2"
5. Many files will be made from this. The one we care about is the ".treefile". It should be named "concatenated_aligned.fasta.treefile".  This takes about 15 minutes to run.
6. Upload the treefile into the interactive tree of life website, (https://itol.embl.de/). 

NOTE: There is a degree of randomness associated with iqtree, so the output may not be exactly the same as what is on the poster. This is to be expected. The general trends we discussed in our analysis will be consistent however.

7. The coloring and text annotations were done manually using the ITOL interface.

Bacteria Tree:
1. Download the Streptomyces 16S rRNA fasta files using the NCBI links found in "16StRNA_Links.txt". Save all these files in their own folder, and make sure they all end in ".fasta"
2. Concatenate all the 16S rRNA files by going into their file and running "cat *.fasta > combine_bac.fasta"
3. With the "combine_bac.fasta" that results from that, in the terminal run "mafft --auto combine_bac.fasta > concatenated_aligned.fasta"  You should already have MAFFT installed from the phage tree process.
4. Install IQtree, I used "brew install iqtree". Then run the following command on the terminal using the "concatenated_aligned.fasta" file made from the MAFFT alignment in the previous step: iqtree2 -s concatenated_aligned.fasta --alrt 1000 -B 1000 -T 4. Note: "iqtree2" will need to be replaced by the path to where you installed the iqtree2 Unix executable file. For me, it was found within the iqtree file that was installed with the following path "iqtree-2.3.6-macOS-arm/bin/iqtree2"
5. Many files will be made from this. The one we care about is the ".treefile". It should be named "concatenated_aligned.fasta.treefile". 
6. Upload the treefile into the interactive tree of life website, (https://itol.embl.de/)
7. The coloring and text annotations were done manually using the ITOL interface. 

NOTE: There is a degree of randomness associated with iqtree, so the output may not be exactly the same as what is on the poster. This is to be expected. The general trends we discussed in our analysis will be consistent however.
