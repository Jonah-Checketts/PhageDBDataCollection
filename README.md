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
## Re-Annotating the Phage Sequences
Note: This analyis takes a LONG time, and so an example folder with the output is given in the **Coding Sequences Data** link. (For TA's, install that folder and skip the steps below.)
<br>
The Phage coding sequences require a new, reannotation as the NCBI annotations are incorrect. Some of our analyses require these reannotated coding sequences.
<br>
We used Pharokka v1.7.5 to conduct this analysis. PHold should be used in conjunction with Pharokka to ensure the most accurate annotations, but the analysis took too long to do both. Pharokka has a large db it needs to run. **installing_pharokka.bash** needs to be ran to install both Pharokka and its db.
<br>
All scripts for using Pharokka are in the **Phage_ReAnnotations** folder. Steps for running the annotations on a list of given phage clusters:

1. Get the names of each of your clusters (BD, BE, CZ, etc.)
2. Run **submit_pharokka.sh** with your clusters as arguments. Example: **submit_pharokka.sh BD BE P CZ**
3. This will run **run_clusters.bash** which runs **run_pharokka_updated** and **process_gbk_to_CDS** for each cluster. 
4. **run_clusters.bash** will use the **api_tool.py** to download the raw data files for each phage cluster directly from NCBI. 
5. The final output will be a folder with the CDS fasta files for each phage in a given cluster. These can then be used in other analyses (like Codon Usage Bias).

Running this process on one cluster usually takes 4 hours long, so the 10 example clusters used in the analysis are avaliable in the following drive:
[Coding Sequences Data](https://drive.google.com/drive/folders/1YSy4Ht9sVDazFc9X_gCZMCKSwuP4yqQR?usp=sharing)
<br>
The host CDS is the *Streptomyces* folder. This folder is not from the analysis above, but is downloaded directly from NCBI (coding sequences genbank file).

## Selection Pressure Analysis
The selection pressure software MEME requires a fasta file with the coding sequence for the target gene for each phage. It also requires the stop codon to be deleted from the sequence. Once you have the correctly annotated sequences following the instructions above, we need to make a file that meets those requirements. <br>

With the conference deadlines coming up, I manually opened each file, found the minor tail protein sequence, and copied the sequence to a new fasta file. After that file was built I manually deleted the stop codon from the end of each sequence. I ended up with the file called "all_phages_minor_tail.fasta". To save you some time, that file is in the SelectionPressure folder. <br>

That file can now be uploaded to the MEME software. To access the software online, go to datamonkey.org > Sites > Episodic > Meme
Now that MEME has been selected, you can upload the file where prompted in the top left. The other parameters can be left to their default setting. At the bottom click run analysis. Your results should appear in less than 10 minutes.


## Running the Codon Usage Bias Analysis, RSCU Figure - CodonUsageBiasAnalysis folder
This analysis requires the reannotated coding sequences (see above). The example files are in the **Coding Sequences Data** folder from before.
<br>
To run the analysis, use the **Calculating_RSCU_Values.py**. It requires command line arguments. The first argument should be the bacterial CDS fasta file, following by the bacterias GC% content. Then you do the same for each phages.
You can download and put the example Coding Sequences Data in the directory and use each clusters folder path. I.e. when the script runs, put in the folder path /CodonUsageBiasAnalysis/Coding_Sequences_Data/cluster_AK_cds with its GC content.
<br>
For the example files, each one has the following GC contents:
* AK : 61.1
* AS : 66.7
* AU : 50.3
* BD : 66.4
* BE : 49.6
* BI : 59.6
* BU : 54.1
* CZ : 66.4
* DJ : 51.5
* P : 67
<br>
Streptomyces bacteria has a GC content of 71.
<br>
The first script will create a .csv file called **RSCU_Differences_vs_Bacteria.csv**. The plotting script **Plotting_Scatterplot.py** will take this .csv file and plot a heatmap scatterplot of the RSCU values. When you run the script, it will ask you for the GC content of the host bacteria, and a list of which phages infect the host.
<br> For the example data, the clusters "BD, BE, BI" infect the bacterial host Streptomyces.
<br> Note: This generated plot is slightly different than the one on our poster presentation, as the clusters are sorted alphabetically instead of randomly.
<br>The plotting script will generate a .png file called "scatterplot_RSCU_difference.png".
<br>This gives a figure for codon usage bias.


# Reproducing Figure 1 and Figure 2
Go into the file "Figures 1 and 2". The necessary files for each figure are in the folders "figure1" and "figure2"
Make sure file names are consistent with this guide.

Figure 1 (GC content scatterplot)

1. Use get_scatterplot_data.py on the command line "python get_scatterplot_data.py" when prompted for a command, put, "other_data BD phage_host_gc"
2. On the command line, run "python add_host_gc.py" to take the data from "streptomyces_csv_data.csv" and add in the proper host gc data. you will get a new file: "phage_host_gc_full.csv". Note: I used the bacteria_tool.py in the "scripts" folder on github to generate the file "streptomyces_csv_data.csv". This file is already in the figure1 folder and the data folder for your convenience, but you can run bacteria_tool.py to get it yourself if you want. Also, for three of the hosts, I had to manually add in the GC percents from NCBI. It is built into the add_host_gc.py file, but if you want to see where the data comes from here are the links:

Streptomyces azureus: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001270025.1/"
Streptomyces venezuelae: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008639165.1/"
Streptomyces coelicolor: "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_008931305.1/"

3. Now run "#Scatterplot.py". This .py file uses the new files that have been created in steps 1 and 2.
4. Done, you can save the figure as a png if you want.


Figure 2 (trees):

Phage tree:
1. Use api_tool.py that is found in the "scripts" folder. On the command line, run "python api_tool.py" when prompted for a command, put, "cluster_gb BD". This will take a minute to fully download the files.
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
