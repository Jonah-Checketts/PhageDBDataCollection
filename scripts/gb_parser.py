import re
import os
import mmap
import pprint

phage_folder_name = "data/Streptomyces_phage_gb"
bac_folder_name = "data/bacteria_genbank"
phage_fasta_file_name = "data/streptomyces_phage_proteins.faa"
bac_fasta_file_name = "data/streptomyces_proteins.faa"
pattern = '/translation="([A-Z\n\r\t ]+)"'

with open(phage_fasta_file_name, "w") as phage_f:
    for file in os.listdir(phage_folder_name):
        curr_label = file[0:file.index(".")]
        file_text = ""
        with open(phage_folder_name + "/" + file, "r") as in_f:
            file_text = in_f.read()
        protein_sequences = re.findall(pattern, file_text, flags=re.MULTILINE)
        i = 0
        for protein in protein_sequences:
            single_line = protein.replace("\n", "").replace(" ", "").replace("'", "")
            phage_f.write(">" + curr_label + " " + str(i) + "\n")
            phage_f.write(single_line + "\n")
            i += 1
with open(bac_fasta_file_name, "w") as bac_f:
    for file in os.listdir(bac_folder_name):
        curr_label = file[0:file.index(".")]
        file_text = ""
        with open(bac_folder_name + "/" + file, "r") as in_f:
            file_text = in_f.read()
        protein_sequences = re.findall(pattern, file_text, flags=re.MULTILINE)
        i = 0
        for protein in protein_sequences:
            single_line = protein.replace("\n", "").replace(" ", "").replace("'", "")
            bac_f.write(">" + curr_label + " " + str(i) + "\n")
            bac_f.write(single_line + "\n")
            i += 1