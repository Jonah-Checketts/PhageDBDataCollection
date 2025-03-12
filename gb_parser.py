import re
import os
import mmap
import pprint

folder_name = "Streptomyces_Phage_gb"
fasta_file_name = "streptomyces_phage_proteins.fasta"
pattern = '/translation="([A-Z\n\r\t ]+)"'

with open(fasta_file_name, "w") as out_f:
    for file in os.listdir(folder_name):
        curr_label = file[0:file.index(".")]
        file_text = ""
        with open(folder_name + "/" + file, "r") as in_f:
            file_text = in_f.read()
        protein_sequences = re.findall(pattern, file_text, flags=re.MULTILINE)
        i = 0
        for protein in protein_sequences:
            single_line = protein.replace("\n", "").replace(" ", "").replace("'", "")
            out_f.write(">" + curr_label + " " + str(i) + "\n")
            out_f.write(single_line + "\n")
            i += 1