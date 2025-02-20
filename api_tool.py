import requests
import csv
import pprint
from Bio import Entrez, SeqIO
import os

phagesdb_api_url = "https://phagesdb.org/api/"
ncbi_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
Entrez.email = "jonahchecketts@gmail.com"
Entrez.tool = "phage_host_specificity_analysis"

while True:
    commands = input("Select Command. Type help for list of commands.\n").split(" ")

    if commands[0] == "help" or commands[0] == "h":
        print("quit - exit the program")
        print("cluster <name> <outfile> - get the genetic code of all phages in cluster <name> and output them as a csv to <outfile>.")
        print("cluster_gb <name> - get the genbank file of all phages in cluster <name> and output the files into a folder named cluster_<name>_gb_files")
    elif commands[0] == "cluster":
        response = requests.get(phagesdb_api_url + "clusters/" + commands[1] + "/phagelist/", verify=False)
        results = response.json()['results']
        phage_access = {i["phage_name"] : i['genbank_accession'] for i in results}
        if ".fasta" not in commands[2]:
            commands[2] += ".fasta"
        with open(commands[2], "w", newline="") as f:
            for name, access in phage_access.items():
                if access != "":
                    try:
                        stream = Entrez.efetch(db="nucleotide", id=access, rettype="fasta", retmode="text")
                        phage_data = stream.read()
                        print(phage_data)
                        f.write(phage_data)
                        break
                    except:
                        continue
    elif commands[0] == "cluster_gb":
        response = requests.get(phagesdb_api_url + "clusters/" + commands[1] + "/phagelist/", verify=False)
        results = response.json()['results']
        phage_access = {i["phage_name"] : i['genbank_accession'] for i in results}
        folder_name = "cluster_" + commands[1] + "_gb_files/"
        os.mkdir(folder_name)
        for name, access in phage_access.items():
            with open(folder_name + name + ".gb", "w", newline="") as f:
                if access != "":
                    try:
                        stream = Entrez.efetch(db="nucleotide", id=access, rettype="gb", retmode="text")
                        phage_data = stream.read()
                        f.write(phage_data)
                    except:
                        continue
    elif commands[0] == "other_data":
        response = requests.get(phagesdb_api_url + "clusters/" + commands[1] + "/phagelist/", verify=False)
        results = response.json()['results']
        if ".csv" not in commands[2]:
            commands[2] += ".csv"
        rows = []
        with open(commands[2], "w", newline="") as f:
            col_names = ["phage_name", "genbank_access_id", "GC_percent", "num_tRNAs", "fasta_url"]
            writer = csv.DictWriter(f, fieldnames=col_names)
            writer.writeheader()
            for result in results:
                rows.append({"phage_name": result['phage_name'], "genbank_access_id": result['genbank_accession'], "GC_percent": result['gcpercent'], "num_tRNAs": result['num_tRNAs'], "fasta_url": result['fasta_file']})
            writer.writerows(rows)
    elif commands[0] == "quit" or commands[0] == "q":
        break