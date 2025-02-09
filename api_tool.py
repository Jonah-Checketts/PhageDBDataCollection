import requests
import csv
from Bio import Entrez, SeqIO

phagesdb_api_url = "https://phagesdb.org/api/"
ncbi_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
Entrez.email = "jonahchecketts@gmail.com"
Entrez.tool = "phage_host_specificity_analysis"

while True:
    commands = input("Select Command. Type help for list of commands.\n").split(" ")

    if commands[0] == "help":
        print("quit - exit the program")
        print("cluster <name> <outfile> - get the genetic code of all phages in cluster <name> and output them as a csv to <outfile>")
    elif commands[0] == "cluster":
        response = requests.get(phagesdb_api_url + "clusters/" + commands[1] + "/phagelist/", verify=False)
        results = response.json()['results']
        phage_access = {i["phage_name"] : i['genbank_accession'] for i in results}
        if ".csv" not in commands[2]:
            commands[2] += ".csv"
        with open(commands[2], "w", newline="") as f:
            w = csv.DictWriter(f, ["phage_name", "dna_sequence"])
            w.writeheader()
            for name, access in phage_access.items():
                if access != "":
                    try:
                        stream = Entrez.efetch(db="nucleotide", id=access, rettype="fasta", retmode="text")
                        phage_data = stream.read()
                        values = {"phage_name" : name , "dna_sequence": phage_data}
                        w.writerow(values)
                    except:
                        continue
    elif commands[0] == "quit":
        break