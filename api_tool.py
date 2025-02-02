import requests
import csv

api_url = "https://phagesdb.org/api/"

while True:
    commands = input("Select Command. Type help for list of commands.\n").split(" ")

    if commands[0] == "help":
        print("quit - exit the program")
        print("cluster <name> <outfile> - get the genetic code of all phages in cluster <name> and output them as a csv to <outfile>")
    elif commands[0] == "cluster":
        response = requests.get(api_url + "clusters/" + commands[1] + "/phagelist/", verify=False)
        results = response.json()['results']
        phage_names = [i['phage_name'] for i in results]
        phages_dna = {}
        for name in phage_names:
            response2 = requests.get(api_url + "genesbyphage/" + name + "/", verify=False)
            results2 = response2.json()['results']
            for result in results2:
                phages_dna[result['GeneID']] = {'phage_name':name, 'gene_name': result['GeneID'], 'start': result['Start'], 'stop': result['Stop'], "protein": result['Translation']}
        if ".csv" not in commands[2]:
            commands[2] += ".csv"
        with open(commands[2], "w", newline="") as f:
            first = True
            for key, values in phages_dna.items():
                if (first):
                    w = csv.DictWriter(f, values.keys())
                    w.writeheader()
                    first = False
                w.writerow(values)             
    elif commands[0] == "quit":
        break

    
    