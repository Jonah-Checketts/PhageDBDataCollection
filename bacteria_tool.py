import csv
import os

dir_name = "C:/Users/schec/Documents/GitHub/PhageDBDataCollection/bacteria_fastas"

with open("streptomyces_csv_data.csv", "w", newline="") as out_f:
    fieldnames = ['strain_name', 'gc_percent']
    writer = csv.DictWriter(out_f, fieldnames=fieldnames)
    writer.writeheader()
    csv_data = []
    for subdir, dirs, files in os.walk(dir_name):
        for file in files:
            with open(os.path.join(subdir, file), "r") as in_f:
                total_count = 0
                gc_count = 0
                first_line = True
                bacteria_name = ""
                for line in in_f:
                    if first_line:
                        strain_info = line.split(" ")
                        bacteria_name = strain_info[0][1:] + " " + strain_info[1] + " " + strain_info[2]
                        first_line = False
                    else:
                        line = line.strip().upper()
                        total_count += len(line)
                        gc_count += line.count("G")
                        gc_count += line.count("C")
                percent = gc_count / total_count
                csv_data.append({'strain_name': bacteria_name, 'gc_percent': percent})
    writer.writerows(csv_data)



                