#!/usr/bin/env python3
# run with python process_gbk_to_CDS.py ./output_pharokka/pharokka.gbk
import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def split_cds_by_phage(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for record in SeqIO.parse(input_file, "genbank"):
        phage_id   = record.id
        fullrecord = record.description.split(" ")
        phage_name = fullrecord[3]
        phage_name = phage_name[0:-1]
        print(f"Writing {phage_name}:")
        cds_records = []
        for feat in record.features:
            if feat.type == "CDS":
                seq     = feat.extract(record.seq)
                tag     = feat.qualifiers.get("locus_tag", feat.qualifiers.get("gene", [f"{phage_id}_CDS"]))[0]
                product = feat.qualifiers.get("product", [""])[0]
                function= feat.qualifiers.get("function", [""])[0]
 
                header_desc = f"phage={phage_name}; product={product}; function={function}"
                cds_records.append(SeqRecord(seq, id=f"{tag}|{phage_name}", description=header_desc))

        if cds_records:
            out_file = os.path.join(output_dir, f"{phage_name}_{phage_id}.fasta")
            SeqIO.write(cds_records, out_file, "fasta")
            print(f"Wrote {len(cds_records)} CDS → {out_file}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python split_cds.py input_file.gbk")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = "./cluster_" + input_file[18:20] + "_cds"
    split_cds_by_phage(input_file, output_dir)
