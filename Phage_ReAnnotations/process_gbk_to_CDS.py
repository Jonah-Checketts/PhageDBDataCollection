#!/usr/bin/env python3
import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def split_cds_by_phage(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for record in SeqIO.parse(input_file, "genbank"):
        phage_id   = record.id
        fullrecord = record.description.split(" ")
        phage_name = fullrecord[3].rstrip(":") 
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
        print("Usage: python process_gbk_to_CDS.py <input_file.gbk>")
        sys.exit(1)

    input_file = sys.argv[1]

    #
    # input_file might be:
    #   ./output_pharokka_P_cluster/pharokka.gbk
    # and the directory name "output_pharokka_P_cluster" includes your cluster name.
    dir_name = os.path.basename(os.path.dirname(input_file))
    prefix = "output_pharokka_"
    suffix = "_cluster"

    # Safely extract the cluster name
    if dir_name.startswith(prefix) and dir_name.endswith(suffix):
        cluster_name = dir_name[len(prefix):-len(suffix)]
    else:
        # fallback if the pattern doesn't match what we expect?
        cluster_name = "unknown"

    output_dir = f"./Coding_Sequence_Data/cluster_{cluster_name}_cds"

    split_cds_by_phage(input_file, output_dir)