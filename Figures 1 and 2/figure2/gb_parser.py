import re
import os

# Define gene keywords with possible variations
gene_keywords = {
    "head_to_tail_adaptor": ["head-to-tail adaptor"],
    "head_to_tail_stopper": ["head-to-tail stopper"],
    "minor_tail_protein": ["minor tail protein"],
    "tail_terminator": ["tail terminator"]
}

folder_name = "cluster_BD_gb_files"
fasta_file_name = "streptomyces_phage_proteins.fasta"
output_path = os.path.join(os.getcwd(), fasta_file_name)

# Regex patterns to extract gene and translation fields
product_pattern = r'/product="([^"]+)"'
translation_pattern = r'/translation="([A-Z\n\r\t ]+)"'

# Dictionary to store phage sequences
phage_sequences = {}

for file in sorted(os.listdir(folder_name)):  # Sort files alphabetically
    curr_label = file.split(".")[0]  # Extract file prefix as label
    with open(os.path.join(folder_name, file), "r") as in_f:
        file_text = in_f.read()

    # Find all product names and corresponding translations
    product_matches = re.findall(product_pattern, file_text, flags=re.MULTILINE)
    translation_matches = re.findall(translation_pattern, file_text, flags=re.MULTILINE)

    if not product_matches or not translation_matches:
        continue  # Skip files with no valid sequences

    # Dictionary to store gene sequences
    gene_sequences = {gene: None for gene in sorted(gene_keywords.keys())}

    for product, translation in zip(product_matches, translation_matches):
        product_clean = product.lower().strip()
        for gene_name, keyword_list in gene_keywords.items():
            if any(keyword.lower() in product_clean for keyword in keyword_list):
                protein_sequence = translation.replace("\n", "").replace(" ", "").replace("'", "")
                gene_sequences[gene_name] = protein_sequence
                break  # Stop checking once a match is found for this product

    # Check if all genes are present
    if None in gene_sequences.values():
        missing_genes = [gene for gene, seq in gene_sequences.items() if seq is None]
        print(f"Skipping {curr_label} due to missing genes: {', '.join(missing_genes)}")
        continue  # Exclude this phage from the FASTA file

    # Concatenate genes in alphabetical order
    concatenated_sequence = "".join(gene_sequences[gene] for gene in sorted(gene_sequences))

    # Store phage sequence
    phage_sequences[curr_label] = concatenated_sequence

# Write phages to FASTA file in alphabetical order
with open(output_path, "w") as out_f:
    for phage in sorted(phage_sequences.keys()):  # Sort phages alphabetically
        out_f.write(f">{phage}\n")
        out_f.write(f"{phage_sequences[phage]}\n")

print(f"FASTA file saved to: {output_path}")
