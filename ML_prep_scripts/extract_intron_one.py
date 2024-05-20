from Bio import SeqIO

# Input and output file paths
input_file = "output.fasta"
output_file = "intron_one.fasta"

# Filter sequences with 'intron1' in the header
with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        if ".intron1" in record.id and record.id.index(".intron1") == len(record.id) - len(".intron1"):
            new_id = record.id.replace(".intron1", "")  # Remove '.intron1' from the ID
            record.id = new_id
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")



