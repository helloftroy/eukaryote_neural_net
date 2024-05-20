
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import torch

class MFDataset(torch.utils.data.Dataset):

    def __init__(self, classes=10):
        self.mfs = pd.read_csv("maturation_eff.csv")
        self.genes = SeqIO.to_dict(SeqIO.parse("genes_of_interest.fasta", "fasta"))
        print("Number of genes before filtering:", len(self.genes))
        filtered_genes = {gene_id: sequence for gene_id, sequence in self.genes.items() if len(sequence) <= 8000 and len(sequence) >= 15}
        self.mfs = self.mfs[self.mfs['intron_id'].isin(filtered_genes.keys())]
        print("Number of genes after filtering:", len(self.mfs))

        self.mfs['class'], class_boundaries = pd.qcut(self.mfs['sa.eff'], classes, labels=False, retbins=True, duplicates='drop')

        print("Class Counts:", self.mfs['class'].value_counts())
        print("Class Boundaries:", class_boundaries)

    def __getitem__(self, idx):
        sample = self.mfs.iloc[idx]
        gene_id = sample["intron_id"]
        eff = sample["class"]
        sa_eff = sample["sa.eff"]
        gene = str(self.genes[gene_id].seq)

        return gene, eff, sa_eff

    def __len__(self):
        return len(self.mfs)

if __name__ == "__main__":
    # Create an instance of MFDataset
    dataset = MFDataset()

    # Provide a key and return/print the fasta sequence
    key = "Chr_5:439311-439434(+)"
    if key in dataset.genes:
        sequence = str(dataset.genes[key].seq)
        print(sequence)
    else:
        print(f"Key '{key}' not found in the dataset.")
