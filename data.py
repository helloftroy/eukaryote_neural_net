
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import torch
import numpy as np

class MFDataset(torch.utils.data.Dataset):

    def __init__(self, classes=10):
        self.max_length = 12000
        self.mfs = pd.read_csv("maturation_eff.csv")
        self.genes = SeqIO.to_dict(SeqIO.parse("genes_of_interest.fasta", "fasta"))
        #files = pd.read_csv('Cum_sum.csv')
        #self.start = dict(zip(files['intron_id'], files['start']))
        #self.length = dict(zip(files['intron_id'], files['length']))
        print("Number of genes before filtering:", len(self.genes))
        filtered_genes = {gene_id: sequence for gene_id, sequence in self.genes.items() if len(sequence) <= self.max_length and len(sequence) >= 15}
        self.mfs = self.mfs[self.mfs['intron_id'].isin(filtered_genes.keys())]
        print("Number of genes after filtering:", len(self.mfs))

        self.mfs['class'], class_boundaries = pd.qcut(self.mfs['sa.eff'], classes, labels=False, retbins=True, duplicates='drop')
    
    # Drop the middle third of genes from the dataset
       # self.mfs = self.mfs[~(self.mfs['class']==1)]
       # self.mfs = self.mfs[~(self.mfs['class']==2)]
       # self.mfs['class'].replace(3, 1, inplace=True)

        #for when you are making 6 classes 
        """class_0_indices = list(self.mfs[self.mfs['class'] == 0].index)
        np.random.shuffle(class_0_indices)
        half = len(class_0_indices) // 2
        self.mfs.loc[class_0_indices[:half], 'class'] = 0
        self.mfs.loc[class_0_indices[half:], 'class'] = classes-1"""
        #for weighted 10 classes
        #class_boundaries = self.mfs['sa.eff'].quantile([0.9, 1.0])  # Calculate the xth percentile
        #self.mfs['class'] = (self.mfs['sa.eff'] > class_boundaries.iloc[0]).astype(int)  # Assign 1 to top 10%, 0 to bottom 90%
        class_counts = self.mfs['class'].value_counts()
        pd.set_option('display.max_rows', None)
        pd.set_option('display.width', None)

        print(self.mfs[['class', 'sa.eff']].head(100))
        print("Class Counts with Boundaries:")
        for class_label, count, boundary in zip(class_counts.index, class_counts, class_boundaries):
            print(f"Class {class_label} - Count: {count}, Boundary: {boundary}")

    def __getitem__(self, idx):
        sample = self.mfs.iloc[idx]
        gene_id = sample["intron_id"]
        eff = sample["class"]
        sa_eff = sample["sa.eff"]
        
        #start_gene_id = gene_id[:7].replace('t', 'g', 1) + gene_id[7:]
        #start = (self.start[start_gene_id]/373245519)*2 - 1
        #length = (self.length[start_gene_id]/self.max_length)*2 - 1
        gene = str(self.genes[gene_id].seq)

        return gene, eff, sa_eff#, start, length

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
