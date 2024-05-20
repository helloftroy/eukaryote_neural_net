import pandas as pd
import numpy as np
# Read the GFF file into a pandas DataFrame
df = pd.read_csv('filtered_output.gff', sep='\t', header=None, comment='#')

# Filter the DataFrame to keep only 'gene' entries
gene_df = df[df[2] == 'gene'].drop(2, axis=1)
gene_df.columns = ['chr', 'ens', 'start_gene', 'end_gene', 'score', 'strand', 'dot', 'gene_info']
# Extract gene IDs from the attributes and create a new column
gene_df['gene_id'] = gene_df['gene_info'].str.extract(r'ID=([^;]+)')
gene_df = gene_df.drop('gene_info', axis=1)
# Filter the DataFrame to keep only 'polypeptide' entries
selected_columns = [3, 4, -1]  # Assuming 0-based indexing
polypeptide_df = df[df[2] == 'polypeptide']
polypeptide_df = polypeptide_df.iloc[:, selected_columns]
polypeptide_df.columns = ['start_polypeptide', 'end_polypeptide', 'gene_info']

# Extract gene IDs from the attributes and create a new column
polypeptide_df['gene_id'] = polypeptide_df['gene_info'].str.extract(r'ID=([^;]+)')
polypeptide_df['polypeptide_id'] = polypeptide_df['gene_info'].str.extract(r'Name=([^;]+)')
polypeptide_df = polypeptide_df.drop('gene_info', axis=1)
print(polypeptide_df)
polypeptide_df.reset_index(drop=True, inplace=True)
gene_df.reset_index(drop=True, inplace=True)


gene_df['gene_id'] = gene_df['gene_id'].astype(str)
polypeptide_df['gene_id'] = polypeptide_df['gene_id'].astype(str)
polypeptide_df['gene_id'] = polypeptide_df['gene_id'].str.replace('t', 'g')
polypeptide_df['polypeptide_id'] = polypeptide_df['polypeptide_id'].str.replace('t', 'g')
#Merge the polypeptide and gene DataFrames based on the gene ID
merged_df = pd.merge(gene_df, polypeptide_df, on='gene_id', how='inner')
#print(merged_df)

# Create a dataframe for 5UTR type
utr_5utr_df = merged_df[['chr', 'ens', 'start_gene', 'start_polypeptide', 'score', 'strand', 'dot', 'gene_id', 'polypeptide_id']].copy()
utr_5utr_df['type'] = '5UTR'
utr_5utr_df['combined_id'] = utr_5utr_df['gene_id'] + ';' + utr_5utr_df['polypeptide_id']

# Reorder the columns to position the 'type' column after the 'ens' column
utr_5utr_df = utr_5utr_df[['chr', 'ens', 'type', 'start_gene', 'start_polypeptide', 'score', 'strand', 'dot','combined_id']]

# Create a dataframe for 3UTR type
utr_3utr_df = merged_df[['chr', 'ens', 'end_polypeptide', 'end_gene', 'score', 'strand', 'dot', 'gene_id','polypeptide_id']].copy()
utr_3utr_df['type'] = '3UTR'
utr_3utr_df['combined_id'] = utr_3utr_df['gene_id'] + ';' + utr_3utr_df['polypeptide_id']

# Reorder the columns to position the 'type' column after the 'ens' column
utr_3utr_df = utr_3utr_df[['chr', 'ens', 'type', 'end_polypeptide', 'end_gene', 'score', 'strand', 'dot', 'combined_id']]
utr_3utr_df.columns = utr_5utr_df.columns
#concatenate
concatenated_df = pd.DataFrame(np.concatenate([utr_5utr_df, utr_3utr_df]))#, ignore_index=True)


# Write the UTR DataFrame to a new GFF file
concatenated_df.to_csv('53utr_output.gff', sep='\t', header=False, index=False)
 

print("UTR GFF file has been created: 53utr_output.gff")
