import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV files into pandas dataframes and set 'tx_id' as the index
df1 = pd.read_csv('batch5_rice.maturation_eff_exp_genes.csv').set_index('tx_id')
df2 = pd.read_csv('rice.maturation_eff_exp_genes.csv').set_index('tx_id')

# Read the genes_list.txt file
with open('genes_list.txt', 'r') as file:
    gene_list = [line.strip() for line in file]

# Perform an inner join on the index (tx_id) to ensure overlap
merged_df = df1.join(df2, how='inner', lsuffix='_df1', rsuffix='_df2')

# Create a scatter plot for all genes
plt.figure(1, figsize=(6, 9))

plt.subplot(211)
scatter_all = plt.scatter(merged_df['maturation_eff_df1'], merged_df['maturation_eff_df2'])
plt.xlabel('maturation_eff for batch5_rice')
plt.ylabel('maturation_eff for rice')
plt.title('Scatter Plot of maturation_eff for All Genes')

# Color the points based on the gene list
colors_all = ['red' if tx_id in gene_list else 'blue' for tx_id in merged_df.index]
scatter_all.set_facecolor(colors_all)

# Filter the merged dataframe based on the gene list
filtered_df = merged_df[merged_df.index.isin(gene_list)]

# Create a scatter plot for genes in gene_list
plt.subplot(212)
scatter_filtered = plt.scatter(filtered_df['maturation_eff_df1'], filtered_df['maturation_eff_df2'])
plt.xlabel('maturation_eff for batch5_rice')
plt.ylabel('maturation_eff for rice')
plt.title('Scatter Plot of maturation_eff for Genes in gene_list')

# Display both plots
plt.tight_layout()
plt.show()
