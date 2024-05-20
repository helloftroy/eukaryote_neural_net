import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the CSV files into pandas dataframes and set 'tx_id' as the index
df1 = pd.read_csv('batch5_rice.maturation_eff_exp_genes.csv').set_index('tx_id')
df2 = pd.read_csv('rice.maturation_eff_exp_genes.csv').set_index('tx_id')

# Read the genes_list.txt file
with open('genes_list.txt', 'r') as file:
    gene_list = [line.strip() for line in file]

# Perform an inner join on the index (tx_id) to ensure overlap
merged_df = df1.join(df2, how='inner', lsuffix='_df1', rsuffix='_df2')

# Create a scatter plot for all genes
plt.figure(figsize=(8, 9))  # Set the figure size
plt.subplot(211)
sns.scatterplot(data=merged_df, x='maturation_eff_df1', y='maturation_eff_df2', hue=merged_df.index.isin(gene_list))
plt.xlabel('maturation_eff for batch 5 rice')
plt.ylabel('maturation_eff for batch 11 rice')
plt.title('Maturation eff for All Genes')

# Filter the merged dataframe based on the gene list
filtered_df = merged_df[merged_df.index.isin(gene_list)]

# Create a scatter plot for genes in gene_list
plt.subplot(212)
sns.scatterplot(data=filtered_df, x='maturation_eff_df1', y='maturation_eff_df2')
plt.xlabel('maturation_eff for batch 5 rice')
plt.ylabel('maturation_eff for batch 11 rice')
plt.title('Maturation_eff for Top Genes')

# Set the spacing between subplots
plt.tight_layout(pad=3.0)

# Display both plots
plt.show()
