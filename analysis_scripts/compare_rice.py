import pandas as pd
from scipy.stats import f_oneway, ttest_ind
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt



# Read the CSV files into pandas DataFrames
df_leaf = pd.read_csv('rice.maturation_eff_exp_genes.csv')
df_callus = pd.read_csv('batch5_rice.maturation_eff_exp_genes.csv')

t_statistic_eff, t_p_value_eff = ttest_ind(df_leaf['maturation_eff'], df_callus['maturation_eff'])
t_statistic_rna, t_p_value_rna = ttest_ind(df_leaf['spkm.rna'], df_callus['spkm.rna'])

# statistics
print("T-statistic for rice maturation factor:", t_statistic_eff)
print("P-value for rice maturation factor:", t_p_value_eff)
print("T-statistic for rice rna spkm:", t_statistic_rna)
print("P-value for rice rna spkm:", t_p_value_rna)

# P-value
# Add suffixes to the 'maturation_eff' column
df_leaf = df_leaf.rename(columns={'maturation_eff': 'maturation_eff_11'})
df_callus = df_callus.rename(columns={'maturation_eff': 'maturation_eff_5'})

# Merge the DataFrames based on the gene names
merged_df = df_leaf.merge(df_callus, on='tx_id')

merged_dt = pd.DataFrame(merged_df)
efficiency_ranges = {
    'Blue': (-10, -2),
    'Green': (-2, -1),
    'Black': (-1, 1),
    'Yellow': (1, 2),
    'Orange': (2, 4),
    'Red': (4, 10)
}

# Categorize genes based on maturation efficiency ranges for each tissue
leaf_categories = pd.cut(merged_dt['maturation_eff_11'], bins=[-10, -2, -1, 1, 2, 4, float('inf')], labels=efficiency_ranges.keys())
callus_categories = pd.cut(merged_dt['maturation_eff_5'], bins=[-10, -2, -1, 1, 2, 4, float('inf')], labels=efficiency_ranges.keys())

# Count the number of genes in each category for each tissue
leaf_counts = leaf_categories.value_counts().reindex(efficiency_ranges.keys(), fill_value=0)
callus_counts = callus_categories.value_counts().reindex(efficiency_ranges.keys(), fill_value=0)

# Create the stacked bar chart
plt.figure(figsize=(10, 6))
barWidth = 0.5
bars = [leaf_counts, callus_counts]
colors = ['blue', 'green', 'yellow', 'orange', 'red']

# Create bars
for i in range(len(bars)):
    if i == 0:
        plt.bar(efficiency_ranges.keys(), bars[i], color=colors[i], edgecolor='grey', width=barWidth, label='batch11')
    else:
        plt.bar(efficiency_ranges.keys(), bars[i], bottom=bars[i-1], color=colors[i], edgecolor='grey', width=barWidth, label=['batch5'][i-1])

# Add labels
plt.xlabel('Maturation Efficiency Ranges')
plt.ylabel('Number of Genes')
plt.title('Distribution of Genes Across Maturation Efficiency Ranges for two batches')
plt.legend()
#plt.savefig('stacked_bar_chart.png')
plt.show()

# Now merged_df contains the efficiency metrics for each gene across all three conditions
# Calculate the difference in efficiency between the conditions for each gene
merged_df['diff_leaf_callus'] = abs(merged_df['maturation_eff_11'] - merged_df['maturation_eff_5'])

# Find genes with very different efficiency between conditions
threshold = 4 # Define a threshold for what is considered "very different"
different_genes = merged_df[(merged_df['diff_leaf_callus'] > threshold)]

# Overall comparison of efficiency between the three files
mean_efficiency = pd.DataFrame({
    'mean_11_efficiency': [df_leaf['maturation_eff_11'].mean()],
    'mean_5_efficiency': [df_callus['maturation_eff_5'].mean()],
})

print("\nGenes with very different efficiency between conditions:")
print(different_genes)

print("\nOverall comparison of efficiency between the three files:")
print(mean_efficiency)
