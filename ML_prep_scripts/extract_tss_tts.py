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
#print(gene_df)
gene_df.reset_index(drop=True, inplace=True)
gene_df['gene_id'] = gene_df['gene_id'].astype(str)

# Create a dataframe for 5UTR type
utr_5utr_df = gene_df.copy()
utr_5utr_df['type'] = 'TSS'
# Conditional assignment based on the 'strand' column
utr_5utr_df['tss_up'] = utr_5utr_df.apply(lambda row: (row['end_gene'] - 500) if row['strand'] == '-' else row['start_gene'] - 1000, axis=1)
utr_5utr_df['tss_down'] = utr_5utr_df.apply(lambda row: (row['end_gene'] + 1000) if row['strand'] == '-' else row['start_gene'] + 500, axis=1)
# Reorder the columns to position the 'type' column after the 'ens' column
utr_5utr_df = utr_5utr_df[['chr', 'ens', 'type', 'tss_up','tss_down', 'score', 'strand', 'dot','gene_id']]
# Create a dataframe for 3UTR type
"""
utr_3utr_df = gene_df.copy()#merged_df[['chr', 'ens', 'end_polypeptide', 'end_gene', 'score', 'strand', 'dot', 'gene_id','polypeptide_id']].copy()
utr_3utr_df['type'] = 'TTS'
utr_3utr_df['tss_up'] = utr_3utr_df['end_gene'] - 500
utr_3utr_df['tss_down'] = utr_3utr_df['end_gene'] + 500
# Reorder the columns to position the 'type' column after the 'ens' column
utr_3utr_df = utr_3utr_df[['chr', 'ens', 'type', 'tss_up','tss_down', 'score', 'strand', 'dot','gene_id']]

#concatenate
concatenated_df = pd.DataFrame(np.concatenate([utr_5utr_df, utr_3utr_df]))#, ignore_index=True)


# Write the UTR DataFrame to a new GFF file
concatenated_df.to_csv('tss_tts_output.gff', sep='\t', header=False, index=False)
"""
utr_5utr_df.to_csv('tss_tts_output.gff', sep='\t', header=False, index=False)
utr_5utr_df[['chr', 'tss_up', 'tss_down', 'gene_id', 'score', 'strand']].to_csv('utr_output.bed', sep='\t', header=False, index=False)
print("tss_tts.GFF file has been created")
