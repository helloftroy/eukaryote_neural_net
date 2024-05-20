import pandas as pd

# Define start codons
start_codons = "ATG"

# Read the CSV file into a DataFrame
df = pd.read_csv('filtered_length_top_utr.csv', sep='\t', header=None)

# Assuming the sequences are in the fourth column (index 3)
sequences = df.iloc[:, 3]

# Create sliding windows without start codons
window_size = 10
subsequence_length = 19
sliding_windows = []
for sequence in sequences:
    i = 0
    while i < len(sequence) - subsequence_length + 1:
        subsequence = sequence[i:i+subsequence_length]
        if start_codons in subsequence:
            i += 3  # Skip the 3 base pairs
        else:
            sliding_windows.append(subsequence)
            i += window_size  # Move to the next window

# Write the sliding windows without start codons to a new CSV file
with open('nineteen.csv', 'w') as file:
    for window in sliding_windows:
        file.write(window + '\n')

