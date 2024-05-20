import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import matplotlib.colors as mcolors
from scipy import stats

# Uses : make a hexplot with the input command (1) the name of the dataframe and (2) the title


parser = argparse.ArgumentParser(description='Plot RNA vs. Maturation Score')
parser.add_argument('dataframes', type=str, nargs=2, help='Name of the dataframe')
parser.add_argument('output_file', type=str, help='Name of the output file')
args = parser.parse_args()

# Create a figure and axis for the subplots
fig, axes = plt.subplots(1, 5, figsize=(25, 5))

for i in range(5):    
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(args.dataframes[i])

    # Use Seaborn to create a plot
    joint = sns.jointplot(data=data, x='maturation_eff', y='signal.rna', kind='kde', cmap= 'rocket',ax=axes[i])#, norm=mcolors.LogNorm())

    # Add regression line with standard error and R-squared value
    sns.regplot(x='maturation_eff', y='signal.rna', data=data, scatter=False, ax=joint.ax_joint, color='r')
    slope, intercept, r_value, p_value, std_err = stats.linregress(data['maturation_eff'], data['signal.rna'])
    joint.ax_joint.text(0.8, 0.98, 'R-squared = {:.2f}'.format(r_value**2), ha='center', va='center', transform=joint.ax_joint.transAxes)
    
    title = args.dataframes[i].split('.')[0]    
    axes[i].set_title(title)

# add legend, and adjust layout to prevent cutoff
#plt.colorbar(label='Density')
plt.tight_layout()

plt.savefig(args.output_file + ".png")

"""
# Add a title and labels to the plot
plt.title('RNA vs. Maturation Score')
plt.xlabel('Maturation Efficiency')
plt.ylabel('Signal (RNA)')
plt.savefig("graph_maturation_raw.png") 

"""
