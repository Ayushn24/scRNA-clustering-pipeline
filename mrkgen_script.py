import sys
import scanpy as sc
import pandas as pd
import os

# Get command-line arguments: input .h5ad file and working directory
input_data = sys.argv[1]
working_dir = sys.argv[2]

# Define paths for results and marker gene output
results_dir = os.path.join(working_dir, "results")
marker_dir = os.path.join(results_dir, "marker_genes")
os.makedirs(marker_dir, exist_ok=True)  # Create marker directory if it doesn't exist

# Load the AnnData object from the provided file
adata = sc.read(input_data)

# Perform differential expression analysis using the Wilcoxon test
# This compares gene expression between clusters defined by "leiden"
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# Plot top 25 marker genes per cluster and save figure
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save="_marker_genes.png")

# Move the generated plot(s) from Scanpy's default 'figures/' folder to 'marker_genes/'
for figs in os.listdir("figures"):
    os.rename(os.path.join("figures", figs), os.path.join(marker_dir, figs))

# Extract ranked gene information (names and p-values) from the results
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names  # Get cluster names

# Build a dataframe of marker genes and p-values for all clusters
ranked_genes = pd.DataFrame(
    {f"{group}_{key}": result[key][group]
     for group in groups for key in ['names', 'pvals']}
)

# Save the ranked marker genes as a CSV file
ranked_genes.to_csv(os.path.join(marker_dir, "ranked_genes.csv"))