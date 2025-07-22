import sys
import scanpy as sc
import os

# Parse input arguments: path to input .h5ad file and working directory
input_data = sys.argv[1]
working_dir = sys.argv[2]

# Define paths for results and clustering output
results_dir = os.path.join(working_dir, "results")
clus_dir = os.path.join(results_dir, "clustering")
os.makedirs(clus_dir, exist_ok=True)  # Create clustering directory if it doesn't exist

# Load the AnnData object from the previous step (after UMAP)
adata = sc.read(input_data)

# Step 1: Perform Leiden clustering at a given resolution
sc.tl.leiden(adata, resolution=0.25)

# Step 2: Plot UMAP colored by Leiden cluster assignment
sc.pl.umap(adata, color=["leiden"], save="_umap_leiden.png")

# Step 3: Optionally perform Louvain clustering for comparison
sc.tl.louvain(adata)

# Step 4: Plot UMAP colored by Louvain cluster assignment
sc.pl.umap(adata, color=["louvain"], save="_umap_louvain.png")

# Step 5: Move generated figures from Scanpy's default 'figures/' folder to clustering directory
for figs in os.listdir("figures"):
    os.rename(os.path.join("figures", figs), os.path.join(clus_dir, figs))

# Step 6: Save the clustered AnnData object
adata.write(os.path.join(clus_dir, "adata_clustered.h5ad"))