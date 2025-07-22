import sys
import scanpy as sc
import os

# Parse input arguments: input .h5ad file and working directory
input_data = sys.argv[1]
working_dir = sys.argv[2]

# Define results and dimensionality reduction output directory
results_dir = os.path.join(working_dir, "results")
dr_dir = os.path.join(results_dir, "dimensionality_reduction")
os.makedirs(dr_dir, exist_ok=True)  # Create directory if it doesn't exist

# Load the AnnData object
adata = sc.read(input_data)

# Step 1: Perform PCA (Principal Component Analysis)
# Uses the ARPACK solver, efficient for large sparse datasets
sc.tl.pca(adata, svd_solver='arpack')

# Step 2: Compute the neighborhood graph based on PCA
# This is required before running UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

# Step 3: Run UMAP (Uniform Manifold Approximation and Projection)
# This embeds the data into 2D space for visualization
sc.tl.umap(adata)

# Save the AnnData object with PCA and UMAP results
adata.write(os.path.join(dr_dir, "adata_umap.h5ad"))

# Step 4: Plot and save PCA variance ratio (scree plot)
sc.pl.pca_variance_ratio(adata, log=True, save="_pca_variance_ratio.png")

# Step 5: Plot and save the UMAP embedding
sc.pl.umap(adata, save="_umap_projection.png")

# Step 6: Move generated plots from Scanpy's default 'figures/' folder to 'dimensionality_reduction/'
for figs in os.listdir("figures"):
    os.rename(os.path.join("figures", figs), os.path.join(dr_dir, figs))