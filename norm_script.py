import sys
import scanpy as sc
import os

# Parse command-line arguments: input .h5ad file and working directory
input_data = sys.argv[1]
working_dir = sys.argv[2]

# Define results and normalization output directories
results_dir = os.path.join(working_dir, "results")
norm_dir = os.path.join(results_dir, "normalized_data")
os.makedirs(norm_dir, exist_ok=True)  # Create the normalization directory if it doesn't exist

# Load the AnnData object (should be pre-filtered and QCed)
adata = sc.read(input_data)

# Step 1: Normalize each cell to have the same total count (library size), target = 10,000
sc.pp.normalize_total(adata, target_sum=1e4)

# Step 2: Apply log1p transformation (natural log(x + 1)) to reduce the effect of outliers
sc.pp.log1p(adata)

# Step 3: Identify highly variable genes (HVGs) based on mean expression and dispersion thresholds
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Step 4: Store the full (log-normalized) data in .raw for future reference
adata.raw = adata

# Step 5: Subset the AnnData object to only include highly variable genes
adata = adata[:, adata.var.highly_variable]

# Step 6: Regress out technical effects: total UMI counts and % mitochondrial gene expression
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Step 7: Scale the data to unit variance and zero mean, clipping values to max of 10
sc.pp.scale(adata, max_value=10)

# Save the processed (normalized, HVG-filtered, scaled) AnnData object
adata.write(os.path.join(norm_dir, "log_normalized.h5ad"))

# Save the list of highly variable genes to CSV
adata.var[adata.var.highly_variable].to_csv(os.path.join(norm_dir, "highly_variable_genes.csv"))
