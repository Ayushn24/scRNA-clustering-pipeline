import sys
import zipfile
import os
import scanpy as sc
import numpy as np

# Arguments: zip_path working_dir
zip_path = sys.argv[1]
working_dir = sys.argv[2]
results_dir = os.path.join(working_dir, "results")

# Create result subdirectories
qc_dir = os.path.join(results_dir, "qc_report")
os.makedirs(qc_dir, exist_ok=True)

# Step 1: Unzip data
with zipfile.ZipFile(zip_path, "r") as zip_ref:
    zip_ref.extractall("data")

# Step 2: Read raw data
adata = sc.read_10x_mtx("data", var_names="gene_symbols", cache=True)

# Save raw h5ad
raw_path = os.path.join(qc_dir, "raw_data.h5ad")
adata.write(raw_path)

# Step 3: Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Step 4: Annotate QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

# Step 5: Filter based on gene count distribution
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.02)
adata = adata[(adata.obs.n_genes_by_counts >= lower_lim) & (adata.obs.n_genes_by_counts <= upper_lim)]
adata = adata[adata.obs.pct_counts_mt < 10]  # Filter out cells with >10% mitochondrial genes

# Step 6: Save filtered matrix
filtered_path = os.path.join(qc_dir, "filtered_matrix.h5ad")
adata.write(filtered_path)

# Step 7: QC Plots
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True,
             save="_violin_qc.png")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_scatter_qc.png")

# Move saved figures to qc_dir (Scanpy saves in `.` with `save=...`)
if os.path.exists("figures"):
    for f in os.listdir("figures"):
        os.rename(os.path.join("figures", f), os.path.join(qc_dir, f))

# Step 8: Export gene type annotations
adata.var[adata.var.mt].to_csv(os.path.join(qc_dir, "mitochondrial_genes.csv"))
adata.var[adata.var.ribo].to_csv(os.path.join(qc_dir, "ribosomal_genes.csv"))
adata.var[adata.var.hb].to_csv(os.path.join(qc_dir, "hemoglobin_genes.csv"))