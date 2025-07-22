# 🧬 scRNA-seq Clustering Pipeline

This repository contains a modular and reproducible pipeline for single-cell RNA-seq (scRNA-seq) clustering analysis using **Nextflow** and **Scanpy**. The pipeline is designed to be executed step-by-step using Python scripts for each process, orchestrated by Nextflow.

## 📂 Pipeline Overview

The pipeline consists of the following stages:

1. **Quality Control** (`quality_control`)
2. **Normalization** (`normalization`)
3. **Dimensionality Reduction** (`dimensionality_reduction`)
4. **Clustering** (`clustering`)
5. **Marker Gene Identification** (`marker_genes_identification`)

Each step generates its own outputs inside the structured `results/` directory.

## 🚀 Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/scrna-clustering-pipeline.git
cd scrna-clustering-pipeline
```

### 2. Load or Build Docker Image

If you have the Docker image file (`scrna_docker_image.tar`):

```bash
docker load -i scrna_docker_image.tar
```

If you change the image name, update the Nextflow `config` file accordingly:

```groovy
process {
  container = 'your_custom_image_name'
  executor = 'local'
}

docker {
  enabled = true
  pullLatest = true
}
```

### 3. Run the Pipeline

```bash
nextflow run scrna_pipeline.nf \
  --zip_file test.zip \
  --qc_script qc_script.py \
  --norm_script norm_script.py \
  --dim_reduc_script dimreduc_script.py \
  --clus_script clus_script.py \
  --mrkgen_script mrkgen_script.py \
  -with-docker
```

## 📁 Output Structure

- `results/qc_report/` – Quality control report and filtered `.h5ad`
- `results/normalized_data/` – Normalized output files
- `results/dimensionality_reduction/` – UMAP and PCA data
- `results/clustering/` – Clustered AnnData
- `results/marker_genes/` – Marker gene tables

## 🎥 Demo Video

A demo video of the pipeline execution is available in the main branch:  
**`scrna_pipe_demo.mp4`**  
You can open it directly to see how the pipeline runs end-to-end.

## 📚 Citations

- **Nextflow**:  
  _P. Di Tommaso, et al. Nextflow enables reproducible computational workflows._  
  *Nature Biotechnology 35, 316–319 (2017)*  
  [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

- **Scanpy**:  
  _F. Alexander Wolf, Philipp Angerer, Fabian J. Theis._  
  *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology (2018).  
  [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

## 🧪 Authors

Developed by Ayush, 2025. Contributions and feedback welcome!
