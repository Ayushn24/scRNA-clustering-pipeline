process quality_control {

    publishDir 'results/qc_report', mode: 'copy'

    input:
    path zip_file
    path qc_script

    output:
    path 'results/qc_report/*'
    path 'results/qc_report/filtered_matrix.h5ad' ,emit: filtered_input


    script:
    """
    python3 ${qc_script} ${zip_file} .
    """
}

process normalization {

    publishDir 'results/normalized_data', mode: 'copy'

    input:
    
    path filtered_input
    path norm_script

    output:
    path 'results/normalized_data/log_normalized.h5ad', emit: normalized_hvg_input 
    path 'results/normalized_data/*' 

    script:
    """
    python ${norm_script} ${filtered_input} .
    """
}

process dimensionality_reduction {

    publishDir 'results/dimensionality_reduced_data', mode: 'copy'

    input:
    
    path normalized_hvg_input
    path dim_reduc_script

    output:
    path 'results/dimensionality_reduction/adata_umap.h5ad' ,emit: dim_reduced_input
    path 'results/dimensionality_reduction/*'

    script:
    """
    python ${dim_reduc_script} ${normalized_hvg_input} .
    """
}

process clustering {

    publishDir 'results/clustering_data', mode: 'copy'

    input:
    path dim_reduced_input
    path clus_script

    output:
    path 'results/clustering/adata_clustered.h5ad' , emit: clus_input
    path 'results/clustering/*'


    script:
    """
    python ${clus_script} ${dim_reduced_input} .
    """
}

process marker_genes_identification {

    publishDir 'results/marker_genes', mode: 'copy'

    input:
    path clus_input
    path mrkgen_script

    output:
    path 'results/marker_genes/*'
    
    script:
    """
    python ${mrkgen_script} ${clus_input} .
    """
}
  workflow {

    // Input script/data channels
    Channel.fromPath(params.zip_file)            .set { zip_channel }
    Channel.fromPath(params.qc_script)           .set { qc_channel }
    Channel.fromPath(params.norm_script)         .set { norm_channel }
    Channel.fromPath(params.dim_reduc_script)    .set { dimreduc_channel }
    Channel.fromPath(params.clus_script)         .set { clus_channel }
    Channel.fromPath(params.mrkgen_script)       .set { mrkgen_channel }

    // Capture emitted output from each process
    qc_output       = quality_control(zip_channel, qc_channel)
    norm_output     = normalization(qc_output.filtered_input, norm_channel)
    dr_output       = dimensionality_reduction(norm_output.normalized_hvg_input, dimreduc_channel)
    clus_output     = clustering(dr_output.dim_reduced_input, clus_channel)
    marker_genes_identification(clus_output.clus_input, mrkgen_channel)
}
workflow.onComplete {
    println "\n✅ scRNA-seq pipeline completed successfully!"
    println "📁 All results are saved under the 'results/' directory."
    println "🧪 Output includes QC reports, normalized data, UMAP plots, clusters, and marker genes.\n"
}