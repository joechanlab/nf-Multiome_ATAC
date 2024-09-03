process ATAC_INSILICOCHIP {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/atac_insilicochip/", mode: 'copy'

    input:
    val name
    path arrow_dir
    val rna_seacells_h5ad
    val rna_seacells_dir
    val rna_h5
    val rna_h5ad

    output:
    val name, emit: name
    path "${name}_tf_binding_scores.ISChIP.csv", emit: output_csv
    path "${name}_tf_peaks.ISChIP_idx.pkl", emit: output_pkl
    path arrow_dir, emit: arrow_dir
    val rna_h5, emit: rna_h5
    val rna_h5ad, emit: rna_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    export MPLCONFIGDIR=\$PWD
    python ${baseDir}/bin/atac_InSilicoCHIP.py \
        ${name} \
        ${arrow_dir} \
        ${rna_seacells_h5ad} \
        ${rna_seacells_dir} \
        ${name}_tf_binding_scores.ISChIP.csv \
        ${name}_tf_peaks.ISChIP_idx.pkl
    """
}
