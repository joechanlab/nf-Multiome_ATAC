process ATAC_INSILICOCHIP {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest'
    publishDir "${params.outdir}/atac_insilicochip/", mode: 'copy'

    input:
    val name
    path arrow_dir
    path rna_seacells_h5ad
    path rna_seacells_dir

    output:
    val name, emit: name
    path arrow_dir, emit: arrow_dir
    path "${name}_tf_binding_scores.ISChIP.csv", emit: output_csv
    path "${name}_tf_peaks.ISChIP_idx.pkl", emit: output_pkl

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
