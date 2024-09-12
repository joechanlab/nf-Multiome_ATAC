process ATAC_NETWORK_ANALYSIS {
    label 'process_medium'
    conda "/usersoftware/chanj3/SnapATAC2"
    publishDir "${params.outdir}/atac_network_analysis/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path atac_dir

    output:
    val name, emit: name
    path "${name}", emit: output_dir

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/atac_network_analysis.py \
        ${atac_dir}/${name}_peak_mtx.h5ad \
        ${atac_dir}/${name}_gene_mtx.h5ad \
        ${atac_dir}/${name}_diff_peaks.csv \
        ${atac_dir}/${name}_motif_enrichment.pkl \
        --output_dir ${name}
    """
}
