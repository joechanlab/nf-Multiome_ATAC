process ATAC_NETWORK_ANALYSIS {
    label 'process_medium'
    conda "/usersoftware/chanj3/SnapATAC2"
    publishDir "${params.outdir}/atac_network_analysis/", mode: 'copy'

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
        ${name} \
        ${atac_dir}/${name}_snapatac2.h5ad \
        ${atac_dir}/${name}_motif_enrichment.pkl
    """
}
