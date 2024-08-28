process ATAC_QC {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/muon.sif:latest'
    publishDir "${params.outdir}/atac_qc/", mode: 'copy'

    input:
    tuple val(name), path(rna_h5ad_path), path(atac_h5ad_path), path(fragment_path), path(fragment_index_path)

    output:
    val name, emit: name
    path "${name}_atac_qc.h5ad", emit: output_h5ad

    script:
    """
    export PATH=/opt/conda/envs/muon/bin/:\$PATH
    export NUMBA_CACHE_DIR=\$PWD
    python ${projectDir}/bin/atac_qc.py \
        ${atac_h5ad_path} \
        ${fragment_path} \
        ${name}_atac_qc.h5ad
    """
}
