process ATAC_QC {
    label 'gpus'
    container 'library://mamie_wang/nf-scrnaseq/muon.sif:latest'
    containerOptions '--nv'
    publishDir "${params.outdir}/atac_qc/", mode: 'copy'

    input:
    tuple val(name), path(h5ad_path), path(fragment_path), path(fragment_index_path)

    output:
    val "${name}", emit: name
    path "${name}_atac_qc.h5ad", emit: output_h5ad

    script:
    """
    export PATH=/opt/conda/envs/muon/bin/:$PATH
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/atac_qc.py \
        ${h5ad_path} \
        ${fragment_path} \
        ${name}_atac_qc.h5ad
    """
}
