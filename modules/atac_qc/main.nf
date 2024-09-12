process ATAC_QC {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/muon.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/atac_qc/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), val(rna_h5), val(rna_h5ad), val(rna_seacells_h5ad), val(rna_seacells_dir), val(atac_h5ad), val(fragment_path), val(fragment_index_path)

    output:
    val name, emit: name
    path "${name}_atac_qc.h5ad", emit: out_h5ad

    script:
    """
    export PATH=/opt/conda/envs/muon/bin/:\$PATH
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/atac_qc.py \
        ${atac_h5ad} \
        ${fragment_path} \
        ${name}_atac_qc.h5ad
    """
}
