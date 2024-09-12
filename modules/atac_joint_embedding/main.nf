process ATAC_JOINT_EMBEDDING {
    label 'process_medium'
    conda "/usersoftware/chanj3/SnapATAC2"
    publishDir "${params.outdir}/atac_joint_embedding/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path atac_dir
    val rna_h5ad

    output:
    val name, emit: name
    path "${name}", emit: output_dir

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/atac_joint_embedding.py \
        ${name} \
        ${atac_dir}/${name}_snapatac2.h5ad \
        ${rna_h5ad}
    """
}
