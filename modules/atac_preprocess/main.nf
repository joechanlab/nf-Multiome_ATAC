process ATAC_PREPROCESS {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/atac_preprocess/", mode: 'copy'

    input:
    tuple val(name), path(rna_h5ad_path), path(atac_h5ad_path), path(fragment_path), path(fragment_index_path)

    output:
    val name, emit: name
    path "${name}_atac_archr", emit: output_dir

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    Rscript ${baseDir}/bin/atac_preprocessing.R \
        -f ${fragment_path} \
        -a ${rna_h5ad_path} \
        -s ${name} \
        -o ${name}_atac_archr
    """
}
