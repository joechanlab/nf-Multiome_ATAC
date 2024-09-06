process ATAC_PREPROCESS_SNAPATAC2 {
    label 'process_medium'
    conda "/usersoftware/chanj3/SnapATAC2"
    publishDir "${params.outdir}/atac_preprocess_snapatac2/", mode: 'copy'

    input:
    tuple val(name), val(rna_h5), val(rna_h5ad), val(rna_seacells_h5ad), val(rna_seacells_dir), val(atac_h5ad), val(fragment_path), val(fragment_index_path)

    output:
    val name, emit: name
    path "${name}", emit: output_dir
    val rna_h5ad, emit: rna_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/atac_preprocessing_snapATAC.py \
        -f ${fragment_path} \
        -a ${rna_h5ad} \
        -s ${name} \
        -o ${name}
    """
}
