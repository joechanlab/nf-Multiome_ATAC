process ATAC_PREPROCESS {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/atac_preprocess/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), val(rna_h5), val(rna_h5ad), val(rna_seacells_h5ad), val(rna_seacells_dir), val(atac_h5ad), val(fragment_path), val(fragment_index_path)

    output:
    val name, emit: name
    path "${name}_atac_archr", emit: output_dir
    val rna_seacells_h5ad, emit: rna_seacells_h5ad
    val rna_seacells_dir, emit: rna_seacells_dir
    val rna_h5, emit: rna_h5
    val rna_h5ad, emit: rna_h5ad

    script:
    """
    Rscript ${baseDir}/bin/atac_preprocessing.R \
        -f ${fragment_path} \
        -a ${rna_h5ad} \
        -s ${name} \
        -o ${name}_atac_archr \
        -m ${params.macs3}
    """
}
