process ATAC_PEAK2GENE {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/atac_peak2gene/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path atac_dir
    val rna_h5
    val rna_h5ad

    output:
    val name, emit: name
    path "${name}_peak2gene", emit: out_dir
    val rna_h5, emit: rna_h5
    val rna_h5ad, emit: rna_h5ad

    script:
    """
    Rscript ${baseDir}/bin/atac_Peak2Genes.R \
        -s ${name} \
        -a ${atac_dir} \
        -r ${rna_h5} \
        -d ${rna_h5ad} \
        -o ${name}_peak2gene
    """
}
