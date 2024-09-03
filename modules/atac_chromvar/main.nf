process ATAC_CHROMVAR {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/atac_chromvar/", mode: 'copy'

    input:
    val name
    path arrow_dir
    path ischip_pkl
    val rna_h5
    val rna_h5ad

    output:
    val name, emit: name
    path "${name}_chromvar", emit: out_dir
    val rna_h5, emit: rna_h5
    val rna_h5ad, emit: rna_h5ad

    script:
    """
    Rscript ${baseDir}/bin/atac_ChromVAR.R \
        -a ${arrow_dir} \
        -i ${ischip_pkl} \
        -o ${name}_chromvar \
        -p ${params.chromvar.python_path}
    """
}
