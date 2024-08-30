process ATAC_CHROMVAR {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'

    publishDir "${params.outdir}/atac_chromvar/", mode: 'copy'

    input:
    val name
    path arrow_dir
    path ischip_pkl

    output:
    val name
    path "${name}_chromvar", emit: out_dir

    script:
    """
    Rscript ${baseDir}/bin/atac_ChromVAR.R \
        -a ${arrow_dir} \
        -i ${ischip_pkl} \
        -o ${name}_chromvar \
        -p ${params.chromvar.python_path}
    """
}
