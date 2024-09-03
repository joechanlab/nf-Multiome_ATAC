process ATAC_GETMARKERS {
    label 'process_medium'
    conda '/usersoftware/chanj3/ArchR'
    publishDir "${params.outdir}/atac_getmarkers/", mode: 'copy'

    input:
    val name
    path atac_dir
    val rna_h5ad

    output:
    val name, emit: name
    path "${name}_markers", emit: markers_dir

    script:
    """
    Rscript ${baseDir}/bin/atac_GetMarkers.R \
        -s ${name} \
        -a ${atac_dir} \
        -d ${rna_h5ad} \
        -g ${params.getmarker.group} \
        -o ${name}_markers
    """
}
