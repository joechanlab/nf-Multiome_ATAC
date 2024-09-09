process ATAC_REPORT {
    label 'process_medium'
    conda "/usersoftware/chanj3/SnapATAC2"
    publishDir "${params.outdir}/atac_report/", mode: 'copy'

    input:
    val name
    path atac_dir

    output:
    path "${name}_report.ipynb", emit: report_ipynb
    path "${name}_report.html", emit: report_html

    script:
    """
    export HOME=\$PWD
    papermill ${baseDir}/bin/atac_report.ipynb ${name}_report.ipynb
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
