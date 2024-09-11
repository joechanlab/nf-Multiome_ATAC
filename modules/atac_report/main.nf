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
    python -m ipykernel install --user --name snapatac2 --display-name "Python (snapatac2)"
    papermill ${baseDir}/bin/atac_report.ipynb ${name}_report.ipynb -p path ${name} -p plots ${params.report.plots}
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
