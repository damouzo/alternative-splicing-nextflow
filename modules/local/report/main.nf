process RENDER_REPORT {
    tag "$comparison_id"
    label 'process_medium'
    
    container 'rocker/verse:4.3.2'
    
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode
    
    input:
    val comparison_id
    path rmats_dir
    path majiq_dir
    path isar_dir
    path report_rmd
    
    output:
    path "${comparison_id}_splicing_report.html", emit: html
    path "versions.yml"                         , emit: versions
    
    script:
    def rmats_arg = rmats_dir.name != 'NO_RMATS' ? rmats_dir : 'NULL'
    def majiq_arg = majiq_dir.name != 'NO_MAJIQ' ? majiq_dir : 'NULL'
    def isar_arg  = isar_dir.name  != 'NO_ISAR'  ? isar_dir  : 'NULL'
    
    """
    # Install required R packages if not present
    Rscript -e "
    if (!require('ggplot2')) install.packages('ggplot2', repos='https://cloud.r-project.org')
    if (!require('dplyr')) install.packages('dplyr', repos='https://cloud.r-project.org')
    if (!require('tidyr')) install.packages('tidyr', repos='https://cloud.r-project.org')
    if (!require('DT')) install.packages('DT', repos='https://cloud.r-project.org')
    "
    
    # Render R Markdown report
    Rscript -e "
    rmarkdown::render(
      input = '${report_rmd}',
      output_format = rmarkdown::html_document(
        toc = TRUE,
        toc_float = TRUE,
        toc_depth = 3,
        code_folding = 'hide',
        theme = 'flatly'
      ),
      output_file = '${comparison_id}_splicing_report.html',
      params = list(
        comparison_id = '${comparison_id}',
        rmats_dir = ${rmats_arg == 'NULL' ? 'NULL' : "'${rmats_arg}'"},
        majiq_dir = ${majiq_arg == 'NULL' ? 'NULL' : "'${majiq_arg}'"},
        isar_dir = ${isar_arg == 'NULL' ? 'NULL' : "'${isar_arg}'"},
        fdr_cutoff = ${params.report_fdr_cutoff},
        dpsi_cutoff = ${params.report_dpsi_cutoff}
      )
    )
    "
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
    END_VERSIONS
    """
}
