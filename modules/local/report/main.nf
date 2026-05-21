process RENDER_REPORT {
    tag "$comparison_id"
    label 'process_medium'
    
    // Container resolved from modules.config (params.report_container or ghcr.io default)
    
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
