process RENDER_REPORT {
    tag "$comparison_id"
    label 'process_medium'
    
    // Container resolved from modules.config (params.report_container or ghcr.io default)
    
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode
    
    // stageAs gives each dir a unique name in the work dir — avoids basename collision
    // when rMATS, MAJIQ and ISAR all emit a directory named after the comparison_id.
    // The original names are passed as vals so NO_* placeholders can still be detected.
    input:
    tuple val(comparison_id),
          val(rmats_name), path(rmats_dir, stageAs: 'rmats_in'),
          val(majiq_name), path(majiq_dir, stageAs: 'majiq_in'),
          val(isar_name),  path(isar_dir,  stageAs: 'isar_in')
    path report_rmd
    
    output:
    path "${comparison_id}_splicing_report.html", emit: html
    path "versions.yml"                         , emit: versions
    
    script:
    def rmats_arg = rmats_name != 'NO_RMATS' ? 'rmats_in' : 'NULL'
    def majiq_arg = majiq_name != 'NO_MAJIQ' ? 'majiq_in' : 'NULL'
    def isar_arg  = isar_name  != 'NO_ISAR'  ? 'isar_in'  : 'NULL'
    
    """
    # Write report params to a JSON file to avoid shell injection from comparison_id
    # or directory paths containing quotes / special characters
    cat > report_params.json << 'JSONEOF'
    {
      "comparison_id": "${comparison_id}",
      "rmats_dir":     ${rmats_arg == 'NULL' ? 'null' : '"' + rmats_arg + '"'},
      "majiq_dir":     ${majiq_arg == 'NULL' ? 'null' : '"' + majiq_arg + '"'},
      "isar_dir":      ${isar_arg  == 'NULL' ? 'null' : '"' + isar_arg  + '"'},
      "fdr_cutoff":            ${params.report_fdr_cutoff},
      "dpsi_cutoff":           ${params.report_dpsi_cutoff},
      "majiq_prob_threshold":  ${params.majiq_probability_threshold},
      "majiq_dpsi_cutoff":     ${params.majiq_delta_psi_threshold}
    }
    JSONEOF

    Rscript -e "
    p <- jsonlite::fromJSON('report_params.json')
    rmarkdown::render(
      input = '${report_rmd}',
      output_format = rmarkdown::html_document(
        toc = TRUE,
        toc_float = TRUE,
        toc_depth = 3,
        code_folding = 'hide',
        theme = 'flatly'
      ),
      output_file = p[['comparison_id']]  |> paste0('_splicing_report.html'),
      params = p
    )
    "
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
    END_VERSIONS
    """
}
