process RENDER_REPORT {
    tag "$comparison_id"
    label 'process_medium'

    // Container resolved from modules.config (params.report_container or default)

    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    // stageAs gives each dir a unique name in the work dir — avoids basename collision
    // when rMATS, MAJIQ, ISAR, sashimi, PEGASAS, and LeafCutter all emit a directory per comparison.
    // Original names are passed as vals so NO_* sentinels can still be detected by the Rmd.
    input:
    tuple val(comparison_id),
          val(rmats_name),       path(rmats_dir,       stageAs: 'rmats_in'),
          val(majiq_name),       path(majiq_dir,       stageAs: 'majiq_in'),
          val(isar_name),        path(isar_dir,        stageAs: 'isar_in'),
          val(sashimi_name),     path(sashimi_dir,     stageAs: 'sashimi_in'),
          val(pegasas_name),     path(pegasas_dir,     stageAs: 'pegasas_in'),
          val(leafcutter_name),  path(leafcutter_dir,  stageAs: 'leafcutter_in'),
          val(group1_sample_ids),
          val(group2_sample_ids)
    path report_rmd

    output:
    tuple val(comparison_id), path("${comparison_id}_splicing_report.html"), emit: html
    path "versions.yml",                          emit: versions

    script:
    def rmats_arg      = rmats_name      != 'NO_RMATS'      ? 'rmats_in'      : 'NULL'
    def majiq_arg      = majiq_name      != 'NO_MAJIQ'      ? 'majiq_in'      : 'NULL'
    def isar_arg       = isar_name       != 'NO_ISAR'       ? 'isar_in'       : 'NULL'
    def sashimi_arg    = sashimi_name    != 'NO_SASHIMI'    ? 'sashimi_in'    : 'NULL'
    def pegasas_arg    = pegasas_name    != 'NO_PEGASAS'    ? 'pegasas_in'    : 'NULL'
    def leafcutter_arg = leafcutter_name != 'NO_LEAFCUTTER' ? 'leafcutter_in' : 'NULL'

    // de_results is an optional external path — passed as a string so
    // the Rmd can use file.exists() without staging it into the work directory
    def de_arg  = params.de_results         ? "\"${params.de_results}\""         : 'null'
    def group1_ids_json = groovy.json.JsonOutput.toJson(group1_sample_ids ?: [])
    def group2_ids_json = groovy.json.JsonOutput.toJson(group2_sample_ids ?: [])
    // Per-tool reliability record (known issues) — mirrors the deliverables manifest
    def tool_reliability_json = groovy.json.JsonOutput.toJson(params.tool_reliability ?: [])

    """
    # Write report params to JSON — avoids shell injection from paths with special chars
    cat > report_params.json << 'JSONEOF'
    {
      "comparison_id":          "${comparison_id}",
      "rmats_dir":              ${rmats_arg      == 'NULL' ? 'null' : '"' + rmats_arg      + '"'},
      "majiq_dir":              ${majiq_arg      == 'NULL' ? 'null' : '"' + majiq_arg      + '"'},
      "isar_dir":               ${isar_arg       == 'NULL' ? 'null' : '"' + isar_arg       + '"'},
      "sashimi_dir":            ${sashimi_arg    == 'NULL' ? 'null' : '"' + sashimi_arg    + '"'},
      "pegasas_dir":            ${pegasas_arg    == 'NULL' ? 'null' : '"' + pegasas_arg    + '"'},
      "leafcutter_dir":         ${leafcutter_arg == 'NULL' ? 'null' : '"' + leafcutter_arg + '"'},
      "fdr_cutoff":             ${params.report_fdr_cutoff},
      "dpsi_cutoff":            ${params.report_dpsi_cutoff},
      "majiq_prob_threshold":   ${params.majiq_probability_threshold},
      "majiq_dpsi_cutoff":      ${params.majiq_delta_psi_threshold},
      "sashimi_png_dpi":        ${params.sashimi_png_dpi},
      "organism":               "${params.organism}",
      "run_de_as":              ${params.run_de_as ? 'true' : 'false'},
      "de_results":             ${de_arg},
      "group1_sample_ids":      ${group1_ids_json},
      "group2_sample_ids":      ${group2_ids_json},
      "publish_deliverables":   ${params.publish_deliverables ? 'true' : 'false'},
      "publish_level":          "${params.publish_level}",
      "results_contract_version": "${params.results_contract_version}",
      "tool_reliability":       ${tool_reliability_json}
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
