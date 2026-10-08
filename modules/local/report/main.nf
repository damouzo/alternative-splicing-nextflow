process RENDER_REPORT {
    tag "$comparison_id"
    label 'process_medium'

    // Container resolved from modules.config (params.report_container or default)

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}", mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    // stageAs gives each dir a unique name in the work dir — avoids basename collision
    // when rMATS, MAJIQ, ISAR, sashimi, PEGASAS, and LeafCutter all emit a directory per comparison.
    // Original names are passed as vals so NO_* sentinels can still be detected by the Rmd.
    input:
    tuple val(comparison_id),
          val(rmats_name),       path(rmats_dir,       stageAs: 'rmats_in'),
          val(majiq_name),       path(majiq_dir,       stageAs: 'majiq_in'),
          val(isar_name),        path(isar_dir,        stageAs: 'isar_in'),
          val(sashimi_name),     path(sashimi_dir,     stageAs: 'sashimi_in'),
          path(sashimi_index,    stageAs: 'sashimi_index_input.tsv'),
          val(pegasas_name),     path(pegasas_dir,     stageAs: 'pegasas_in'),
          val(leafcutter_name),  path(leafcutter_dir,  stageAs: 'leafcutter_in'),
          val(group1_sample_ids),
          val(group2_sample_ids),
          val(group1_name),
          val(group2_name),
          val(effective_test_method)
    path report_rmd

    output:
    tuple val(comparison_id), path("${comparison_id}_splicing_report.html"), emit: html
    path "plots",                                 emit: plots
    path "versions.yml",                          emit: versions

    script:
    def rmats_arg      = rmats_name      != 'NO_RMATS'      ? 'rmats_in'      : 'NULL'
    def majiq_arg      = majiq_name      != 'NO_MAJIQ'      ? 'majiq_in'      : 'NULL'
    def isar_arg       = isar_name       != 'NO_ISAR'       ? 'isar_in'       : 'NULL'
    def sashimi_arg    = sashimi_name    != 'NO_SASHIMI'    ? 'sashimi_in'    : 'NULL'
    def pegasas_arg    = pegasas_name    != 'NO_PEGASAS'    ? 'pegasas_in'    : 'NULL'
    def leafcutter_arg = leafcutter_name != 'NO_LEAFCUTTER' ? 'leafcutter_in' : 'NULL'

    // dge_dirs is a list of optional external DGE directories, passed as JSON so
    // the Rmd can file.exists()/list.dirs() them without staging into work.
    def dge_dirs = []
    if (params.dge_dirs) {
        def v = params.dge_dirs
        dge_dirs = (v instanceof List || v instanceof Object[])
            ? v.collect { it.toString() }
            : v.toString().split(',')*.trim()
    }
    if (params.de_results) { dge_dirs << params.de_results.toString() }
    dge_dirs = dge_dirs.findAll { it }.unique()
    def dge_dirs_json = groovy.json.JsonOutput.toJson(dge_dirs)
    // Deprecated single-dir alias kept for backward compatibility.
    def de_arg  = params.de_results ? groovy.json.JsonOutput.toJson(params.de_results.toString()) : 'null'
    def de_map_json = groovy.json.JsonOutput.toJson(params.de_results_map ?: [:])
    def group1_ids_json = groovy.json.JsonOutput.toJson(group1_sample_ids ?: [])
    def group2_ids_json = groovy.json.JsonOutput.toJson(group2_sample_ids ?: [])

    """
    # Write report params to JSON — avoids shell injection from paths with special chars
    cat > report_params.json << 'JSONEOF'
    {
      "comparison_id":          "${comparison_id}",
      "plots_dir":              "plots",
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
      "dge_dirs":               ${dge_dirs_json},
      "de_results":             ${de_arg},
      "de_results_map":         ${de_map_json},
      "group1_sample_ids":      ${group1_ids_json},
      "group2_sample_ids":      ${group2_ids_json},
      "group1_name":            "${group1_name}",
      "group2_name":            "${group2_name}",
      "effective_test_method":  "${effective_test_method}",
      "report_min_reads":       ${params.report_min_reads},
      "publish_deliverables":   ${params.publish_deliverables ? 'true' : 'false'},
      "publish_level":          "${params.publish_level}",
      "results_contract_version": "${params.results_contract_version}"
    }
    JSONEOF

    # Declared output: the report writes each plot here as a rendering side
    # effect (params\$plots_dir); created up front so the output always exists,
    # even for contrasts where every plot is skipped.
    mkdir -p plots

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

    # Fold the sashimi deliverables (PDFs + index) into plots/ so this process
    # stays the only publisher of the plots/ tree. A separate publisher under
    # plots/ gets wiped when Nextflow re-creates the parent directory on publish.
    if [ "${sashimi_name}" != "NO_SASHIMI" ]; then
        for ETYPE in SE A5SS A3SS MXE RI; do
            if compgen -G "sashimi_in/\${ETYPE}/Sashimi_plot/*.pdf" >/dev/null; then
                mkdir -p "plots/rmats/sashimi/\${ETYPE}"
                cp -L sashimi_in/\${ETYPE}/Sashimi_plot/*.pdf "plots/rmats/sashimi/\${ETYPE}/"
            elif compgen -G "sashimi_in/\${ETYPE}/*.pdf" >/dev/null; then
                mkdir -p "plots/rmats/sashimi/\${ETYPE}"
                cp -L sashimi_in/\${ETYPE}/*.pdf "plots/rmats/sashimi/\${ETYPE}/"
            fi
        done
        if [ -f sashimi_index_input.tsv ]; then
            mkdir -p plots/rmats/sashimi
            cp -L sashimi_index_input.tsv plots/rmats/sashimi/sashimi_index.tsv
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        rmarkdown: \$(Rscript -e "cat(as.character(packageVersion('rmarkdown')))")
    END_VERSIONS
    """
}
