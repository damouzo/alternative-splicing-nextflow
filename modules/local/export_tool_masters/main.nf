/*
 * ========================================================================================
 *  EXPORT_TOOL_MASTERS: standardised masters for MAJIQ / ISAR / LeafCutter / PEGASAS
 * ========================================================================================
 *  Forms:
 *    - EXPORT_TOOL_MASTERS(comp, tool, tool_dir, tool_params) — per-tool tables
 *    - CROSS_TOOL_MASTER(comp, master_files)                  — gene x tool + gene summary
 */

process EXPORT_TOOL_MASTERS {
    tag "${comparison_id}_${tool}"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}/tables",
        mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    input:
    tuple val(comparison_id), val(tool), path(tool_dir), val(tool_params_json)
    // Script as input so content edits invalidate the cache on -resume
    path exporter_script

    output:
    tuple val(comparison_id),
          val(tool),
          path("${comparison_id}.${tool}.master.tsv"),
          path("${comparison_id}.${tool}.significant.tsv"),
          path("${comparison_id}.${tool}.summary.tsv"),
          emit: tables
    path "versions.yml", emit: versions

    script:
    // base64 keeps tool-specific values (paths with special chars) out of the shell
    def tool_params_b64 = tool_params_json.bytes.encodeBase64().toString()

    """
    python3 ${exporter_script} \\
        --comparison-id "${comparison_id}" \\
        --tool          ${tool} \\
        --input-dir     ${tool_dir} \\
        --tool-params   ${tool_params_b64} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}

process CROSS_TOOL_MASTER {
    tag "$comparison_id"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}/tables",
        mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    input:
    tuple val(comparison_id), path(master_files), val(tool_params_json)
    // Script as input so content edits invalidate the cache on -resume
    path exporter_script

    output:
    tuple val(comparison_id),
          path("${comparison_id}.cross_tool.master.tsv"),
          path("${comparison_id}.cross_tool.gene_summary.tsv"),
          emit: master
    path "versions.yml",                           emit: versions

    script:
    def tool_params_b64 = tool_params_json.bytes.encodeBase64().toString()
    """
    python3 ${exporter_script} \\
        --comparison-id "${comparison_id}" \\
        --tool cross_tool \\
        --master-files ${master_files} \\
        --tool-params   ${tool_params_b64} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}