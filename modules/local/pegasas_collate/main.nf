/*
 * PEGASAS_COLLATE — per-contrast gather + pathway score collation.
 *
 * Produces the single per-contrast pegasas directory handed to the report:
 *   pegasas_out/all_pathways_scores.tsv  — long (pathway, sample, group, KS, p, median_rank), contrast subset
 *   pegasas_out/correlation_out/         — gathered correlation results
 * and, separately, the cross-contrast summary:
 *   <comp_id>_sig_pathways.tsv           — pathway, n_sig_events (for the upset/heatmap)
 */
process PEGASAS_COLLATE {
    tag "$comparison_id"
    label 'process_low'

    container 'local/pegasas:latest'

    // raw/ layer — unpublished in 'core' mode, except *_sig_pathways.tsv which
    // feeds the pegasas master (n_sig_events) and must survive at every level.
    // versions.yml goes to run_info/.
    publishDir "${params.outdir}/raw/${params.tool_ids.pegasas}/${comparison_id}", mode: params.publish_dir_mode,
        saveAs: { f ->
            def name = f.toString()
            if (name == 'versions.yml') {
                return null
            }
            if (params.publish_level == 'core' || !params.publish_raw) {
                return name.endsWith('_sig_pathways.tsv') ? f : null
            }
            f.startsWith('pegasas_out/') ? f.substring('pegasas_out/'.length()) : f
        }

    input:
    tuple val(comparison_id),
          path(correlation_out),
          path(shared_scores),
          val(g1_ids),
          val(g2_ids)
    // Script as input so content edits invalidate the cache on -resume
    path collate_script

    output:
    tuple val(comparison_id), path("pegasas_out/"), emit: results
    tuple val(comparison_id), path("${comparison_id}_sig_pathways.tsv"), emit: sig_pathways
    path "versions.yml", emit: versions

    script:
    def g1_arg = g1_ids instanceof List ? g1_ids.join(',') : g1_ids
    def g2_arg = g2_ids instanceof List ? g2_ids.join(',') : g2_ids
    """
    mkdir -p pegasas_out
    python3 ${collate_script} \\
        ${shared_scores} \\
        ${correlation_out} \\
        --g1-ids "${g1_arg}" \\
        --g2-ids "${g2_arg}" \\
        --comp-id "${comparison_id}" \\
        --fdr-cutoff ${params.report_fdr_cutoff} \\
        --out-dir pegasas_out \\
        --sig-out .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
