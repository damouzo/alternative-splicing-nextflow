/*
 * PEGASAS_PREPARE_CONTRAST — per-contrast PSI matrix + 2-group order.
 *
 * Builds the contrast-specific inputs for the correlation step:
 *   PSI_bySample.tsv  — SE events x samples (IncLevel1/IncLevel2)
 *   group_order.txt   — comma-separated, the 2 contrast groups only
 *
 * Everything is scoped to g1_ids + g2_ids so the score matrix and PSI matrix
 * share exactly the same samples (fixes the upstream KeyError that broke the
 * whole correlation branch).
 */
process PEGASAS_PREPARE_CONTRAST {
    tag "$comparison_id"
    label 'process_low'

    container 'local/pegasas:latest'

    input:
    tuple val(comparison_id),
          path(rmats_se),
          path(b1_samples),   // column->sample order as written by RMATS_POST
          path(b2_samples),
          path(group_info),
          val(g1_ids),
          val(g2_ids)

    output:
    tuple val(comparison_id),
          path("pegasas_inputs/PSI_bySample.tsv"),
          path("pegasas_inputs/group_order.txt"),
          path("pegasas_inputs/contrast_samples.txt"),
          emit: results
    path "versions.yml", emit: versions

    script:
    def g1_arg = g1_ids instanceof List ? g1_ids.join(',') : g1_ids
    def g2_arg = g2_ids instanceof List ? g2_ids.join(',') : g2_ids
    def fdr_flag  = params.pegasas_fdr_cutoff  > 0 ? "--fdr-cutoff  ${params.pegasas_fdr_cutoff}"  : ""
    def dpsi_flag = params.pegasas_dpsi_cutoff > 0 ? "--dpsi-cutoff ${params.pegasas_dpsi_cutoff}" : ""
    // Order-identity guard: the Python script uses the materialized bN_samples.txt
    // files to verify the channel-provided ids match the rMATS column order,
    // aborting on any mismatch instead of silently misassigning PSI columns.
    // The files are guaranteed to exist (checked upstream in PEGASAS_ANALYSIS);
    // the conditional keeps the script usable standalone without them.
    def s1_flag = b1_samples ? "--g1-samples-file ${b1_samples}" : ""
    def s2_flag = b2_samples ? "--g2-samples-file ${b2_samples}" : ""
    """
    ${projectDir}/bin/prepare_pegasas_inputs.py \\
        ${rmats_se} \\
        ${group_info} \\
        --g1-ids "${g1_arg}" \\
        --g2-ids "${g2_arg}" \\
        ${s1_flag} ${s2_flag} \\
        --out-dir pegasas_inputs/ \\
        ${fdr_flag} ${dpsi_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
