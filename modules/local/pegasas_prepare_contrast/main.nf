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
    """
    ${projectDir}/bin/prepare_pegasas_inputs.py \\
        ${rmats_se} \\
        ${group_info} \\
        --g1-ids "${g1_arg}" \\
        --g2-ids "${g2_arg}" \\
        --out-dir pegasas_inputs/ \\
        ${fdr_flag} ${dpsi_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
