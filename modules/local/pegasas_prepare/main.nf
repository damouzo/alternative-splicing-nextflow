process PEGASAS_PREPARE {
    tag "$comparison_id"
    label 'process_low'

    input:
    tuple val(comparison_id),
          path(salmon_tpm),
          path(rmats_se),
          path(group_info),
          val(g1_ids),
          val(g2_ids)

    output:
    tuple val(comparison_id),
          path("pegasas_inputs/gene_exp_bySample.tsv"),
          path("pegasas_inputs/PSI_bySample.tsv"),
          path("pegasas_inputs/group_info.tsv"),
          path("pegasas_inputs/group_order.txt"),
          emit: inputs
    path "versions.yml", emit: versions

    script:
    // Pass sample_ids as comma-separated lists; the script splits and validates
    // them against the per-event column counts in SE.MATS.JC.txt.
    def g1_arg = g1_ids instanceof List ? g1_ids.join(',') : g1_ids
    def g2_arg = g2_ids instanceof List ? g2_ids.join(',') : g2_ids
    """
    prepare_pegasas_inputs.py \\
        ${salmon_tpm} \\
        ${rmats_se} \\
        ${group_info} \\
        --g1-ids "${g1_arg}" \\
        --g2-ids "${g2_arg}" \\
        --out-dir pegasas_inputs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
