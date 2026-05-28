process PEGASAS_PREPARE {
    tag "$comparison_id"
    label 'process_low'

    input:
    tuple val(comparison_id),
          path(salmon_tpm),
          path(rmats_se),
          path(group_info)

    output:
    tuple val(comparison_id),
          path("pegasas_inputs/gene_exp_bySample.tsv"),
          path("pegasas_inputs/PSI_bySample.tsv"),
          path("pegasas_inputs/group_info.tsv"),
          path("pegasas_inputs/group_order.txt"),
          emit: inputs
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/prepare_pegasas_inputs.py \\
        ${salmon_tpm} \\
        ${rmats_se} \\
        ${group_info} \\
        --out-dir pegasas_inputs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
