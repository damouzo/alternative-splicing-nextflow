process PEGASAS_PATHWAY {
    tag "$comparison_id"
    label 'process_high'

    container 'local/pegasas:latest'

    publishDir "${params.outdir}/pegasas/${comparison_id}/pathway", mode: params.publish_dir_mode

    input:
    tuple val(comparison_id),
          path(gene_exp),
          path(psi_matrix),
          path(group_info),
          path(group_order)
    path  gmt_file

    output:
    tuple val(comparison_id),
          path("pathway_out/"),
          path(psi_matrix),
          path(group_order),
          emit: results
    path "versions.yml", emit: versions

    script:
    """
    PEGASAS pathway \\
        ${gene_exp} \\
        ${gmt_file} \\
        ${group_info} \\
        -o pathway_out/ \\
        -n ${params.pegasas_num_interval}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PEGASAS: \$(PEGASAS --version 2>&1 | head -1)
    END_VERSIONS
    """
}
