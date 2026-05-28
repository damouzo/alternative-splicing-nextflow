process PEGASAS_CORRELATION {
    tag "$comparison_id"
    label 'process_medium'

    container 'local/pegasas:latest'

    publishDir "${params.outdir}/pegasas/${comparison_id}/correlation", mode: params.publish_dir_mode

    input:
    // pathway_out/ contains one subdirectory per GMT gene set, each with *.scores.txt
    tuple val(comparison_id),
          path(pathway_out),
          path(psi_matrix),
          path(group_order)

    output:
    tuple val(comparison_id), path("correlation_out/"), emit: results
    path "versions.yml",                                emit: versions

    script:
    """
    mkdir -p correlation_out

    # Run correlation for each gene set scores file produced by the pathway step
    for SCORES in \$(find pathway_out/ -name '*.scores.txt'); do
        PEGASAS correlation \\
            "\$SCORES" \\
            ${psi_matrix} \\
            ${group_order} \\
            -o correlation_out/ \\
        || echo "[WARN] Correlation failed for \$SCORES — continuing"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PEGASAS: \$(PEGASAS --version 2>&1 | head -1)
    END_VERSIONS
    """
}
