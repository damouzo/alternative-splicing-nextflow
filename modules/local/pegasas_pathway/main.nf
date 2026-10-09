/*
 * PEGASAS_PATHWAY — shared KS pathway enrichment (computed once, all samples).
 *
 * KS enrichment is a deterministic function of (sample expression, gene set),
 * independent of the contrast, so this runs once for the whole cohort and the
 * per-contrast correlation subsets the resulting scores.
 */
process PEGASAS_PATHWAY {
    tag "shared"
    label 'process_high'

    container 'local/pegasas:latest'

    // Shared, cohort-level output (not tied to a contrast) → raw/pegasas/_shared/
    publishDir "${params.outdir}/raw/${params.tool_ids.pegasas}/_shared", mode: params.publish_dir_mode,
        saveAs: { f -> (f.toString() == 'versions.yml') ? null :
                      (f.toString().startsWith('pathway_out/') ? f.toString().substring('pathway_out/'.length()) : f) }

    input:
    path gene_exp
    path group_info
    path gmt_file
    // Run the repo copy of PEGASAS: the published image predates the per-sample
    // min-TPM filter, so -m is unknown to the baked-in CLI. Staging the source
    // also invalidates the cache on -resume when the package changes.
    path pegasas_src

    output:
    path "pathway_out/", emit: results
    path "versions.yml", emit: versions

    script:
    """
    export PYTHONPATH="${pegasas_src}:\${PYTHONPATH:-}"

    python3 ${pegasas_src}/bin/PEGASAS pathway \\
        ${gene_exp} \\
        ${gmt_file} \\
        ${group_info} \\
        -o pathway_out/ \\
        -n ${params.pegasas_num_interval} \\
        -m ${params.pegasas_min_tpm}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PEGASAS: \$(python3 ${pegasas_src}/bin/PEGASAS --version 2>&1 | head -1)
    END_VERSIONS
    """
}
