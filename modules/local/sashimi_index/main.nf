/*
 * ========================================================================================
 *  SASHIMI_INDEX: navigable PDF index for sashimi plots (Fase D)
 * ========================================================================================
 *  Index entries point at the published outdir paths so the table stays valid
 *  after the run finishes (addresses PDF discoverability).
 */

process SASHIMI_INDEX {
    tag "$comparison_id"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}/plots/sashimi",
        mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    input:
    tuple val(comparison_id), path(sashimi_dir)
    // Script as input so content edits invalidate the cache on -resume
    // (plain `${projectDir}/bin/...` references are not hashed).
    path index_builder

    output:
    tuple val(comparison_id), path("sashimi_index.tsv"), emit: index
    path "versions.yml", emit: versions

    script:
    """
    python3 ${index_builder} \\
        --comparison-id  "${comparison_id}" \\
        --sashimi-dir    ${sashimi_dir} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}