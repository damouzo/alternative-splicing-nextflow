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
        mode: params.publish_dir_mode

    input:
    tuple val(comparison_id), path(sashimi_dir)
    val  publish_root    // absolute outdir — index paths are relative to it

    output:
    tuple val(comparison_id), path("sashimi_index.tsv"), emit: index
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/build_sashimi_index.py \\
        --comparison-id "${comparison_id}" \\
        --sashimi-dir   ${sashimi_dir} \\
        --publish-root  ${publish_root} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}