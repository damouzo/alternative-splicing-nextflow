/*
 * PEGASAS_CROSS_CONTRAST — cross-contrast comparison of significant pathways.
 *
 * Collects every <comp_id>_sig_pathways.tsv produced by PEGASAS_COLLATE and builds
 * an UpSet plot + a pathway x contrast heatmap of significantly correlated
 * splicing events. Standalone output (not part of any per-contrast report),
 * published to report/contrast_comparison/pegasas/ so other cross-contrast
 * comparisons can be added there later.
 */
process PEGASAS_CROSS_CONTRAST {
    tag "cross_contrast"
    label 'process_low'

    // Container resolved from modules.config (report image: ggplot2 + UpSetR).

    publishDir "${params.outdir}/report/contrast_comparison/pegasas", mode: params.publish_dir_mode,
        saveAs: { f -> f.startsWith('cross_contrast/') ? f.substring('cross_contrast/'.length()) : f }

    input:
    path sig_files

    output:
    path "cross_contrast/", emit: results
    path "versions.yml",    emit: versions

    script:
    """
    mkdir -p cross_contrast
    pegasas_cross_contrast.R cross_contrast

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
    END_VERSIONS
    """
}
