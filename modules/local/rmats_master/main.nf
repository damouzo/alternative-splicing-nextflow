/*
 * ========================================================================================
 *  RMATS_MASTER: rMATS CORE deliverables (Fase C)
 * ========================================================================================
 *  Builds <comparison_id>.rmats.{master,significant,summary}.tsv from the raw
 *  rMATS JC tables, flagging the --cstat numeric floor (fdr_floor_flag).
 */

process RMATS_MASTER {
    tag "$comparison_id"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}/tables",
        mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    input:
    tuple val(comparison_id), path(rmats_dir)
    val  fdr_cutoff
    val  dpsi_cutoff
    // Declared as input so its CONTENT is part of the task hash: a plain
    // `${projectDir}/bin/...` reference is not hashed and edits would be
    // silently served from cache on -resume.
    path master_builder

    output:
    tuple val(comparison_id),
          path("${comparison_id}.rmats.master.tsv"),
          path("${comparison_id}.rmats.significant.tsv"),
          path("${comparison_id}.rmats.summary.tsv"),
          emit: tables
    path "versions.yml", emit: versions

    script:
    """
    python3 ${master_builder} \\
        --comparison-id "${comparison_id}" \\
        --rmats-dir     ${rmats_dir} \\
        --fdr-cutoff    ${fdr_cutoff} \\
        --dpsi-cutoff   ${dpsi_cutoff} \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}