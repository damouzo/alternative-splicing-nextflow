/*
 * ========================================================================================
 *  CONTRAST_MANIFEST: per-comparison deliverable manifest (Fase A/B)
 * ========================================================================================
 *  Lists every CORE deliverable of one comparison (data tables, sashimi index,
 *  report) with paths relative to the run outdir, plus thresholds and the
 *  tool_reliability record (known issues).
 *
 *  The manifest reflects the published deliverables directly (deterministic
 *  layout), which removes any dependency on channel aggregation order.
 *  Ordering is guaranteed upstream by the join-gated trigger channel.
 */

process CONTRAST_MANIFEST {
    tag "$comparison_id"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/contrasts/${comparison_id}/metadata",
        mode: params.publish_dir_mode

    input:
    tuple val(comparison_id),
          val(pipeline_json),
          val(outdir_root),
          val(group1_name),
          val(group2_name)

    output:
    tuple val(comparison_id), path("contrast_manifest.yaml"), emit: manifest
    path "versions.yml",           emit: versions

    script:
    def payload = groovy.json.JsonOutput.toJson([
        comparison_id: comparison_id,
        group1_name:   group1_name,
        group2_name:   group2_name,
        pipeline:      new groovy.json.JsonSlurper().parseText(pipeline_json),
        outdir:        outdir_root
    ])
    def payload_b64 = payload.bytes.encodeBase64().toString()

    """
    echo '${payload_b64}' | base64 -d > payload.json
    python3 ${projectDir}/bin/build_results_manifest.py \\
        --mode contrast_disk \\
        --payload-file payload.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}