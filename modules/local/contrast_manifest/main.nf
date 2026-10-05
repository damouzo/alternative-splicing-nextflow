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

    // Manifests are run metadata, not deliverables: publish to results/run_info/
    // (sibling of raw/), outside the shipped deliverables/ tree.
    publishDir "${params.outdir}/run_info/contrast_manifests",
        mode: params.publish_dir_mode,
        saveAs: { f ->
            f.toString() == 'versions.yml' ? null :
                "${comparison_id}_contrast_manifest.yaml"
        }

    input:
    tuple val(comparison_id),
          val(pipeline_json),
          val(outdir_root),
          val(deliverables_root),
          val(group1_name),
          val(group2_name),
          val(effective_test_method)
    // Script as input so content edits invalidate the cache on -resume
    // (plain `${projectDir}/bin/...` references are not hashed).
    path manifest_builder

    output:
    tuple val(comparison_id), path("contrast_manifest.yaml"), emit: manifest
    path "versions.yml",           emit: versions

    script:
    def payload = groovy.json.JsonOutput.toJson([
        comparison_id:    comparison_id,
        group1_name:      group1_name,
        group2_name:      group2_name,
        effective_test_method: effective_test_method,
        pipeline:         new groovy.json.JsonSlurper().parseText(pipeline_json),
        outdir:           outdir_root,
        deliverables_root: deliverables_root
    ])
    def payload_b64 = payload.bytes.encodeBase64().toString()

    """
    echo '${payload_b64}' | base64 -d > payload.json
    python3 ${manifest_builder} \\
        --mode contrast_disk \\
        --payload-file payload.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}