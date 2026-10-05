/*
 * ========================================================================================
 *  EXPORT_METADATA: global publication metadata (Fase A)
 * ========================================================================================
 *  Writes the shippable run metadata:
 *    run_manifest.yaml  — pipeline/params/tools/reliability record for the run
 *    sample_index.tsv   — sample x comparison traceability (no cluster paths)
 *    sample_paths.tsv   — absolute input paths, routed to raw/_internal/ only
 *    README.md          — layout + how to read the TSVs
 */

process EXPORT_METADATA {
    tag "run-info"
    label 'process_low'

    // Shippable run metadata. sample_paths.tsv carries the absolute cluster
    // paths and is routed to raw/_internal/ only (never shipped).
    // run_manifest.yaml lives outside deliverables/ (results/run_info/) so the
    // shipment carries no internal run metadata; it is audit/QA only.
    publishDir "${params.outdir}/deliverables/run_info",
        mode: params.publish_dir_mode,
        saveAs: { f -> (f.toString() in ['sample_paths.tsv', 'run_manifest.yaml', 'versions.yml']) ? null : f }
    publishDir "${params.outdir}/run_info",
        mode: params.publish_dir_mode,
        saveAs: { f -> (f.toString() == 'run_manifest.yaml') ? f : null }
    publishDir "${params.outdir}/raw/_internal",
        mode: params.publish_dir_mode,
        saveAs: { f -> (f.toString() == 'sample_paths.tsv' &&
                        params.publish_level != 'core' && params.publish_raw) ? f : null }

    input:
    val payload_json    // JSON string: {samples: [...], pipeline: {...}} built in the workflow
    // Script as input so content edits invalidate the cache on -resume
    // (plain `${projectDir}/bin/...` references are not hashed).
    path manifest_builder

    output:
    path "run_manifest.yaml", emit: run_manifest
    path "sample_index.tsv",  emit: sample_index
    path "sample_paths.tsv",  emit: sample_paths
    path "README.md",         emit: readme
    path "versions.yml",      emit: versions

    script:
    // base64 keeps arbitrary characters in sample ids/paths out of the shell
    def payload_b64 = payload_json.bytes.encodeBase64().toString()

    """
    echo '${payload_b64}' | base64 -d > payload.json
    python3 ${manifest_builder} \\
        --mode run \\
        --payload-file payload.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}