/*
 * ========================================================================================
 *  EXPORT_METADATA: global publication metadata (Fase A)
 * ========================================================================================
 *  Writes the deliverables layer metadata:
 *    run_manifest.yaml  — pipeline/params/tools/reliability record for the run
 *    sample_index.tsv   — sample x comparison traceability
 *    tools_matrix.tsv   — tool enablement + container record
 */

process EXPORT_METADATA {
    tag "deliverables-metadata"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/metadata", mode: params.publish_dir_mode

    input:
    val payload_json    // JSON string: {samples: [...], pipeline: {...}} built in the workflow

    output:
    path "run_manifest.yaml", emit: run_manifest
    path "sample_index.tsv",  emit: sample_index
    path "tools_matrix.tsv",  emit: tools_matrix
    path "versions.yml",      emit: versions

    script:
    // base64 keeps arbitrary characters in sample ids/paths out of the shell
    def payload_b64 = payload_json.bytes.encodeBase64().toString()

    """
    echo '${payload_b64}' | base64 -d > payload.json
    python3 ${projectDir}/bin/build_results_manifest.py \\
        --mode run \\
        --payload-file payload.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}