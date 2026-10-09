/*
 * ========================================================================================
 *  VALIDATE_RESULTS: results contract QA gate + run_info finalisation
 * ========================================================================================
 *  Runs after the shippable layer is published and:
 *    - merges every per-process versions.yml into run_info/software_versions.yml
 *    - writes run_info/qa_report.txt (structural failures abort the pipeline)
 *
 *  Structural failures (missing shippable deliverables) abort the run; content
 *  warnings are reported but do not block.
 */

process VALIDATE_RESULTS {
    tag "results-contract"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/run_info",
        mode: params.publish_dir_mode,
        saveAs: { f -> f.toString() == 'versions.yml' ? null : f }

    input:
    tuple val(n_contrast_manifests), val(n_global_manifests), val(version_files)
    val  outdir_root
    val  deliverables_root
    // Scripts and schema as inputs so content edits invalidate the cache on
    // -resume (plain `${projectDir}/...` references are not hashed).
    path merge_versions_script
    path contract_validator
    path results_schema

    output:
    path "software_versions.yml", emit: software_versions
    path "qa_report.txt",         emit: qa_report
    path "versions.yml",          emit: versions

    script:
    """
    python3 ${merge_versions_script} \\
        --inputs '${version_files}' \\
        --output software_versions.yml

    # warnings (exit 1) must not abort the run; only structural failures (>= 2) do
    set +e
    python3 ${contract_validator} \\
        --outdir            ${outdir_root} \\
        --deliverables-root ${deliverables_root} \\
        --schema            ${results_schema} \\
        --ignore-missing    software_versions.yml \\
        --report-file       qa_report.txt
    RC=\$?
    set -e
    if [ \$RC -ge 2 ]; then
        echo "[VALIDATE_RESULTS] results contract broken — see qa_report.txt" >&2
        exit \$RC
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
