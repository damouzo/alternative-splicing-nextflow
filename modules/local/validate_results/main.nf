/*
 * ========================================================================================
 *  VALIDATE_RESULTS: results contract QA gate (Fase E)
 * ========================================================================================
 *  Runs validate_results_contract.py against the published outdir. Structural
 *  failures (missing CORE deliverables) abort the pipeline; content warnings
 *  are reported but do not block the run.
 */

process VALIDATE_RESULTS {
    tag "results-contract"
    label 'process_low'

    publishDir "${params.outdir}/deliverables/metadata", mode: params.publish_dir_mode

    input:
    tuple val(n_contrast_manifests), val(n_global_manifests)
    val  outdir_root

    output:
    path "results_contract_report.txt", emit: report
    path "versions.yml",                emit: versions

    script:
    """
    # warnings (exit 1) must not abort the run; only structural failures (>= 2) do
    set +e
    python3 ${projectDir}/bin/validate_results_contract.py \\
        --outdir      ${outdir_root} \\
        --report-file results_contract_report.txt
    RC=\$?
    set -e
    if [ \$RC -ge 2 ]; then
        echo "[VALIDATE_RESULTS] results contract broken — see results_contract_report.txt" >&2
        exit \$RC
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}