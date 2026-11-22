/*
 * ========================================================================================
 *  SASHIMI_MERGE: fold targeted sashimi plots into the automatic plot tree
 * ========================================================================================
 *  Combines the SASHIMI_PLOTS output (top-N) with the SASHIMI_TARGETED output
 *  (requested genes) into a single sashimi_out/ directory, so SASHIMI_INDEX and
 *  RENDER_REPORT keep their single-directory contract. Targeted PDFs and the
 *  targeted_classes.tsv manifest are copied next to the automatic ones.
 * ========================================================================================
 */

process SASHIMI_MERGE {
    tag "$comparison_id"
    label 'process_low'

    input:
    // Both inputs are directories named sashimi_out; stage them apart so they do
    // not collide in the task work dir.
    tuple val(comparison_id),
          path(auto_dir,     stageAs: 'auto_sashimi'),
          path(targeted_dir, stageAs: 'targeted_sashimi')

    output:
    tuple val(comparison_id), path("sashimi_out/"), emit: results
    path "versions.yml",                            emit: versions

    script:
    """
    cp -r auto_sashimi sashimi_out

    for ETYPE in SE A5SS A3SS MXE RI; do
        if compgen -G "targeted_sashimi/\${ETYPE}/Sashimi_plot/*.pdf" >/dev/null; then
            mkdir -p "sashimi_out/\${ETYPE}/Sashimi_plot"
            cp -L "targeted_sashimi/\${ETYPE}/Sashimi_plot/"*.pdf "sashimi_out/\${ETYPE}/Sashimi_plot/"
        elif compgen -G "targeted_sashimi/\${ETYPE}/*.pdf" >/dev/null; then
            mkdir -p "sashimi_out/\${ETYPE}"
            cp -L "targeted_sashimi/\${ETYPE}/"*.pdf "sashimi_out/\${ETYPE}/"
        fi
    done

    if [ -f "targeted_sashimi/targeted_classes.tsv" ]; then
        cp -L "targeted_sashimi/targeted_classes.tsv" sashimi_out/targeted_classes.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: "\$(cp --version 2>/dev/null | head -1 | grep -oEm1 '[0-9]+\\.[0-9]+' || echo unknown)"
    END_VERSIONS
    """
}
