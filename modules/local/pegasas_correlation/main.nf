/*
 * PEGASAS_CORRELATION — per-contrast pathway x splicing correlation.
 *
 * Reuses the SHARED pathway scores (pathway_out/, all samples) and subsets them
 * to the contrast via the per-contrast group_order.txt (2 groups). The two
 * upstream Python steps are invoked through the installed PEGASAS package
 * (resolved via PEGASAS.config resource paths); the final R correlation step
 * uses the project-bin cor_matrix_direct_perm.R (fixed version: lazy PDF, all
 * events in the global matrix). This avoids rebuilding the container.
 *
 * Output: correlation_out/ — one subdir per pathway with:
 *   <pathway>_global_cor_matrix.txt  (ALL events: r, p, significant)
 *   <pathway>_high_cor_matrix.txt    (significant events only)
 *   <pathway>_high_cor_scatterplots.pdf (only if >=1 significant event)
 */
process PEGASAS_CORRELATION {
    tag "$comparison_id"
    label 'process_medium'

    container 'local/pegasas:latest'

    input:
    tuple val(comparison_id),
          path(pathway_out),
          path(psi_matrix),
          path(group_order),
          path(contrast_samples)
    // Script as input so content edits invalidate the cache on -resume
    path cor_script

    output:
    tuple val(comparison_id), path("correlation_out/"), emit: results
    path "versions.yml",                            emit: versions

    script:
    // Export the allocated CPUs so the vendored R script can parallelise the
    // permutation test with mclapply (apptainer does not propagate NSLOTS).
    def ncores = Math.min(task.cpus as int, 8)
    """
    mkdir -p correlation_out
    export NSLOTS=${ncores}

    # Resolve the packaged upstream Python helpers once (unchanged PEGASAS scripts).
    REORDER_PY=\$(python3 -c "import PEGASAS.config as c; print(c.MAT_REORDER)")
    GEN_PY=\$(python3 -c "import PEGASAS.config as c; print(c.MAT_GENERATE)")

    # -L: pathway_out is staged as a symlink to the shared pathway run; find
    # must follow it. Scores live at pathway_out/<pathway>/<pathway>.scores.txt.
    for SCORES in \$(find -L ${pathway_out} -name '*.scores.txt' | sort); do
        SIG=\$(basename "\$SCORES" .scores.txt)
        mkdir -p "correlation_out/\$SIG"

        # Step 1: reorder signature scores by group (contrast group_order subsets
        # the shared scores to the contrast samples).
        python3 "\$REORDER_PY" "\$SCORES" "correlation_out" "${group_order}"

        # Step 2: align the PSI matrix to the same sample order.
        python3 "\$GEN_PY" "${psi_matrix}" "correlation_out/\$SIG/\$SIG.sorted.txt"

        # Step 3: Pearson correlation + permutation test (vendored fixed R script).
        REFINED=\$(ls "correlation_out/\$SIG"/refinedBySample.*.txt 2>/dev/null | head -1)
        if [ -n "\$REFINED" ]; then
            Rscript ${cor_script} \\
                "correlation_out/\$SIG/\$SIG.sorted.txt" \\
                "\$REFINED" \\
                "correlation_out/\$SIG" \\
                ${params.pegasas_perms} \\
                || echo "[WARN] correlation failed for \$SIG — continuing"
        else
            echo "[WARN] no refined PSI matrix for \$SIG — skipping"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PEGASAS: \$(PEGASAS --version 2>&1 | head -1)
    END_VERSIONS
    """
}
