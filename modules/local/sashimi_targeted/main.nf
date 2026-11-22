/*
 * ========================================================================================
 *  SASHIMI_TARGETED: sashimi plots for user-requested genes (--sashimi_gene_must)
 * ========================================================================================
 *  Separate process from SASHIMI_PLOTS so that adding genes to the list and
 *  re-running with -resume only rebuilds the targeted plots (and the report),
 *  never the top-N automatic plots.
 *
 *  Plots every significant rMATS event (same FDR / |dPSI| cutoffs as the top-N)
 *  whose geneSymbol is in the list, for every event type. Output PDFs are named
 *  A_<GENE>_<EVENT_TYPE>_<coords>.pdf so two events of the same gene stay
 *  distinct.
 *
 *  No publishDir: the PDFs are staged into RENDER_REPORT, which is the single
 *  owner of the shippable plots/ tree.
 * ========================================================================================
 */

process SASHIMI_TARGETED {
    tag "$comparison_id"
    label 'process_medium'

    input:
    tuple val(comparison_id),
          path(rmats_dir,   stageAs: 'rmats_results'),
          path(b1_bams,     stageAs: 'b1_bams/*'),
          path(b1_bais,     stageAs: 'b1_bams/*'),
          val(b1_ids),
          path(b2_bams,     stageAs: 'b2_bams/*'),
          path(b2_bais,     stageAs: 'b2_bams/*'),
          val(b2_ids),
          val(group1_name),
          val(group2_name)
    path gene_must
    // Scripts as inputs so content edits invalidate the cache on -resume
    path filter_script
    path finalize_script

    output:
    tuple val(comparison_id), path("sashimi_out/"), emit: results
    path  "versions.yml",                          emit: versions

    script:
    def fdr          = params.report_fdr_cutoff
    def dpsi         = params.report_dpsi_cutoff
    def exon_scale   = params.sashimi_exon_scale
    def intron_scale = params.sashimi_intron_scale
    def g1_label     = group1_name ?: params.sashimi_group1_label
    def g2_label     = group2_name ?: params.sashimi_group2_label
    def b1_ids_lines = b1_ids.join('\n')
    def b2_ids_lines = b2_ids.join('\n')

    """
    # Step 1: select every significant event of the requested genes
    python3 ${filter_script} \\
        rmats_results/ \\
        targeted_events/ \\
        --fdr   ${fdr} \\
        --dpsi  ${dpsi} \\
        --gene-must ${gene_must}

    mkdir -p sashimi_out

    HAS_EVENTS=0
    for ETYPE in SE A5SS A3SS MXE RI; do
        if [ -f "targeted_events/\${ETYPE}.top.txt" ] \\
           && [ "\$(tail -n +2 "targeted_events/\${ETYPE}.top.txt" | wc -l)" -gt 0 ]; then
            HAS_EVENTS=1
        fi
    done

    if [ "\$HAS_EVENTS" -eq 0 ]; then
        echo "[INFO] No significant targeted events — recording requested genes in the manifest only"
        mkdir -p sashimi_out
        python3 ${finalize_script} \\
            --targeted-events targeted_events/ \\
            --sashimi-out sashimi_out/
        printf '%s\\n' \\
            '"${task.process}":' \\
            '    rmats2sashimiplot: "not run"' \\
            > versions.yml
        exit 0
    fi

    # Resolve rmats2sashimiplot executable across image layouts
    SASHIMI_BIN=\$(command -v rmats2sashimiplot || true)
    if [ -z "\$SASHIMI_BIN" ] && [ -x /rmats2sashimiplot/conda_env/bin/rmats2sashimiplot ]; then
        SASHIMI_BIN=/rmats2sashimiplot/conda_env/bin/rmats2sashimiplot
    fi
    if [ -z "\$SASHIMI_BIN" ]; then
        echo "[ERROR] rmats2sashimiplot executable not found in container" >&2
        exit 1
    fi
    if ! command -v samtools >/dev/null 2>&1; then
        echo "[ERROR] samtools is required by rmats2sashimiplot but is missing" >&2
        exit 1
    fi

    # Step 2: build comma-separated BAM lists for each group
    shopt -s nullglob
    B1_ARR=(b1_bams/*.bam)
    B2_ARR=(b2_bams/*.bam)
    if [ "\${#B1_ARR[@]}" -eq 0 ] || [ "\${#B2_ARR[@]}" -eq 0 ]; then
        echo "[ERROR] No BAM files found in one or both groups" >&2
        exit 1
    fi
    B1_BAMS="\$(printf '%s,' "\${B1_ARR[@]}")"; B1_BAMS="\${B1_BAMS%,}"
    B2_BAMS="\$(printf '%s,' "\${B2_ARR[@]}")"; B2_BAMS="\${B2_BAMS%,}"

    cat <<'B1_IDS_EOF' > b1_ids.txt
${b1_ids_lines}
B1_IDS_EOF
    cat <<'B2_IDS_EOF' > b2_ids.txt
${b2_ids_lines}
B2_IDS_EOF
    mapfile -t B1_IDS < b1_ids.txt
    mapfile -t B2_IDS < b2_ids.txt

    GROUP_FILE="sample_groups.gf"
    {
        for i in "\${!B1_IDS[@]}"; do
            printf '%s: %d\\n' "\${B1_IDS[i]}" "\$(( i + 1 ))"
        done
        for i in "\${!B2_IDS[@]}"; do
            printf '%s: %d\\n' "\${B2_IDS[i]}" "\$(( \${#B1_IDS[@]} + i + 1 ))"
        done
    } > "\$GROUP_FILE"

    # Step 3: run rmats2sashimiplot for each event type with targeted events
    for ETYPE in SE A5SS A3SS MXE RI; do
        EFILE="targeted_events/\${ETYPE}.top.txt"
        [ -f "\$EFILE" ] || continue
        NEVENTS=\$(tail -n +2 "\$EFILE" | wc -l)
        [ "\$NEVENTS" -gt 0 ] || continue

        echo "[INFO] Plotting \$NEVENTS targeted \$ETYPE events..."
        OUT_DIR="sashimi_out/\${ETYPE}"
        mkdir -p "\$OUT_DIR"
        if ! "\$SASHIMI_BIN" \\
            --b1 "\$B1_BAMS" --b2 "\$B2_BAMS" \\
            --event-type "\$ETYPE" \\
            -e   "\$EFILE" \\
            --l1 "${g1_label}" --l2 "${g2_label}" \\
            --exon_s  ${exon_scale} --intron_s ${intron_scale} \\
            --group-info "\$GROUP_FILE" \\
            -o   "\$OUT_DIR/"; then
            echo "[WARN] \$ETYPE targeted sashimi plot failed — continuing" >&2
        fi
    done

    # Step 4: rename to A_<GENE>_<TYPE>_<coords>.pdf and write the manifest
    python3 ${finalize_script} \\
        --targeted-events targeted_events/ \\
        --sashimi-out sashimi_out/

    SASHIMI_VERSION="\$("\$SASHIMI_BIN" -h 2>&1 | head -1 || echo unknown)"
    printf '%s\\n' \\
        '"${task.process}":' \\
        "    rmats2sashimiplot: \\"\$SASHIMI_VERSION\\"" \\
        > versions.yml
    """
}
