process SASHIMI_PLOTS {
    tag "$comparison_id"
    label 'process_medium'

    publishDir "${params.outdir}/sashimi/${comparison_id}", mode: params.publish_dir_mode

    input:
    tuple val(comparison_id),
          path(rmats_dir,   stageAs: 'rmats_results'),
          path(b1_bams,     stageAs: 'b1_bams/*'),
          path(b1_bais,     stageAs: 'b1_bams/*'),
          val(b1_ids),
          path(b2_bams,     stageAs: 'b2_bams/*'),
          path(b2_bais,     stageAs: 'b2_bams/*'),
          val(b2_ids)

    output:
    tuple val(comparison_id), path("sashimi_out/"), emit: results
    path  "versions.yml",                           emit: versions

    script:
    def top_n       = params.sashimi_top_n
    def fdr         = params.report_fdr_cutoff
    def dpsi        = params.report_dpsi_cutoff
    def exon_scale  = params.sashimi_exon_scale
    def intron_scale = params.sashimi_intron_scale
    def b1_ids_quoted = b1_ids.collect { '"' + it.replace('"', '\\"') + '"' }.join(' ')
    def b2_ids_quoted = b2_ids.collect { '"' + it.replace('"', '\\"') + '"' }.join(' ')

    """
    # Step 1: filter rMATS output to top-N events per event type
    filter_rmats_for_sashimi.py \\
        rmats_results/ \\
        filtered_events/ \\
        --top-n ${top_n} \\
        --fdr   ${fdr} \\
        --dpsi  ${dpsi}

    mkdir -p sashimi_out

    # Resolve rmats2sashimiplot executable across image layouts
    SASHIMI_BIN=\$(command -v rmats2sashimiplot || true)
    if [ -z "\$SASHIMI_BIN" ] && [ -x /rmats2sashimiplot/conda_env/bin/rmats2sashimiplot ]; then
        SASHIMI_BIN=/rmats2sashimiplot/conda_env/bin/rmats2sashimiplot
    fi
    if [ -z "\$SASHIMI_BIN" ]; then
        echo "[ERROR] rmats2sashimiplot executable not found in container" >&2
        exit 1
    fi

    SASHIMI_HELP="\$("\$SASHIMI_BIN" -h 2>&1 || true)"
    if ! printf '%s\n' "\$SASHIMI_HELP" | grep -q -- '--event-type'; then
        echo "[ERROR] Incompatible rmats2sashimiplot CLI (missing --event-type)" >&2
        echo "[ERROR] Resolved executable: \$SASHIMI_BIN" >&2
        printf '%s\n' "\$SASHIMI_HELP" | head -20 >&2
        exit 1
    fi

    SASHIMI_VERSION_LINE="\$(printf '%s\n' "\$SASHIMI_HELP" | head -1 || echo unknown)"
    echo "[INFO] Using rmats2sashimiplot: \$SASHIMI_BIN"
    echo "[INFO] rmats2sashimiplot help: \$SASHIMI_VERSION_LINE"

    if ! command -v samtools >/dev/null 2>&1; then
        echo "[ERROR] samtools is required by rmats2sashimiplot but is missing from container" >&2
        echo "[ERROR] Configure --sashimi_container with an image that includes samtools" >&2
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

    for BAM in "\${B1_ARR[@]}" "\${B2_ARR[@]}"; do
        if [ ! -f "\${BAM}.bai" ] && [ ! -f "\${BAM%.bam}.bai" ]; then
            echo "[ERROR] Missing BAM index (.bai) for \$BAM" >&2
            exit 1
        fi
    done

    B1_BAMS="\$(printf '%s,' "\${B1_ARR[@]}")"
    B2_BAMS="\$(printf '%s,' "\${B2_ARR[@]}")"
    B1_BAMS="\${B1_BAMS%,}"
    B2_BAMS="\${B2_BAMS%,}"

    B1_IDS=(${b1_ids_quoted})
    B2_IDS=(${b2_ids_quoted})

    if [ "\${#B1_IDS[@]}" -ne "\${#B1_ARR[@]}" ] || [ "\${#B2_IDS[@]}" -ne "\${#B2_ARR[@]}" ]; then
        echo "[ERROR] Sample IDs and BAMs are inconsistent for sashimi plotting" >&2
        echo "[ERROR] B1 ids/bams: \${#B1_IDS[@]}/\${#B1_ARR[@]}" >&2
        echo "[ERROR] B2 ids/bams: \${#B2_IDS[@]}/\${#B2_ARR[@]}" >&2
        exit 1
    fi

    GROUP_FILE="sample_groups.gf"
    {
        for i in "\${!B1_IDS[@]}"; do
            printf '%s: %d\n' "\${B1_IDS[i]}" "$((i+1))"
        done
        for i in "\${!B2_IDS[@]}"; do
            printf '%s: %d\n' "\${B2_IDS[i]}" "$((\${#B1_IDS[@]} + i + 1))"
        done
    } > "\$GROUP_FILE"

    # Step 3: run rmats2sashimiplot for each event type that has filtered events
    HAD_EVENTS=0
    PLOT_OK=0
    for ETYPE in SE A5SS A3SS MXE RI; do
        EFILE="filtered_events/\${ETYPE}.top.txt"
        [ -f "\$EFILE" ] || continue

        # Skip if file contains only a header line
        NEVENTS=\$(tail -n +2 "\$EFILE" | wc -l)
        [ "\$NEVENTS" -gt 0 ] || continue
        HAD_EVENTS=1

        echo "[INFO] Plotting \$NEVENTS \$ETYPE events..."

        OUT_DIR="sashimi_out/\${ETYPE}"
        mkdir -p "\$OUT_DIR"

        if "\$SASHIMI_BIN" \\
            --b1 "\$B1_BAMS" \\
            --b2 "\$B2_BAMS" \\
            --event-type "\$ETYPE" \\
            -e   "\$EFILE" \\
            --l1 "group1" \\
            --l2 "group2" \\
            --exon_s  ${exon_scale} \\
            --intron_s ${intron_scale} \\
            --group-info "\$GROUP_FILE" \\
            -o   "\$OUT_DIR/"; then
            if compgen -G "\$OUT_DIR/*.pdf" >/dev/null || compgen -G "\$OUT_DIR/Sashimi_plot/*.pdf" >/dev/null; then
                PLOT_OK=1
            else
                echo "[WARN] \$ETYPE completed but no PDF plots were generated — continuing"
            fi
        else
            echo "[WARN] \$ETYPE sashimi plot failed — continuing"
        fi
    done

    if [ "\$HAD_EVENTS" -eq 1 ] && [ "\$PLOT_OK" -eq 0 ]; then
        echo "[ERROR] Failed to generate sashimi plots for all event types" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats2sashimiplot: "\${SASHIMI_VERSION_LINE:-unknown}"
    END_VERSIONS
    """
}
