process SASHIMI_PLOTS {
    tag "$comparison_id"
    label 'process_medium'

    container 'docker.io/xinglab/rmats2sashimiplot:v2.0.4'

    publishDir "${params.outdir}/sashimi/${comparison_id}", mode: params.publish_dir_mode

    input:
    tuple val(comparison_id),
          path(rmats_dir,   stageAs: 'rmats_results'),
          path(b1_bams,     stageAs: 'b1_bams/*'),
          path(b1_bais,     stageAs: 'b1_bams/*'),
          path(b2_bams,     stageAs: 'b2_bams/*'),
          path(b2_bais,     stageAs: 'b2_bams/*')

    output:
    tuple val(comparison_id), path("sashimi_out/"), emit: results
    path  "versions.yml",                           emit: versions

    script:
    def top_n       = params.sashimi_top_n
    def fdr         = params.report_fdr_cutoff
    def dpsi        = params.report_dpsi_cutoff
    def label1      = params.sashimi_group1_label
    def label2      = params.sashimi_group2_label
    def exon_scale  = params.sashimi_exon_scale
    def intron_scale = params.sashimi_intron_scale

    """
    # Step 1: filter rMATS output to top-N events per event type
    python3 ${projectDir}/bin/filter_rmats_for_sashimi.py \\
        rmats_results/ \\
        filtered_events/ \\
        --top-n ${top_n} \\
        --fdr   ${fdr} \\
        --dpsi  ${dpsi}

    mkdir -p sashimi_out

    # Step 2: build comma-separated BAM lists for each group
    B1_BAMS=\$(ls b1_bams/*.bam 2>/dev/null | tr '\\n' ',' | sed 's/,\$//')
    B2_BAMS=\$(ls b2_bams/*.bam 2>/dev/null | tr '\\n' ',' | sed 's/,\$//')

    if [ -z "\$B1_BAMS" ] || [ -z "\$B2_BAMS" ]; then
        echo "[ERROR] No BAM files found in one or both groups" >&2
        exit 1
    fi

    # Step 3: run rmats2sashimiplot for each event type that has filtered events
    for ETYPE in SE A5SS A3SS MXE RI; do
        EFILE="filtered_events/\${ETYPE}.top.txt"
        [ -f "\$EFILE" ] || continue

        # Skip if file contains only a header line
        NEVENTS=\$(tail -n +2 "\$EFILE" | wc -l)
        [ "\$NEVENTS" -gt 0 ] || continue

        echo "[INFO] Plotting \$NEVENTS \$ETYPE events..."

        python3 /usr/local/bin/rmats2sashimiplot \\
            --b1 "\$B1_BAMS" \\
            --b2 "\$B2_BAMS" \\
            -t   "\$ETYPE" \\
            -e   "\$EFILE" \\
            --l1 "${label1}" \\
            --l2 "${label2}" \\
            --exon_s  ${exon_scale} \\
            --intron_s ${intron_scale} \\
            -o   "sashimi_out/\${ETYPE}/" \\
        || echo "[WARN] \$ETYPE sashimi plot failed — continuing"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats2sashimiplot: \$(python3 /usr/local/bin/rmats2sashimiplot --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
