process RMATS_POST {
    tag "$comparison_id"
    label 'process_high'
    label 'process_long'
    
    container 'docker.io/xinglab/rmats:v4.3.0'
    
    publishDir "${params.outdir}/rmats/${comparison_id}", mode: params.publish_dir_mode
    
    input:
    val  comparison_id
    path rmats_files_g1   // staged .rmats files for group 1
    path rmats_files_g2   // staged .rmats files for group 2
    path bams_g1          // staged BAM files for group 1 (symlinked, not copied)
    path bams_g2          // staged BAM files for group 2
    val  sample_ids_g1    // sample ids in the same order as bams_g1 -> order of b1.txt
    val  sample_ids_g2    // sample ids in the same order as bams_g2 -> order of b2.txt
    path gtf
    
    output:
    tuple val(comparison_id), path("${comparison_id}"), emit: results
    path "versions.yml"                               , emit: versions
    
    script:
    def rmats_lib_type = params.strandedness == 'unstranded' ? 'fr-unstranded' :
                         params.strandedness == 'forward'    ? 'fr-secondstrand' :
                         params.strandedness == 'reverse'    ? 'fr-firststrand' : 'fr-unstranded'
    def g1_bams_staged = (bams_g1 instanceof List ? bams_g1 : [bams_g1]).join(' ')
    def g2_bams_staged = (bams_g2 instanceof List ? bams_g2 : [bams_g2]).join(' ')
    // Materialize the column -> sample mapping next to the rMATS output so
    // downstream consumers (report, PEGASAS) can validate id order instead of
    // trusting a channel that could drift from the actual BAM list order.
    def g1_ids_staged = (sample_ids_g1 instanceof List ? sample_ids_g1 : [sample_ids_g1]).join("\n")
    def g2_ids_staged = (sample_ids_g2 instanceof List ? sample_ids_g2 : [sample_ids_g2]).join("\n")
    // novelSS options mirror what was used in PREP — controlled via params.rmats_novel_ss
    def novelss_opt = params.rmats_novel_ss ? '--novelSS' : ''
    def mil_opt     = params.rmats_novel_ss ? "--mil ${params.rmats_min_intron_length}" : ''
    def mel_opt     = params.rmats_novel_ss ? "--mel ${params.rmats_max_exon_length}" : ''
    """
    # Create BAM lists — resolve staged symlinks to absolute paths via realpath
    for bam in ${g1_bams_staged}; do realpath "\$bam"; done | paste -sd ',' - > b1.txt
    for bam in ${g2_bams_staged}; do realpath "\$bam"; done | paste -sd ',' - > b2.txt

    mkdir -p ${comparison_id}

    # Same order, one sample id per line — defines the column -> sample mapping.
    # Written inside the published output dir so downstream consumers (PEGASAS,
    # manual QC) can resolve it next to the rMATS counts tables.
    cat > ${comparison_id}/b1_samples.txt << 'EOF1'
${g1_ids_staged}
EOF1
    cat > ${comparison_id}/b2_samples.txt << 'EOF2'
${g2_ids_staged}
EOF2

    # Collect all staged .rmats files into merged_tmp
    mkdir -p merged_tmp
    for f in ${rmats_files_g1} ${rmats_files_g2}; do
        cp "\$f" merged_tmp/
    done

    # Run rMATS POST
    python /rmats/rmats.py \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --gtf ${gtf} \\
        -t ${params.rmats_read_type} \\
        --readLength ${params.read_length} \\
        --variable-read-length \\
        --allow-clipping \\
        --libType ${rmats_lib_type} \\
        --nthread ${task.cpus} \\
        --tstat ${params.rmats_tstat_threads} \\
        --cstat ${params.rmats_cstat} \\
        --individual-counts \\
        --od ${comparison_id} \\
        --tmp merged_tmp \\
        --task post \\
        ${novelss_opt} ${mil_opt} ${mel_opt}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(python /rmats/rmats.py --version 2>&1 | grep -oP 'v\\d+\\.\\d+\\.\\d+' || echo "4.3.0")
    END_VERSIONS
    """
}
