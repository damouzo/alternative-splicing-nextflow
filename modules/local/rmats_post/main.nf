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
    """
    # Create BAM lists — resolve staged symlinks to absolute paths via realpath
    for bam in ${g1_bams_staged}; do realpath "\$bam"; done | paste -sd ',' - > b1.txt
    for bam in ${g2_bams_staged}; do realpath "\$bam"; done | paste -sd ',' - > b2.txt

    # Collect all staged .rmats files into merged_tmp
    mkdir -p merged_tmp
    for f in ${rmats_files_g1} ${rmats_files_g2}; do
        cp "\$f" merged_tmp/
    done

    mkdir -p ${comparison_id}
    
    # Run rMATS POST
    python /rmats/rmats.py \\
        --b1 b1.txt \\
        --b2 b2.txt \\
        --gtf ${gtf} \\
        -t paired \\
        --readLength ${params.read_length} \\
        --variable-read-length \\
        --allow-clipping \\
        --libType ${rmats_lib_type} \\
        --nthread ${task.cpus} \\
        --tstat ${params.rmats_tstat_threads} \\
        --cstat ${params.rmats_cstat} \\
        --novelSS \\
        --mil ${params.rmats_min_intron_length} \\
        --mel ${params.rmats_max_exon_length} \\
        --individual-counts \\
        --od ${comparison_id} \\
        --tmp merged_tmp \\
        --task post
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(python /rmats/rmats.py --version 2>&1 | grep -oP 'v\\d+\\.\\d+\\.\\d+' || echo "4.3.0")
    END_VERSIONS
    """
}
