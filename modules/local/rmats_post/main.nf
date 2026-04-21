process RMATS_POST {
    tag "$comparison_id"
    label 'process_high'
    label 'process_long'
    
    container 'docker.io/xinglab/rmats:v4.3.0'
    
    publishDir "${params.outdir}/rmats/${comparison_id}", mode: params.publish_dir_mode
    
    input:
    val comparison_id
    val rmats_files_g1   // list of absolute path strings — not staged
    val rmats_files_g2
    val bams_g1
    val bams_g2
    path gtf
    
    output:
    tuple val(comparison_id), path("${comparison_id}"), emit: results
    path "versions.yml"                                , emit: versions
    
    script:
    def strandedness = params.strandedness
    def rmats_lib_type = strandedness == 'unstranded' ? 'fr-unstranded' :
                         strandedness == 'forward'    ? 'fr-secondstrand' :
                         strandedness == 'reverse'    ? 'fr-firststrand' : 'fr-unstranded'
    
    """
    # Create BAM lists for group1 and group2
    echo "${bams_g1.join(',')}" > b1.txt
    echo "${bams_g2.join(',')}" > b2.txt

    # Create merged tmp directory with all .rmats files
    mkdir -p merged_tmp

    # Copy all group1 .rmats files (absolute paths from PREP work dirs)
    for f in ${rmats_files_g1.join(' ')}; do
        cp "\$f" merged_tmp/
    done

    # Copy all group2 .rmats files
    for f in ${rmats_files_g2.join(' ')}; do
        cp "\$f" merged_tmp/
    done
    
    # Create output directory
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
