process RMATS_PREP {
    tag "$meta.id"
    label 'process_medium'
    
    container 'docker.io/xinglab/rmats:v4.3.0'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    
    output:
    tuple val(meta), path("${meta.id}_*.rmats"), emit: rmats_files
    path "versions.yml"                       , emit: versions
    
    script:
    def strandedness = params.strandedness
    def rmats_lib_type = strandedness == 'unstranded' ? 'fr-unstranded' :
                         strandedness == 'forward'    ? 'fr-secondstrand' :
                         strandedness == 'reverse'    ? 'fr-firststrand' : 'fr-unstranded'
    // novelSS options must be identical in PREP and POST — controlled via params.rmats_novel_ss
    def novelss_opt = params.rmats_novel_ss ? '--novelSS' : ''
    def mil_opt     = params.rmats_novel_ss ? "--mil ${params.rmats_min_intron_length}" : ''
    def mel_opt     = params.rmats_novel_ss ? "--mel ${params.rmats_max_exon_length}" : ''
    
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    # Create single-sample BAM list — use realpath so the path is valid
    # regardless of stageInMode (copy/symlink/rellink)
    echo "\$(realpath "${bam}")" > bam_list.txt
    
    # Create output and temp directories
    mkdir -p prep_output
    mkdir -p prep_tmp
    
    # Run rMATS PREP
    python /rmats/rmats.py \\
        --b1 bam_list.txt \\
        --gtf ${gtf} \\
        -t ${params.rmats_read_type} \\
        --readLength ${params.read_length} \\
        --variable-read-length \\
        --allow-clipping \\
        --libType ${rmats_lib_type} \\
        --nthread ${task.cpus} \\
        --od prep_output \\
        --tmp prep_tmp \\
        --task prep \\
        ${novelss_opt} ${mil_opt} ${mel_opt}
    
    # Copy .rmats files to CWD with sample prefix (avoids staging subdirectory collision in POST)
    for f in prep_tmp/*.rmats; do
        if [ -f "\$f" ]; then
            basename=\$(basename "\$f")
            cp "\$f" "${meta.id}_\${basename}"
        fi
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(python /rmats/rmats.py --version 2>&1 | grep -oP 'v\\d+\\.\\d+\\.\\d+' || echo "4.3.0")
    END_VERSIONS
    """
}
