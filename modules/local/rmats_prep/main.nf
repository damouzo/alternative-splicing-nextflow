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
    
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    # Create single-sample BAM list for PREP
    echo "${bam}" > bam_list.txt
    
    # Create output and temp directories
    mkdir -p prep_output
    mkdir -p prep_tmp
    
    # Run rMATS PREP
    python /rmats/rmats.py \\
        --b1 bam_list.txt \\
        --gtf ${gtf} \\
        -t paired \\
        --readLength ${params.read_length} \\
        --variable-read-length \\
        --allow-clipping \\
        --libType ${rmats_lib_type} \\
        --nthread ${task.cpus} \\
        --od prep_output \\
        --tmp prep_tmp \\
        --task prep
    
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
