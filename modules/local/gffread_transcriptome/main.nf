process GFFREAD_TRANSCRIPTOME {
    tag "transcript_extraction"
    label 'process_medium'
    
    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_3'
    
    input:
    path gtf
    path genome_fasta
    
    output:
    path "transcript_sequences.fa", emit: transcript_fasta
    path "versions.yml"           , emit: versions
    
    script:
    """
    gffread \\
        ${gtf} \\
        -g ${genome_fasta} \\
        -w transcript_sequences.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1 | sed 's/gffread //g')
    END_VERSIONS
    """
}
