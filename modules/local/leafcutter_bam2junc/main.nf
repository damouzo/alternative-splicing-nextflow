process LEAFCUTTER_BAM2JUNC {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.junc"), emit: junc
    path "versions.yml"                     , emit: versions

    script:
    // Map strandedness to regtools -s flag: 0=unstranded, 1=first-strand, 2=second-strand
    def strand_flag = params.strandedness == 'forward'  ? '1' :
                      params.strandedness == 'reverse'  ? '2' : '0'
    """
    regtools junctions extract \\
        -s ${strand_flag} \\
        -a 8 \\
        -m 50 \\
        -M 500000 \\
        ${bam} \\
        -o ${meta.id}.junc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(regtools --version 2>&1 | head -1 | sed 's/regtools //')
    END_VERSIONS
    """
}
