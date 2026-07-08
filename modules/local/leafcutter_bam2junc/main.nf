process LEAFCUTTER_BAM2JUNC {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.junc"), emit: junc
    path "versions.yml"                     , emit: versions

    script:
    // Map strandedness to regtools -s flag:
    //   XS=unstranded, RF=first-strand (nf-core 'reverse'),
    //   FR=second-strand (nf-core 'forward').
    // Source: regtools src/junctions/junctions_extractor.{h,cc} (RF/FR first/second-strand).
    // Keep in sync with --libType mapping in modules/local/rmats_prep/main.nf.
    def strand_flag = params.strandedness == 'forward'  ? 'FR' :
                      params.strandedness == 'reverse'  ? 'RF' : 'XS'
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
