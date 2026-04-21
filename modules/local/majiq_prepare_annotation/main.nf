// Convert GTF to GFF3 format required by MAJIQ v3.
// The input GTF lacks gene records; gtf_to_gff3.py derives them from
// transcript coordinates and emits a gene -> transcript -> exon hierarchy.
process MAJIQ_PREPARE_ANNOTATION {
    tag "gtf_to_gff3"
    label 'process_low'

    container 'quay.io/biocontainers/pysam:0.22.1--py312hcfdcdd7_2'

    input:
    path gtf

    output:
    path "annotation.gff3", emit: gff3
    path "versions.yml"   , emit: versions

    script:
    """
    gtf_to_gff3.py ${gtf} annotation.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
