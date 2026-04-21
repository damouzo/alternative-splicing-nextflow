process VALIDATE_INPUT {
    tag "samplesheet"
    label 'process_single'
    
    container 'quay.io/biocontainers/pysam:0.22.1--py312hcfdcdd7_2'
    
    input:
    path samplesheet
    path gtf
    val read_length
    
    output:
    path "samplesheet_validated.csv", emit: csv
    path "versions.yml"             , emit: versions
    
    script:
    """
    validate_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --gtf ${gtf} \\
        --read_length ${read_length} \\
        --output samplesheet_validated.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
