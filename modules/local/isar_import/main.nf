process ISAR_IMPORT {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    container 'oras://community.wave.seqera.io/library/bioconductor-dexseq_bioconductor-isoformswitchanalyzer_bioconductor-saturn_bioconductor-tximeta:63f4f88fcb3af3a9'
    
    input:
    val comparison_id
    path samplesheet
    path gtf
    path transcript_fasta
    
    output:
    tuple val(comparison_id), path("${comparison_id}_raw.rds"), emit: rds
    path "versions.yml"                                        , emit: versions
    
    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    isar_import.R \\
        --samplesheet ${samplesheet} \\
        --gtf ${gtf} \\
        --transcript_fasta ${transcript_fasta} \\
        --output ${comparison_id}_raw.rds
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //; s/ .*//')
        isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """
}
