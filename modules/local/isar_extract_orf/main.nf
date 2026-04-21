process ISAR_EXTRACT_ORF {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    container 'oras://community.wave.seqera.io/library/bioconductor-dexseq_bioconductor-isoformswitchanalyzer_bioconductor-saturn_bioconductor-tximeta:63f4f88fcb3af3a9'
    
    input:
    tuple val(comparison_id), path(rds_input)
    path gtf
    
    output:
    tuple val(comparison_id), path("${comparison_id}_orf.rds"), emit: rds
    path "isoformSwitchAnalyzeR_*.fasta"                      , emit: fastas
    path "versions.yml"                                       , emit: versions
    
    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    isar_extract_orf.R \\
        --input ${rds_input} \\
        --gtf ${gtf} \\
        --output ${comparison_id}_orf.rds \\
        --output_dir .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """
}
