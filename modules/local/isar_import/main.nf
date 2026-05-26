process ISAR_IMPORT {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    // Container resolved from modules.config (params.isar_container or ghcr.io default)
    
    input:
    val  comparison_id
    path samplesheet    // partial CSV: sample,condition,replicate (no salmon_dir)
    path gtf
    path transcript_fasta
    path salmon_dirs    // staged Salmon output directories — one per sample, same order as CSV rows
    
    output:
    tuple val(comparison_id), path("${comparison_id}_raw.rds"), emit: rds
    path "versions.yml"                                        , emit: versions
    
    script:
    // Build space-separated list of staged salmon dir names for bash
    def staged_dirs = (salmon_dirs instanceof List ? salmon_dirs : [salmon_dirs])
                        .collect { it.name }.join(' ')
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    # Append staged salmon_dir paths to the partial samplesheet CSV.
    # paste merges the two files column-by-column; both have the same row order.
    staged_dirs=(${staged_dirs})
    printf '%s\\n' "\${staged_dirs[@]}" > salmon_dirs.txt
    paste -d',' ${samplesheet} salmon_dirs.txt > full_samplesheet.csv

    isar_import.R \\
        --samplesheet full_samplesheet.csv \\
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
