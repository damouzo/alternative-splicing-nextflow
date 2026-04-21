process ISAR_SWITCH_TEST {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    container 'oras://community.wave.seqera.io/library/bioconductor-dexseq_bioconductor-isoformswitchanalyzer_bioconductor-saturn_bioconductor-tximeta:63f4f88fcb3af3a9'
    
    input:
    tuple val(comparison_id), path(rds_input)
    
    output:
    tuple val(comparison_id), path("${comparison_id}_tested.rds"), emit: rds
    path "versions.yml"                                          , emit: versions
    
    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    isar_switch_test.R \\
        --input ${rds_input} \\
        --output ${comparison_id}_tested.rds \\
        --alpha ${params.isar_alpha} \\
        --dif_cutoff ${params.isar_dif_cutoff}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saturn: \$(Rscript -e "library(satuRn); cat(as.character(packageVersion('satuRn')))")
    END_VERSIONS
    """
}
