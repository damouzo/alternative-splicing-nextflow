process ISAR_SWITCH_TEST {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    // Container resolved from modules.config (params.isar_container or ghcr.io default)
    
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
        --dif_cutoff ${params.isar_dif_cutoff} \\
        --gene_expr_cutoff ${params.isar_gene_expr_cutoff} \\
        --iso_expr_cutoff ${params.isar_iso_expr_cutoff}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saturn: \$(Rscript -e "library(satuRn); cat(as.character(packageVersion('satuRn')))")
    END_VERSIONS
    """
}
