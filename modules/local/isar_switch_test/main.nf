process ISAR_SWITCH_TEST {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'
    
    // Container resolved from modules.config (params.isar_container or ghcr.io default)
    
    input:
    tuple val(comparison_id), path(rds_input)
    
    output:
    tuple val(comparison_id), path("${comparison_id}_tested.rds"), emit: rds
    path "effective_test_method.txt"                             , emit: effective_method
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
        --iso_expr_cutoff ${params.isar_iso_expr_cutoff} \\
        --test_method ${params.isar_test_method}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        test_method: "\$(cat effective_test_method.txt)"
        requested_test_method: "${params.isar_test_method}"
        saturn: \$(Rscript -e "library(satuRn); cat(as.character(packageVersion('satuRn')))" 2>/dev/null || echo NA)
        dexseq: \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))" 2>/dev/null || echo NA)
    END_VERSIONS
    """
}
