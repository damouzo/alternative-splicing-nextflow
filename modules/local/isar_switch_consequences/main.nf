process ISAR_SWITCH_CONSEQUENCES {
    tag "$comparison_id"
    label 'process_medium'
    label 'process_high_memory'
    
    // Container resolved from modules.config (params.isar_container or ghcr.io default)
    
    publishDir "${params.outdir}/isoformswitchr/${comparison_id}", mode: params.publish_dir_mode
    
    input:
    tuple val(comparison_id), path(rds_input)
    
    output:
    tuple val(comparison_id), path("${comparison_id}"), emit: results
    path "versions.yml"                               , emit: versions
    
    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    # Create output directory
    mkdir -p ${comparison_id}
    
    isar_switch_consequences.R \\
        --input ${rds_input} \\
        --output ${comparison_id}/${comparison_id}_final.rds \\
        --output_dir ${comparison_id} \\
        --alpha ${params.isar_alpha} \\
        --dif_cutoff ${params.isar_dif_cutoff} \\
        --top_n_plots 25
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """
}
