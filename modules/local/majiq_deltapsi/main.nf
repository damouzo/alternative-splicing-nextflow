process MAJIQ_DELTAPSI {
    tag "$comparison_id"
    label 'process_high'

    container 'local/majiq:3.0'

    input:
    val  comparison_id
    path sj_files_g1    // .sj files for group 1
    path sj_files_g2    // .sj files for group 2
    path splicegraph    // built_sg.zarr directory
    val  group1_name
    val  group2_name

    output:
    tuple val(comparison_id), path("${comparison_id}.dpsicov"), path("${comparison_id}.tsv"), emit: deltapsi
    path "versions.yml"                                                                      , emit: versions

    script:
    def sj1 = sj_files_g1 instanceof List ? sj_files_g1.join(' ') : sj_files_g1
    def sj2 = sj_files_g2 instanceof List ? sj_files_g2.join(' ') : sj_files_g2
    // Guard against identical group names (would silently overwrite the first psicov file)
    def g1_label = group1_name
    def g2_label = group1_name == group2_name ? "${group2_name}_g2" : group2_name

    """
    # Set license only when a non-null, non-empty path is provided
    ${(params.majiq_license && params.majiq_license != 'null') ? "export MAJIQ_LICENSE_FILE=\"${params.majiq_license}\"" : "# MAJIQ_LICENSE_FILE not set — rely on environment"}
    # Prevent OpenBLAS/OpenMP from spawning excessive threads on shared systems
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

    # PSI coverage per group — apply read/bin thresholds at this step
    majiq-v3 psi-coverage ${splicegraph} ${g1_label}.psicov \
        --minreads ${params.majiq_min_reads} \
        --minbins ${params.majiq_min_nonzero} \
        ${sj1}
    majiq-v3 psi-coverage ${splicegraph} ${g2_label}.psicov \
        --minreads ${params.majiq_min_reads} \
        --minbins ${params.majiq_min_nonzero} \
        ${sj2}

    # DeltaPSI — outputs both voila file and TSV directly
    majiq-v3 deltapsi \\
        --splicegraph ${splicegraph} \\
        -psi1 ${g1_label}.psicov \\
        -psi2 ${g2_label}.psicov \\
        --output-voila ${comparison_id}.dpsicov \\
        --output-tsv ${comparison_id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        majiq: \$(majiq-v3 --version 2>&1 | head -1 || echo "3.0")
    END_VERSIONS
    """
}
