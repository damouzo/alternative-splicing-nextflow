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

    """
    export MAJIQ_LICENSE_FILE="${params.majiq_license}"
    # Prevent OpenBLAS/OpenMP from spawning excessive threads on shared systems
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

    # PSI coverage per group
    majiq-v3 psi-coverage ${splicegraph} ${group1_name}.psicov ${sj1}
    majiq-v3 psi-coverage ${splicegraph} ${group2_name}.psicov ${sj2}

    # DeltaPSI — outputs both voila file and TSV directly
    majiq-v3 deltapsi \\
        --splicegraph ${splicegraph} \\
        -psi1 ${group1_name}.psicov \\
        -psi2 ${group2_name}.psicov \\
        --output-voila ${comparison_id}.dpsicov \\
        --output-tsv ${comparison_id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        majiq: \$(majiq-v3 --version 2>&1 | head -1 || echo "3.0")
    END_VERSIONS
    """
}
