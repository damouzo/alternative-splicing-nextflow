process MAJIQ_BUILD {
    tag "$comparison_id"
    label 'process_high'
    label 'process_high_memory'

    container 'local/majiq:3.0'

    input:
    val  comparison_id
    path bams
    path bais
    path gff3         // annotation.gff3 — converted from GTF by MAJIQ_PREPARE_ANNOTATION
    val  sample_info  // [[sample_id], ...] — order must match staged bams list

    output:
    tuple val(comparison_id), path("sj/*.sj"), path("built_sg.zarr"), emit: majiq_build
    path "versions.yml"                                              , emit: versions

    script:
    // Derive sample IDs in order — must match the order of staged bams
    def sample_ids     = sample_info.collect { it[0] }
    def sample_ids_str = sample_ids.join(' ')

    // Build groups TSV: only 'group' and 'sj' columns are required by majiq-v3 build
    def tsv_rows = sample_ids.collect { sid ->
        "all\tsj/${sid}.sj"
    }.join('\n')

    """
    # Set license only when a non-null, non-empty path is provided
    ${(params.majiq_license && params.majiq_license != 'null') ? "export MAJIQ_LICENSE_FILE=\"${params.majiq_license}\"" : "# MAJIQ_LICENSE_FILE not set — rely on environment"}
    # Prevent OpenBLAS/OpenMP from spawning excessive threads on shared systems
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

    # Step 1: Build zarr splicegraph from GFF3 annotation
    majiq-v3 gff3 ${gff3} ann_sg.zarr

    # Step 2: Extract splice junctions per sample using staged BAM paths.
    # bams are staged in the work dir — iterate using the declared sample order.
    mkdir -p sj
    bam_array=(${(bams instanceof List ? bams : [bams]).join(' ')})
    sid_array=(${sample_ids_str})
    for i in "\${!bam_array[@]}"; do
        majiq-v3 sj "\${bam_array[\$i]}" ann_sg.zarr "sj/\${sid_array[\$i]}.sj"
    done

    # Step 3: Build splicegraph across all samples
    printf 'group\\tsj\\n${tsv_rows}\\n' > build_config.tsv
    majiq-v3 build ann_sg.zarr built_sg.zarr \\
        --groups-tsv build_config.tsv \\
        --minreads ${params.majiq_min_reads} \\
        --minpos ${params.majiq_min_intronic_cov} \\
        --mindenovo ${params.majiq_min_denovo_reads} \\
        -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        majiq: \$(majiq-v3 --version 2>&1 | head -1 || echo "3.0")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
