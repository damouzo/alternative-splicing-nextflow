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
    val  sample_info  // [[sample_id, bam_path], ...]

    output:
    tuple val(comparison_id), path("sj/*.sj"), path("built_sg.zarr"), emit: majiq_build
    path "versions.yml"                                              , emit: versions

    script:
    // Build per-sample sj extraction commands (one per line)
    def sj_cmds = sample_info.collect { sample_id, bam_path ->
        "majiq-v3 sj ${bam_path} ann_sg.zarr sj/${sample_id}.sj"
    }.join('\n')

    // Build groups TSV content (all samples under a single group for the build step)
    def tsv_rows = sample_info.collect { sample_id, _ ->
        "all\t${sample_id}\tsj/${sample_id}.sj"
    }.join('\n')

    """
    export MAJIQ_LICENSE_FILE="${params.majiq_license}"
    # Prevent OpenBLAS/OpenMP from spawning excessive threads on shared systems
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

    # Step 1: Build zarr splicegraph from GFF3 annotation
    majiq-v3 gff3 ${gff3} ann_sg.zarr

    # Step 2: Extract splice junctions per sample (parallelisable in v3)
    mkdir -p sj
    ${sj_cmds}

    # Step 3: Build splicegraph across all samples
    printf 'group\\tprefix\\tsj\\n${tsv_rows}\\n' > build_config.tsv
    majiq-v3 build ann_sg.zarr built_sg.zarr \\
        --groups-tsv build_config.tsv \\
        -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        majiq: \$(majiq-v3 --version 2>&1 | head -1 || echo "3.0")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
