process LEAFCUTTER_CLUSTER {
    tag "$comparison_id"
    label 'process_medium'

    input:
    tuple val(comparison_id), val(sample_ids), path(junc_files)

    output:
    tuple val(comparison_id), path("${comparison_id}_perind_numers.counts.gz"), emit: counts
    tuple val(comparison_id), path("${comparison_id}_perind.counts.gz")       , emit: perind
    path "versions.yml"                                                        , emit: versions

    script:
    """
    # Write the list of junction files
    printf '%s\n' ${junc_files.join(' ')} | tr ' ' '\n' > junc_file_list.txt

    # leafcutter_cluster_regtools.py in the pinned container expects
    # short/legacy options (-j/-o/-m).
    python3 /opt/leafcutter-src/clustering/leafcutter_cluster_regtools.py \
        -j junc_file_list.txt \
        -o ${comparison_id} \
        -m 1

    # Rename outputs to include comparison_id prefix
    [ -f "${comparison_id}_perind_numers.counts.gz" ] || \
        mv perind_numers.counts.gz ${comparison_id}_perind_numers.counts.gz 2>/dev/null || true
    [ -f "${comparison_id}_perind.counts.gz" ] || \
        mv perind.counts.gz ${comparison_id}_perind.counts.gz 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(python3 -c "import leafcutter; print(getattr(leafcutter, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
