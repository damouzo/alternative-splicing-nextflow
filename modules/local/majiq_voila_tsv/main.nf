// In MAJIQ v3, the deltapsi TSV is generated directly by majiq-v3 deltapsi.
// This module organises the outputs into a results directory for the report.
process MAJIQ_VOILA_TSV {
    tag "$comparison_id"
    label 'process_low'

    container 'local/majiq:3.0'

    publishDir "${params.outdir}/majiq/${comparison_id}", mode: params.publish_dir_mode

    input:
    val  comparison_id
    path splicegraph     // built_sg.zarr directory
    path dpsicov         // .dpsicov voila file
    path tsv             // deltapsi TSV

    output:
    tuple val(comparison_id), path("${comparison_id}"), emit: results
    path "versions.yml"                               , emit: versions

    script:
    """
    mkdir -p ${comparison_id}
    cp ${tsv} ${comparison_id}/
    cp -r ${dpsicov}    ${comparison_id}/
    cp -r ${splicegraph} ${comparison_id}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        majiq: \$(majiq-v3 --version 2>&1 | head -1 || echo "3.0")
    END_VERSIONS
    """
}
