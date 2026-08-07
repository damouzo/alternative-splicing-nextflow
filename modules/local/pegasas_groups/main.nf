process PEGASAS_GROUPS {
    tag "pegasas_groups"
    label 'process_low'

    input:
    path samplesheet

    output:
    path "pegasas_groups.tsv", emit: groups

    script:
    """
    tail -n +2 "${samplesheet}" | cut -d',' -f1,2 | tr ',' '\t' > pegasas_groups.tsv
    """
}
