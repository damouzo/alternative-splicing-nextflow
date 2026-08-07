/*
 * PEGASAS_PREPARE — shared gene expression matrix for the pathway (KS) step.
 *
 * The KS enrichment score of a sample is a function of that sample's expression
 * and the gene set only, so pathway scores are computed ONCE for the whole
 * cohort and subset per contrast downstream. This process builds the shared
 * gene expression matrix (samples x genes) and the global group info.
 */
process PEGASAS_PREPARE {
    tag "shared"
    label 'process_low'

    container 'local/pegasas:latest'

    input:
    path salmon_tpm
    path group_info

    output:
    path "pegasas_inputs/gene_exp_bySample.tsv", emit: gene_exp
    path "pegasas_inputs/group_info.tsv",        emit: group_info
    path "versions.yml",                          emit: versions

    script:
    """
    prepare_pegasas_gene_matrix.py \\
        ${salmon_tpm} \\
        ${group_info} \\
        --out-dir pegasas_inputs/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
