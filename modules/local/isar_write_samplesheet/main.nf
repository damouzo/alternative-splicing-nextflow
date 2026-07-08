process ISAR_WRITE_SAMPLESHEET {
    tag "$comparison_id"
    label 'process_single'

    // versions.yml is intentionally omitted for this file-writing helper.

    input:
    val comparison_id
    val sample_rows   // List<String> — one "sample_id,condition,replicate" per sample (no salmon_dir)

    output:
    tuple val(comparison_id), path("${comparison_id}_samplesheet.csv"), emit: samplesheet

    script:
    // Write partial CSV (without salmon_dir) in the task sandbox.
    // Salmon directories are passed as staged path inputs to ISAR_IMPORT so that
    // Nextflow mounts them correctly inside containers (Docker and Singularity).
    def rows = sample_rows instanceof List ? sample_rows : [sample_rows]
    def csv_lines = (["sample,condition,replicate"] + rows).join('\n')
    def csv_escaped = csv_lines.replace("'", "'\"'\"'")
    """
    printf '%s\n' '${csv_escaped}' > "${comparison_id}_samplesheet.csv"
    """
}
