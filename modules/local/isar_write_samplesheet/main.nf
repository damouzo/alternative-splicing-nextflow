process ISAR_WRITE_SAMPLESHEET {
    tag "$comparison_id"
    label 'process_single'

    // exec: block — no shell script, no tool version to capture.
    // versions.yml is intentionally omitted for this file-writing helper.

    input:
    val comparison_id
    val sample_rows   // List<String> — one "sample_id,condition,replicate" per sample (no salmon_dir)

    output:
    tuple val(comparison_id), path("${comparison_id}_samplesheet.csv"), emit: samplesheet

    exec:
    // Write partial CSV (without salmon_dir) to the task work dir.
    // Salmon directories are passed as staged path inputs to ISAR_IMPORT so that
    // Nextflow mounts them correctly inside containers (Docker and Singularity).
    def rows = sample_rows instanceof List ? sample_rows : [sample_rows]
    task.workDir
        .resolve("${comparison_id}_samplesheet.csv")
        .text = "sample,condition,replicate\n" + rows.join('\n') + '\n'
}
