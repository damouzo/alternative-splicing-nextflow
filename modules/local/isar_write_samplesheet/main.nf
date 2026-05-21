process ISAR_WRITE_SAMPLESHEET {
    tag "$comparison_id"
    label 'process_single'

    input:
    val comparison_id
    val sample_rows   // List<String> — one "sample_id,condition,replicate,salmon_dir" per sample

    output:
    tuple val(comparison_id), path("${comparison_id}_samplesheet.csv"), emit: samplesheet

    exec:
    // Write CSV to the task's work directory — tracked by Nextflow, never touches the launch dir
    def rows = sample_rows instanceof List ? sample_rows : [sample_rows]
    task.workDir
        .resolve("${comparison_id}_samplesheet.csv")
        .text = "sample,condition,replicate,salmon_dir\n" + rows.join('\n') + '\n'
}
