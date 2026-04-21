/*
 * ========================================================================================
 *  ALTERNATIVE_SPLICING: Main workflow
 * ========================================================================================
 */

// Import subworkflows
include { INPUT_CHECK           } from '../subworkflows/local/input_check/main'
include { RMATS_ANALYSIS        } from '../subworkflows/local/rmats_analysis/main'
include { MAJIQ_ANALYSIS        } from '../subworkflows/local/majiq_analysis/main'
include { ISOFORMSWITCHR_ANALYSIS } from '../subworkflows/local/isoformswitchr_analysis/main'

// Import modules
include { RENDER_REPORT } from '../modules/local/report/main'

workflow ALTERNATIVE_SPLICING {
    
    // Create input channels
    ch_samplesheet = channel.fromPath(params.input, checkIfExists: true)
    ch_comparisons = channel.fromPath(params.comparisons, checkIfExists: true)
    ch_gtf         = channel.fromPath(params.gtf, checkIfExists: true).first()
    ch_report_rmd  = channel.fromPath("${projectDir}/bin/render_report.Rmd", checkIfExists: true).first()

    def no_rmats_dir = file("${workflow.projectDir}/assets/empty/NO_RMATS")
    def no_majiq_dir = file("${workflow.projectDir}/assets/empty/NO_MAJIQ")
    def no_isar_dir  = file("${workflow.projectDir}/assets/empty/NO_ISAR")
    
    /*
     * SUBWORKFLOW: Input validation and channel creation
     */
    INPUT_CHECK(
        ch_samplesheet,
        ch_comparisons,
        ch_gtf
    )
    
    // Channels from input_check
    ch_samples_bam    = INPUT_CHECK.out.samples_bam
    ch_samples_salmon = INPUT_CHECK.out.samples_salmon
    ch_comparisons_meta = INPUT_CHECK.out.comparisons

    // Split comparison metadata for independent consumers (DSL2: multiMap replaces into)
    ch_comparisons_meta
        .multiMap { meta ->
            rmats:    meta
            majiq:    meta
            isar:     meta
            for_ids:  meta
        }
        .set { ch_comp_split }

    // Comparison IDs are needed to assemble report inputs even when a branch is disabled
    ch_comp_split.for_ids
        .map { comp_meta -> comp_meta.id }
        .multiMap { comp_id ->
            rmats: comp_id
            majiq: comp_id
            isar:  comp_id
        }
        .set { ch_ids_split }
    
    /*
     * SUBWORKFLOW: rMATS-turbo analysis (optional)
     */
    ch_rmats_for_report = channel.empty()
    if (params.run_rmats) {
        RMATS_ANALYSIS(
            ch_samples_bam,
            ch_comp_split.rmats,
            ch_gtf
        )

        ch_ids_split.rmats
            .join(RMATS_ANALYSIS.out.results, by: 0)
            .map { comp_id, rmats_dir -> [comp_id, rmats_dir] }
            .set { ch_rmats_for_report }
    } else {
        ch_ids_split.rmats
            .map { comp_id -> [comp_id, no_rmats_dir] }
            .set { ch_rmats_for_report }
    }
    
    /*
     * SUBWORKFLOW: MAJIQ analysis (optional)
     */
    ch_majiq_for_report = channel.empty()
    if (params.run_majiq) {
        MAJIQ_ANALYSIS(
            ch_samples_bam,
            ch_comp_split.majiq,
            ch_gtf
        )

        ch_ids_split.majiq
            .join(MAJIQ_ANALYSIS.out.results, by: 0)
            .map { comp_id, majiq_dir -> [comp_id, majiq_dir] }
            .set { ch_majiq_for_report }
    } else {
        ch_ids_split.majiq
            .map { comp_id -> [comp_id, no_majiq_dir] }
            .set { ch_majiq_for_report }
    }
    
    /*
     * SUBWORKFLOW: IsoformSwitchAnalyzeR analysis (optional)
     */
    ch_isar_for_report = channel.empty()
    if (params.run_isar) {
        ISOFORMSWITCHR_ANALYSIS(
            ch_samples_salmon,
            ch_comp_split.isar,
            ch_gtf
        )

        ch_ids_split.isar
            .join(ISOFORMSWITCHR_ANALYSIS.out.results, by: 0)
            .map { comp_id, isar_dir -> [comp_id, isar_dir] }
            .set { ch_isar_for_report }
    } else {
        ch_ids_split.isar
            .map { comp_id -> [comp_id, no_isar_dir] }
            .set { ch_isar_for_report }
    }
    
    /*
     * MODULE: Render final R Markdown report per comparison
     */
    // Join all results for each comparison
    ch_rmats_for_report
        .join(ch_majiq_for_report, by: 0)
        .join(ch_isar_for_report, by: 0)
        .set { ch_report_inputs }
    
    RENDER_REPORT(
        ch_report_inputs.map { comparison_id, _rmats_dir, _majiq_dir, _isar_dir -> comparison_id },
        ch_report_inputs.map { _comparison_id, rmats_dir, _majiq_dir, _isar_dir -> rmats_dir },
        ch_report_inputs.map { _comparison_id, _rmats_dir, majiq_dir, _isar_dir -> majiq_dir },
        ch_report_inputs.map { _comparison_id, _rmats_dir, _majiq_dir, isar_dir -> isar_dir },
        ch_report_rmd                    // report Rmd template
    )
    
    /*
     * TODO: MODULE: MultiQC - aggregate QC reports
     * Future enhancement: collect logs from all processes
     */
    // ch_multiqc_files = channel.empty()
    // MULTIQC(ch_multiqc_files.collect())
}
