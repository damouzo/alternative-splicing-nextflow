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
        }
        .set { ch_comp_split }

    // Comparison IDs for the disabled-tool branches (emit [comp_id, placeholder_dir])
    ch_comparisons_meta
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
        // Output already carries [comparison_id, dir] — no join needed
        ch_rmats_for_report = RMATS_ANALYSIS.out.results
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
        ch_majiq_for_report = MAJIQ_ANALYSIS.out.results
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
        ch_isar_for_report = ISOFORMSWITCHR_ANALYSIS.out.results
    } else {
        ch_ids_split.isar
            .map { comp_id -> [comp_id, no_isar_dir] }
            .set { ch_isar_for_report }
    }
    
    /*
     * MODULE: Render final R Markdown report per comparison
     */
    // Join all results for each comparison into a single tuple.
    // Pass original dir names as vals so RENDER_REPORT can detect NO_* placeholders
    // even after stageAs renames the paths in the work directory.
    ch_rmats_for_report
        .join(ch_majiq_for_report, by: 0)
        .join(ch_isar_for_report, by: 0)
        .map { comp_id, rdir, mdir, idir ->
            [comp_id, rdir.name, rdir, mdir.name, mdir, idir.name, idir]
        }
        .set { ch_report_inputs }
    
    // Pass the full tuple — avoids multiple queue-channel consumers causing desync
    RENDER_REPORT(ch_report_inputs, ch_report_rmd)
    
    /*
     * TODO: MODULE: MultiQC - aggregate QC reports
     * Future enhancement: collect logs from all processes
     */
    // ch_multiqc_files = channel.empty()
    // MULTIQC(ch_multiqc_files.collect())
}
