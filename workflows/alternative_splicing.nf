/*
 * ========================================================================================
 *  ALTERNATIVE_SPLICING: Main workflow
 * ========================================================================================
 */

// Import subworkflows
include { INPUT_CHECK             } from '../subworkflows/local/input_check/main'
include { RMATS_ANALYSIS          } from '../subworkflows/local/rmats_analysis/main'
include { MAJIQ_ANALYSIS          } from '../subworkflows/local/majiq_analysis/main'
include { ISOFORMSWITCHR_ANALYSIS } from '../subworkflows/local/isoformswitchr_analysis/main'
include { SASHIMI_ANALYSIS        } from '../subworkflows/local/sashimi_analysis/main'
include { PEGASAS_ANALYSIS        } from '../subworkflows/local/pegasas_analysis/main'
include { LEAFCUTTER_ANALYSIS     } from '../subworkflows/local/leafcutter_analysis/main'

// Import modules
include { RENDER_REPORT } from '../modules/local/report/main'

workflow ALTERNATIVE_SPLICING {

    // Input channels
    ch_samplesheet = channel.fromPath(params.input,       checkIfExists: true)
    ch_comparisons = channel.fromPath(params.comparisons, checkIfExists: true)
    ch_gtf         = channel.fromPath(params.gtf,         checkIfExists: true).first()
    ch_report_rmd  = channel.fromPath("${projectDir}/bin/render_report.Rmd", checkIfExists: true).first()

    // Placeholder directories for disabled tools — sentinel names checked by name in RENDER_REPORT
    def no_rmats_dir      = file("${workflow.projectDir}/assets/empty/NO_RMATS")
    def no_majiq_dir      = file("${workflow.projectDir}/assets/empty/NO_MAJIQ")
    def no_isar_dir       = file("${workflow.projectDir}/assets/empty/NO_ISAR")
    def no_sashimi_dir    = file("${workflow.projectDir}/assets/empty/NO_SASHIMI")
    def no_pegasas_dir    = file("${workflow.projectDir}/assets/empty/NO_PEGASAS")
    def no_leafcutter_dir = file("${workflow.projectDir}/assets/empty/NO_LEAFCUTTER")

    /*
     * SUBWORKFLOW: Input validation and channel creation
     */
    INPUT_CHECK(
        ch_samplesheet,
        ch_comparisons,
        ch_gtf
    )

    ch_samples_bam      = INPUT_CHECK.out.samples_bam
    ch_samples_salmon   = INPUT_CHECK.out.samples_salmon
    ch_comparisons_meta = INPUT_CHECK.out.comparisons

    // Fan-out comparison IDs for disabled-tool fallback branches
    ch_comparisons_meta
        .map { meta -> meta.id }
        .multiMap { comp_id ->
            rmats:      comp_id
            majiq:      comp_id
            isar:       comp_id
            sashimi:    comp_id
            pegasas:    comp_id
            leafcutter: comp_id
        }
        .set { ch_ids_split }

    /*
     * SUBWORKFLOW: rMATS-turbo (optional)
     */
    ch_rmats_for_report   = channel.empty()
    ch_rmats_for_sashimi  = channel.empty()
    ch_rmats_for_pegasas  = channel.empty()

    if (params.run_rmats) {
        RMATS_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_rmats_for_report  = RMATS_ANALYSIS.out.results
        ch_rmats_for_sashimi = RMATS_ANALYSIS.out.results
        ch_rmats_for_pegasas = RMATS_ANALYSIS.out.results
    } else {
        ch_ids_split.rmats
            .map { comp_id -> [comp_id, no_rmats_dir] }
            .set { ch_rmats_for_report }
    }

    /*
     * SUBWORKFLOW: MAJIQ (optional)
     */
    ch_majiq_for_report = channel.empty()
    if (params.run_majiq) {
        MAJIQ_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_majiq_for_report = MAJIQ_ANALYSIS.out.results
    } else {
        ch_ids_split.majiq
            .map { comp_id -> [comp_id, no_majiq_dir] }
            .set { ch_majiq_for_report }
    }

    /*
     * SUBWORKFLOW: IsoformSwitchAnalyzeR (optional)
     */
    ch_isar_for_report = channel.empty()
    if (params.run_isar) {
        ISOFORMSWITCHR_ANALYSIS(
            ch_samples_salmon,
            ch_gtf
        )
        ch_isar_for_report = ISOFORMSWITCHR_ANALYSIS.out.results
    } else {
        ch_ids_split.isar
            .map { comp_id -> [comp_id, no_isar_dir] }
            .set { ch_isar_for_report }
    }

    /*
     * SUBWORKFLOW: Sashimi plots (optional; requires run_rmats = true)
     */
    ch_sashimi_for_report = channel.empty()
    if (params.run_sashimi && params.run_rmats) {
        SASHIMI_ANALYSIS(
            ch_rmats_for_sashimi,
            ch_samples_bam,
            ch_comparisons_meta
        )
        ch_sashimi_for_report = SASHIMI_ANALYSIS.out.results
    } else {
        ch_ids_split.sashimi
            .map { comp_id -> [comp_id, no_sashimi_dir] }
            .set { ch_sashimi_for_report }
    }

    /*
     * SUBWORKFLOW: PEGASAS pathway-splicing correlation (optional; requires run_rmats = true)
     */
    ch_pegasas_for_report = channel.empty()
    if (params.run_pegasas && params.run_rmats && params.salmon_merged_tpm && params.pegasas_groups) {
        ch_salmon_tpm = channel.fromPath(params.salmon_merged_tpm, checkIfExists: true).first()
        ch_gmt        = channel.fromPath(params.pathway_gmt ?: error(
            '--pathway_gmt is required when --run_pegasas is true'
        ), checkIfExists: true).first()

        PEGASAS_ANALYSIS(
            ch_rmats_for_pegasas,
            ch_salmon_tpm,
            ch_gmt,
            RMATS_ANALYSIS.out.sample_ids
        )
        ch_pegasas_for_report = PEGASAS_ANALYSIS.out.results
    } else {
        ch_ids_split.pegasas
            .map { comp_id -> [comp_id, no_pegasas_dir] }
            .set { ch_pegasas_for_report }
    }

    /*
     * SUBWORKFLOW: LeafCutter intron excision differential splicing (optional)
     */
    ch_leafcutter_for_report = channel.empty()
    if (params.run_leafcutter) {
        LEAFCUTTER_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_leafcutter_for_report = LEAFCUTTER_ANALYSIS.out.results
    } else {
        ch_ids_split.leafcutter
            .map { comp_id -> [comp_id, no_leafcutter_dir] }
            .set { ch_leafcutter_for_report }
    }

    /*
     * MODULE: Render per-comparison HTML report
     */
    ch_rmats_for_report
        .join(ch_majiq_for_report,      by: 0)
        .join(ch_isar_for_report,       by: 0)
        .join(ch_sashimi_for_report,    by: 0)
        .join(ch_pegasas_for_report,    by: 0)
        .join(ch_leafcutter_for_report, by: 0)
        .map { comp_id, rdir, mdir, idir, sdir, pdir, ldir ->
            [comp_id,
             rdir.name, rdir,
             mdir.name, mdir,
             idir.name, idir,
             sdir.name, sdir,
             pdir.name, pdir,
             ldir.name, ldir]
        }
        .set { ch_report_inputs }

    RENDER_REPORT(ch_report_inputs, ch_report_rmd)
}
