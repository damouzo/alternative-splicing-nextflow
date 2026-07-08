/*
 * ========================================================================================
 *  SASHIMI_ANALYSIS subworkflow
 *  Generates sashimi plots for top significant rMATS events using rmats2sashimiplot.
 *
 *  Takes:
 *   - ch_rmats_results: [comp_id, rmats_dir]   (from RMATS_ANALYSIS)
 *   - ch_samples_bam:   [meta, bam, bai]        (from INPUT_CHECK)
 *
 *  Emits:
 *   - results: [comp_id, sashimi_dir]
 * ========================================================================================
 */

include { SASHIMI_PLOTS } from '../../../modules/local/sashimi_plots/main'

workflow SASHIMI_ANALYSIS {

    take:
    ch_rmats_results  // [comp_id, rmats_dir]
    ch_samples_bam    // [meta, bam, bai]

    main:
    // Build per-comparison BAM channels from validated metadata.
    // INPUT_CHECK already assigns comparison_id and group (1/2) to each sample.

    // Flatten to [comp_id, group_label, bam, bai]
    ch_comp_bam = ch_samples_bam
        .map { sample_meta, bam, bai ->
            def which_group = sample_meta.group == 1 ? "b1" : "b2"
            [sample_meta.comparison_id, which_group, bam, bai]
        }

    // Collect all b1 BAMs + BAIs per comparison
    ch_b1 = ch_comp_bam
        .filter { comp_id, grp, bam, bai -> grp == "b1" }
        .map    { comp_id, grp, bam, bai -> [comp_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, bams, bais ->
            def uniq_bams = bams.flatten().unique { it.toString() }
            def uniq_bais = bais.flatten().unique { it.toString() }
            [comp_id, uniq_bams, uniq_bais]
        }

    // Collect all b2 BAMs + BAIs per comparison
    ch_b2 = ch_comp_bam
        .filter { comp_id, grp, bam, bai -> grp == "b2" }
        .map    { comp_id, grp, bam, bai -> [comp_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, bams, bais ->
            def uniq_bams = bams.flatten().unique { it.toString() }
            def uniq_bais = bais.flatten().unique { it.toString() }
            [comp_id, uniq_bams, uniq_bais]
        }

    // Join rMATS results with BAMs for each comparison
    ch_sashimi_input = ch_rmats_results
        .join(ch_b1, by: 0)
        .join(ch_b2, by: 0)
        .map { comp_id, rmats_dir, b1_bams, b1_bais, b2_bams, b2_bais ->
            [comp_id, rmats_dir, b1_bams, b1_bais, b2_bams, b2_bais]
        }

    SASHIMI_PLOTS(ch_sashimi_input)

    emit:
    results = SASHIMI_PLOTS.out.results
}
