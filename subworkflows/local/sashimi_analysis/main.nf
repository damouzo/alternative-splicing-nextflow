/*
 * ========================================================================================
 *  SASHIMI_ANALYSIS subworkflow
 *  Generates sashimi plots for top significant rMATS events using rmats2sashimiplot.
 *
 *  Takes:
 *   - ch_rmats_results: [comp_id, rmats_dir]   (from RMATS_ANALYSIS)
 *   - ch_samples_bam:   [meta, bam, bai]        (from INPUT_CHECK)
 *   - ch_comparisons:   [comp_meta]             (from INPUT_CHECK)
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
    ch_comparisons    // [comp_meta]  — meta has .id, .group1, .group2

    main:
    // Build per-comparison BAM channels by matching samples to their comparison group.
    // ch_samples_bam carries meta with .condition; comp_meta has .group1 / .group2.
    // Strategy: combine, filter on condition match, then group by comparison.

    // Flatten to [comp_id, group_label, bam, bai]
    ch_comp_bam = ch_comparisons
        .combine(ch_samples_bam)
        .filter { comp_meta, sample_meta, bam, bai ->
            sample_meta.condition == comp_meta.group1 ||
            sample_meta.condition == comp_meta.group2
        }
        .map { comp_meta, sample_meta, bam, bai ->
            def which_group = (sample_meta.condition == comp_meta.group1) ? "b1" : "b2"
            [comp_meta.id, which_group, bam, bai]
        }

    // Collect all b1 BAMs + BAIs per comparison
    ch_b1 = ch_comp_bam
        .filter { comp_id, grp, bam, bai -> grp == "b1" }
        .map    { comp_id, grp, bam, bai -> [comp_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, bams, bais -> [comp_id, bams.flatten(), bais.flatten()] }

    // Collect all b2 BAMs + BAIs per comparison
    ch_b2 = ch_comp_bam
        .filter { comp_id, grp, bam, bai -> grp == "b2" }
        .map    { comp_id, grp, bam, bai -> [comp_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, bams, bais -> [comp_id, bams.flatten(), bais.flatten()] }

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
