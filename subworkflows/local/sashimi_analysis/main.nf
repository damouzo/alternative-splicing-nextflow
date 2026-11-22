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

include { SASHIMI_PLOTS }    from '../../../modules/local/sashimi_plots/main'
include { SASHIMI_TARGETED } from '../../../modules/local/sashimi_targeted/main'
include { SASHIMI_MERGE }    from '../../../modules/local/sashimi_merge/main'

workflow SASHIMI_ANALYSIS {

    take:
    ch_rmats_results  // [comp_id, rmats_dir]
    ch_samples_bam    // [meta, bam, bai]
    ch_comp_names     // [comp_id, group1_name, group2_name]

    main:
    // Build per-comparison BAM channels from validated metadata.
    // INPUT_CHECK already assigns comparison_id and group (1/2) to each sample.

    // Flatten to [comp_id, group_label, sample_id, bam, bai]
    ch_comp_bam = ch_samples_bam
        .map { sample_meta, bam, bai ->
            def which_group = sample_meta.group == 1 ? "b1" : "b2"
            [sample_meta.comparison_id, which_group, sample_meta.id, bam, bai]
        }

    // Collect all b1 BAMs + BAIs + sample IDs per comparison
    ch_b1 = ch_comp_bam
        .filter { comp_id, grp, sample_id, bam, bai -> grp == "b1" }
        .map    { comp_id, grp, sample_id, bam, bai -> [comp_id, sample_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, sample_ids, bams, bais ->
            def rows = [sample_ids, bams, bais].transpose()
            def uniq_rows = rows.unique { row -> row[1].toString() }
            // Sort ids/bams/bais together by BAM file name so the order matches
            // the alphabetical shell glob (b1_bams/*.bam) in SASHIMI_PLOTS;
            // otherwise --group-info tags replicate labels onto the wrong track.
            uniq_rows = uniq_rows.sort { row -> row[1].name }
            [
                comp_id,
                uniq_rows.collect { row -> row[1] },
                uniq_rows.collect { row -> row[2] },
                uniq_rows.collect { row -> row[0] }
            ]
        }

    // Collect all b2 BAMs + BAIs + sample IDs per comparison
    ch_b2 = ch_comp_bam
        .filter { comp_id, grp, sample_id, bam, bai -> grp == "b2" }
        .map    { comp_id, grp, sample_id, bam, bai -> [comp_id, sample_id, bam, bai] }
        .groupTuple(by: 0)
        .map    { comp_id, sample_ids, bams, bais ->
            def rows = [sample_ids, bams, bais].transpose()
            def uniq_rows = rows.unique { row -> row[1].toString() }
            // Sort ids/bams/bais together by BAM file name so the order matches
            // the alphabetical shell glob (b1_bams/*.bam) in SASHIMI_PLOTS;
            // otherwise --group-info tags replicate labels onto the wrong track.
            uniq_rows = uniq_rows.sort { row -> row[1].name }
            [
                comp_id,
                uniq_rows.collect { row -> row[1] },
                uniq_rows.collect { row -> row[2] },
                uniq_rows.collect { row -> row[0] }
            ]
        }

    // Join rMATS results with BAMs and real condition labels for each comparison
    ch_sashimi_input = ch_rmats_results
        .join(ch_b1, by: 0)
        .join(ch_b2, by: 0)
        .join(ch_comp_names, by: 0)
        .map { comp_id, rmats_dir, b1_bams, b1_bais, b1_ids, b2_bams, b2_bais, b2_ids, g1, g2 ->
            [comp_id, rmats_dir, b1_bams, b1_bais, b1_ids, b2_bams, b2_bais, b2_ids, g1, g2]
        }

    SASHIMI_PLOTS(ch_sashimi_input, file("${projectDir}/bin/filter_rmats_for_sashimi.py"))

    /*
     * Optional: targeted sashimi plots for a user-provided gene list.
     * A separate process keeps -resume selective: adding genes rebuilds only
     * the targeted plots (and the report), not the top-N automatic ones.
     * When the list is not set the process is skipped entirely.
     */
    if (params.sashimi_gene_must) {
        ch_gene_must = channel.fromPath(params.sashimi_gene_must, checkIfExists: true).first()

        SASHIMI_TARGETED(
            ch_sashimi_input,
            ch_gene_must,
            file("${projectDir}/bin/filter_rmats_for_sashimi.py"),
            file("${projectDir}/bin/finalize_targeted_sashimi.py")
        )

        SASHIMI_MERGE(
            SASHIMI_PLOTS.out.results
                .join(SASHIMI_TARGETED.out.results, by: 0)
        )

        ch_sashimi_results = SASHIMI_MERGE.out.results
        ch_sashimi_versions = SASHIMI_PLOTS.out.versions
            .mix(SASHIMI_TARGETED.out.versions)
            .mix(SASHIMI_MERGE.out.versions)
    } else {
        ch_sashimi_results  = SASHIMI_PLOTS.out.results
        ch_sashimi_versions = SASHIMI_PLOTS.out.versions
    }

    emit:
    results  = ch_sashimi_results
    versions = ch_sashimi_versions
}
