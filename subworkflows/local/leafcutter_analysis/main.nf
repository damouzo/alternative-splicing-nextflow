/*
 * ========================================================================================
 *  LEAFCUTTER_ANALYSIS: LeafCutter intron excision differential splicing
 * ========================================================================================
 *  Workflow:
 *    1. BAM2JUNC: Per-sample BAM → junction file
 *    2. CLUSTER: Cluster introns across all samples per comparison
 *    3. DS: Differential splicing test
 */

include { LEAFCUTTER_BAM2JUNC } from '../../../modules/local/leafcutter_bam2junc/main'
include { LEAFCUTTER_CLUSTER  } from '../../../modules/local/leafcutter_cluster/main'
include { LEAFCUTTER_DS       } from '../../../modules/local/leafcutter_ds/main'

workflow LEAFCUTTER_ANALYSIS {
    take:
    samples_bam    // channel: [meta, bam, bai] — meta has .comparison_id, .id, .condition
    gtf            // path: annotation GTF

    main:

    /*
     * Per-sample: extract splice junctions from BAM
     */
    LEAFCUTTER_BAM2JUNC(samples_bam)

    /*
     * Group junction files by comparison_id, carrying sample IDs, conditions
     * and junc files in a single groupTuple so their internal order is
     * guaranteed to stay aligned (two separate groupTuple calls are not).
     */
    LEAFCUTTER_BAM2JUNC.out.junc
        .map { meta, junc ->
            [meta.comparison_id, meta.id, meta.condition, junc]
        }
        .groupTuple(by: 0)
        .set { ch_leaf_grouped }  // [comp_id, [ids], [conds], [juncs]] — same order by construction

    /*
     * Per-comparison: cluster introns
     */
    LEAFCUTTER_CLUSTER(
        ch_leaf_grouped.map { comparison_id, sample_ids, _conditions, junc_files ->
            [comparison_id, sample_ids, junc_files]
        }
    )

    /*
     * Per-comparison: differential splicing
     * Join cluster counts with sample/condition metadata — ids and conditions
     * come from the same tuple as the juncs that fixed the count columns.
     */
    LEAFCUTTER_CLUSTER.out.counts
        .join(
            ch_leaf_grouped.map { comparison_id, sample_ids, conditions, _junc_files ->
                [comparison_id, sample_ids, conditions]
            }
        )
        .set { ch_ds_input }

    LEAFCUTTER_DS(ch_ds_input, gtf)

    emit:
    results  = LEAFCUTTER_DS.out.results
    versions = LEAFCUTTER_BAM2JUNC.out.versions
        .mix(LEAFCUTTER_CLUSTER.out.versions)
        .mix(LEAFCUTTER_DS.out.versions)
}
