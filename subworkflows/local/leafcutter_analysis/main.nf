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
     * Group junction files by comparison_id, carry sample IDs and conditions
     * for groups.txt generation in LEAFCUTTER_DS
     */
    LEAFCUTTER_BAM2JUNC.out.junc
        .map { meta, junc ->
            [meta.comparison_id, meta.id, meta.condition, junc]
        }
        .groupTuple(by: 0)
        .map { comparison_id, sample_ids, _conditions, junc_files ->
            [comparison_id, sample_ids, junc_files]
        }
        .set { ch_grouped_juncs }

    // Carry conditions alongside for DS step
    LEAFCUTTER_BAM2JUNC.out.junc
        .map { meta, _junc ->
            [meta.comparison_id, meta.id, meta.condition]
        }
        .groupTuple(by: 0)
        .map { comparison_id, sample_ids, conditions ->
            [comparison_id, sample_ids, conditions]
        }
        .set { ch_grouped_meta }

    /*
     * Per-comparison: cluster introns
     */
    LEAFCUTTER_CLUSTER(ch_grouped_juncs)

    /*
     * Per-comparison: differential splicing
     * Join cluster counts with sample/condition metadata
     */
    LEAFCUTTER_CLUSTER.out.counts
        .join(ch_grouped_meta)
        .set { ch_ds_input }

    LEAFCUTTER_DS(ch_ds_input, gtf)

    emit:
    results  = LEAFCUTTER_DS.out.results
    versions = LEAFCUTTER_BAM2JUNC.out.versions
        .mix(LEAFCUTTER_CLUSTER.out.versions)
        .mix(LEAFCUTTER_DS.out.versions)
}
