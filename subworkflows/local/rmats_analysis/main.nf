/*
 * ========================================================================================
 *  RMATS_ANALYSIS: rMATS-turbo differential alternative splicing analysis
 * ========================================================================================
 *  Workflow:
 *    1. PREP: Run per-sample preprocessing in parallel
 *    2. POST: Aggregate all samples per comparison and run statistical testing
 */

include { RMATS_PREP } from '../../../modules/local/rmats_prep/main'
include { RMATS_POST } from '../../../modules/local/rmats_post/main'

workflow RMATS_ANALYSIS {
    take:
    samples_bam  // channel: [meta, bam, bai] with meta.comparison_id and meta.group
    gtf          // path: annotation.gtf
    
    main:
    
    /*
     * Run RMATS PREP for each sample
     */
    RMATS_PREP(
        samples_bam,
        gtf
    )
    
    /*
     * Group .rmats files by comparison and group
     */
    RMATS_PREP.out.rmats_files
        .map { meta, rmats_files ->
            [meta.comparison_id, meta.group, meta, rmats_files]
        }
        .set { ch_rmats_prep_grouped }
    
    // Also group BAM files for POST
    samples_bam
        .map { meta, bam, _bai ->
            [meta.comparison_id, meta.group, meta.id, bam]
        }
        .set { ch_bams_grouped }
    
    /*
     * Collect all files per comparison
     * Group by comparison_id, then separate group1 and group2
     */
    ch_rmats_prep_grouped
        .groupTuple(by: 0)  // Group by comparison_id
        .map { comparison_id, groups, metas, rmats_files_list ->
            // Separate group 1 and group 2, tracking samples that did not
            // produce .rmats files (RMATS_PREP output is optional: true
            // to avoid Nextflow aborting when an empty BAM slips through).
            def g1_files = []
            def g2_files = []
            def missing  = []

            groups.eachWithIndex { group, idx ->
                def files = rmats_files_list[idx]
                if (!files) {
                    missing.add(metas[idx].id)
                    return
                }
                // rmats_files_list[idx] may be a single Path or a List<Path>
                // depending on how many files the glob matched — normalise to list
                def fileList = files instanceof List ? files : [files]
                if (group == 1) {
                    g1_files.addAll(fileList)
                } else {
                    g2_files.addAll(fileList)
                }
            }

            if (missing) {
                log.warn "[RMATS_ANALYSIS] ${comparison_id}: no .rmats files for samples ${missing} — excluding from POST"
            }
            if (!g1_files || !g2_files) {
                error "[RMATS_ANALYSIS] ${comparison_id}: missing .rmats files for group 1 or group 2 — cannot run rMATS POST. Check that the affected BAMs have reads aligned to splice junctions and match the GTF annotation."
            }

            // Keep as Path objects — Nextflow will stage them properly in RMATS_POST work dir
            [comparison_id, g1_files, g2_files, missing]
        }
        .set { ch_rmats_files_by_comparison_with_missing }

    // The samples that produced no .rmats files must be excluded from the BAM
    // lists and sample ids too — otherwise b1.txt/b1_samples.txt would declare
    // samples that have no counts in the rMATS merge.
    ch_rmats_files_by_comparison_with_missing
        .map { comparison_id, g1_files, g2_files, missing ->
            [comparison_id, g1_files, g2_files]
        }
        .set { ch_rmats_files_by_comparison }

    ch_rmats_files_by_comparison_with_missing
        .map { comparison_id, _g1, _g2, missing ->
            [comparison_id, missing]
        }
        .set { ch_missing_by_comparison }

    /*
     * Single groupTuple for BAMs and sample_ids: the two lists share the same
     * internal order by construction. Nextflow only guarantees alignment
     * *within* one groupTuple; two independent groupTuple calls (as used
     * before) may order their lists differently, silently crossing the
     * column labels of IncLevel1/IncLevel2 against the real rMATS columns.
     */
    ch_bams_grouped
        .groupTuple(by: 0)  // Group by comparison_id
        .join(ch_missing_by_comparison)
        .map { comparison_id, groups, sample_ids, bams, missing_ids ->
            def g1_bams = []
            def g2_bams = []
            def g1_ids  = []
            def g2_ids  = []

            groups.eachWithIndex { group, idx ->
                if (sample_ids[idx] in missing_ids) return
                if (group == 1) {
                    g1_bams.add(bams[idx])
                    g1_ids.add(sample_ids[idx])
                } else {
                    g2_bams.add(bams[idx])
                    g2_ids.add(sample_ids[idx])
                }
            }

            // b1.txt/b2.txt order and sample ids come from the same tuple,
            // so column index i in IncLevelN maps to gN_ids[i] by construction
            [comparison_id, g1_bams, g2_bams, g1_ids, g2_ids]
        }
        .set { ch_bams_and_ids }
    
    ch_sample_ids_by_comparison = ch_bams_and_ids.map { c, _g1, _g2, i1, i2 -> [c, i1, i2] }
    
    // Join .rmats files with BAMs and their ids for each comparison
    ch_rmats_files_by_comparison
        .join(ch_bams_and_ids)
        .map { comparison_id, g1_rmats, g2_rmats, g1_bams, g2_bams, g1_ids, g2_ids ->
            [comparison_id, g1_rmats, g2_rmats, g1_bams, g2_bams, g1_ids, g2_ids]
        }
        .set { ch_rmats_post_input }
    
    /*
     * Run RMATS POST for each comparison
     */
    RMATS_POST(
        ch_rmats_post_input.map { comparison_id, _g1_rmats, _g2_rmats, _g1_bams, _g2_bams, _g1_ids, _g2_ids -> comparison_id },
        ch_rmats_post_input.map { _comparison_id, g1_rmats, _g2_rmats, _g1_bams, _g2_bams, _g1_ids, _g2_ids -> g1_rmats },  // List<Path> — staged
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, g2_rmats, _g1_bams, _g2_bams, _g1_ids, _g2_ids -> g2_rmats },  // List<Path> — staged
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, g1_bams, _g2_bams, _g1_ids, _g2_ids -> g1_bams },  // List<Path> — symlinked
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, _g1_bams, g2_bams, _g1_ids, _g2_ids -> g2_bams },  // List<Path> — symlinked
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, _g1_bams, _g2_bams, g1_ids, _g2_ids -> g1_ids },  // order matches b1.txt
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, _g1_bams, _g2_bams, _g1_ids, g2_ids -> g2_ids },  // order matches b2.txt
        gtf
    )
    
    emit:
    results    = RMATS_POST.out.results        // [comparison_id, results_dir]
    sample_ids = ch_sample_ids_by_comparison  // [comparison_id, g1_ids, g2_ids]
    versions   = RMATS_PREP.out.versions.mix(RMATS_POST.out.versions)
}
