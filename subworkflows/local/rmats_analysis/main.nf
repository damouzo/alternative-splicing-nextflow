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
    _comparisons // channel: [comparison_meta]
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
            [meta.comparison_id, meta.group, bam]
        }
        .set { ch_bams_grouped }
    
    /*
     * Collect all files per comparison
     * Group by comparison_id, then separate group1 and group2
     */
    ch_rmats_prep_grouped
        .groupTuple(by: 0)  // Group by comparison_id
        .map { comparison_id, groups, _metas, rmats_files_list ->
            // Separate group 1 and group 2
            def g1_files = []
            def g2_files = []
            
            groups.eachWithIndex { group, idx ->
                // rmats_files_list[idx] may be a single Path or a List<Path>
                // depending on how many files the glob matched — normalise to list
                def files = rmats_files_list[idx]
                def fileList = files instanceof List ? files : [files]
                if (group == 1) {
                    g1_files.addAll(fileList)
                } else {
                    g2_files.addAll(fileList)
                }
            }
            
            // Keep as Path objects — Nextflow will stage them properly in RMATS_POST work dir
            [comparison_id, g1_files, g2_files]
        }
        .set { ch_rmats_files_by_comparison }
    
    // Similarly for BAMs
    ch_bams_grouped
        .groupTuple(by: 0)  // Group by comparison_id
        .map { comparison_id, groups, bams ->
            def g1_bams = []
            def g2_bams = []
            
            groups.eachWithIndex { group, idx ->
                if (group == 1) {
                    g1_bams.add(bams[idx])
                } else {
                    g2_bams.add(bams[idx])
                }
            }
            
            // Keep as Path objects — staged as symlinks, no large file copies
            [comparison_id, g1_bams, g2_bams]
        }
        .set { ch_bams_by_comparison }
    
    // Join .rmats files and BAMs for each comparison
    ch_rmats_files_by_comparison
        .join(ch_bams_by_comparison)
        .map { comparison_id, g1_rmats, g2_rmats, g1_bams, g2_bams ->
            [comparison_id, g1_rmats, g2_rmats, g1_bams, g2_bams]
        }
        .set { ch_rmats_post_input }
    
    /*
     * Run RMATS POST for each comparison
     */
    RMATS_POST(
        ch_rmats_post_input.map { comparison_id, _g1_rmats, _g2_rmats, _g1_bams, _g2_bams -> comparison_id },
        ch_rmats_post_input.map { _comparison_id, g1_rmats, _g2_rmats, _g1_bams, _g2_bams -> g1_rmats },  // List<Path> — staged
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, g2_rmats, _g1_bams, _g2_bams -> g2_rmats },  // List<Path> — staged
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, g1_bams, _g2_bams -> g1_bams },  // List<Path> — symlinked
        ch_rmats_post_input.map { _comparison_id, _g1_rmats, _g2_rmats, _g1_bams, g2_bams -> g2_bams },  // List<Path> — symlinked
        gtf
    )
    
    emit:
    results  = RMATS_POST.out.results  // [comparison_id, results_dir]
    versions = RMATS_PREP.out.versions.mix(RMATS_POST.out.versions)
}
