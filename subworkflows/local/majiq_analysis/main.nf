/*
 * ========================================================================================
 *  MAJIQ_ANALYSIS: MAJIQ V3 Local Splicing Variation analysis
 * ========================================================================================
 *  Workflow:
 *    1. BUILD: Build splice graph from all BAMs
 *    2. DELTAPSI: Compute posterior deltaPSI distributions per comparison
 *    3. VOILA_TSV: Export results to TSV format
 */

include { MAJIQ_PREPARE_ANNOTATION } from '../../../modules/local/majiq_prepare_annotation/main'
include { MAJIQ_BUILD              } from '../../../modules/local/majiq_build/main'
include { MAJIQ_DELTAPSI           } from '../../../modules/local/majiq_deltapsi/main'
include { MAJIQ_ORGANIZE_RESULTS   } from '../../../modules/local/majiq_organize_results/main'

workflow MAJIQ_ANALYSIS {
    take:
    samples_bam  // channel: [meta, bam, bai] with meta.comparison_id and meta.group
    gtf          // path: annotation.gtf — will be converted to GFF3 internally
    
    main:

    // Convert GTF to GFF3 — MAJIQ v3 requires GFF3 with Parent= hierarchy
    MAJIQ_PREPARE_ANNOTATION(gtf)
    // .first() converts the queue channel to a value channel so all MAJIQ_BUILD
    // invocations (one per comparison) can consume the same GFF3 without blocking.
    ch_gff3 = MAJIQ_PREPARE_ANNOTATION.out.gff3.first()

    /*
     * Group samples by comparison_id
     */
    samples_bam
        .map { meta, bam, bai ->
            [meta.comparison_id, meta.group, meta.id, bam, bai]
        }
        .groupTuple(by: 0)  // Group by comparison_id
        .set { ch_samples_by_comparison }

    /*
     * Prepare inputs for MAJIQ BUILD
     * Need: comparison_id, all_bams, all_bais, sample_info
     */
    ch_samples_by_comparison
        .map { comparison_id, _groups, sample_ids, bams, bais ->
            // sample_info carries only sample IDs — their order matches the staged bams list.
            // bam_path is no longer passed to avoid using pre-staging absolute paths in the script.
            def sample_info = sample_ids.collect { sid -> [sid] }
            [comparison_id, bams, bais, sample_info]
        }
        .set { ch_majiq_build_input }

    /*
     * Run MAJIQ BUILD
     */
    MAJIQ_BUILD(
        ch_majiq_build_input.map { comparison_id, _bams, _bais, _sample_info -> comparison_id },
        ch_majiq_build_input.map { _comparison_id, bams, _bais, _sample_info -> bams },
        ch_majiq_build_input.map { _comparison_id, _bams, bais, _sample_info -> bais },
        ch_gff3,
        ch_majiq_build_input.map { _comparison_id, _bams, _bais, sample_info -> sample_info }
    )
    
    /*
     * Separate sj files by group for DELTAPSI
     */
    MAJIQ_BUILD.out.majiq_build
        .map { comparison_id, sj_files, splicegraph ->
            [comparison_id, sj_files, splicegraph]
        }
        .set { ch_majiq_built }
    
    // Get group assignments for each sample
    samples_bam
        .map { meta, _bam, _bai ->
            [meta.comparison_id, meta.id, meta.group, meta.condition]
        }
        .groupTuple(by: 0)
        .set { ch_sample_groups }
    
    // Join build output (sj files + zarr splicegraph) with sample group metadata
    ch_majiq_built
        .join(ch_sample_groups)
        .map { comparison_id, sj_files, splicegraph, sample_ids, groups, conditions ->
            // Separate sj files by group using exact filename match — substring
            // matching is unsafe (e.g. sample1 vs sample10 would collide).
            def g1_files = []
            def g2_files = []
            def g1_name  = null
            def g2_name  = null

            sample_ids.eachWithIndex { sid, idx ->
                def sj_file = sj_files.find { it.name == "${sid}.sj" }
                if (!sj_file) {
                    error "[MAJIQ_ANALYSIS] No .sj file found for sample '${sid}' (expected: sj/${sid}.sj) in ${comparison_id}. MAJIQ_BUILD did not produce a junction file for this sample — check upstream BAM/annotation inputs."
                }
                if (groups[idx] == 1) {
                    g1_files.add(sj_file)
                    if (!g1_name) g1_name = conditions[idx]
                } else {
                    g2_files.add(sj_file)
                    if (!g2_name) g2_name = conditions[idx]
                }
            }

            [comparison_id, g1_files, g2_files, splicegraph, g1_name, g2_name]
        }
        .set { ch_majiq_deltapsi_input }
    
    /*
     * Run MAJIQ DELTAPSI
     */
    MAJIQ_DELTAPSI(
        ch_majiq_deltapsi_input.map { comparison_id, _g1_files, _g2_files, _splicegraph, _g1_name, _g2_name -> comparison_id },
        ch_majiq_deltapsi_input.map { _comparison_id, g1_files, _g2_files, _splicegraph, _g1_name, _g2_name -> g1_files },
        ch_majiq_deltapsi_input.map { _comparison_id, _g1_files, g2_files, _splicegraph, _g1_name, _g2_name -> g2_files },
        ch_majiq_deltapsi_input.map { _comparison_id, _g1_files, _g2_files, splicegraph, _g1_name, _g2_name -> splicegraph },
        ch_majiq_deltapsi_input.map { _comparison_id, _g1_files, _g2_files, _splicegraph, g1_name, _g2_name -> g1_name },
        ch_majiq_deltapsi_input.map { _comparison_id, _g1_files, _g2_files, _splicegraph, _g1_name, g2_name -> g2_name }
    )
    
    /*
     * Organise deltapsi outputs into a results directory for the report
     */
    MAJIQ_DELTAPSI.out.deltapsi
        .map { comparison_id, dpsicov, tsv ->
            [comparison_id, dpsicov, tsv]
        }
        .join(
            ch_majiq_deltapsi_input.map { comparison_id, _g1, _g2, splicegraph, _g1n, _g2n ->
                [comparison_id, splicegraph]
            }
        )
        .set { ch_organize_input }

    MAJIQ_ORGANIZE_RESULTS(
        ch_organize_input.map { comparison_id, _dpsicov, _tsv, _splicegraph -> comparison_id },
        ch_organize_input.map { _comparison_id, _dpsicov, _tsv, splicegraph -> splicegraph },
        ch_organize_input.map { _comparison_id, dpsicov, _tsv, _splicegraph -> dpsicov },
        ch_organize_input.map { _comparison_id, _dpsicov, tsv, _splicegraph -> tsv }
    )
    
    emit:
    results  = MAJIQ_ORGANIZE_RESULTS.out.results  // [comparison_id, results_dir]
    versions = MAJIQ_PREPARE_ANNOTATION.out.versions
        .mix(MAJIQ_BUILD.out.versions)
        .mix(MAJIQ_DELTAPSI.out.versions)
        .mix(MAJIQ_ORGANIZE_RESULTS.out.versions)
}
