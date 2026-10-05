/*
 * ========================================================================================
 *  PEGASAS_ANALYSIS subworkflow
 *  Pathway-splicing correlations using PEGASAS (Python 3 port).
 *
 *  Architecture (matches the original PEGASAS design):
 *    - SHARED: pathway KS enrichment is computed ONCE for all samples, because
 *      a sample's KS score depends only on its expression and the gene set.
 *    - PER CONTRAST: PSI matrix + 2-group order, then correlation subsets the
 *      shared scores to the contrast samples and correlates them with PSI.
 *    - COLLATE: one all_pathways_scores.tsv per contrast + a significant-pathway
 *      summary, gathered into a single per-contrast pegasas directory.
 *    - CROSS CONTRAST: UpSet + heatmap of significant pathways across contrasts.
 *
 *  Takes:
 *   - ch_rmats_results: [comp_id, rmats_dir]
 *   - ch_salmon_tpm:    path (single, shared)
 *   - ch_gmt:           path (single, shared)
 *   - ch_sample_ids:    [comp_id, g1_ids, g2_ids] (rMATS BAM-list order)
 *   - ch_group_info:    path (global: sample<TAB>group)
 *
 *  Emits:
 *   - results: [comp_id, pegasas_out/] per contrast (for the report)
 * ========================================================================================
 */

include { PEGASAS_PREPARE          } from '../../../modules/local/pegasas_prepare/main'
include { PEGASAS_PREPARE_CONTRAST } from '../../../modules/local/pegasas_prepare_contrast/main'
include { PEGASAS_PATHWAY          } from '../../../modules/local/pegasas_pathway/main'
include { PEGASAS_CORRELATION      } from '../../../modules/local/pegasas_correlation/main'
include { PEGASAS_COLLATE          } from '../../../modules/local/pegasas_collate/main'
include { PEGASAS_CROSS_CONTRAST   } from '../../../modules/local/pegasas_cross_contrast/main'

workflow PEGASAS_ANALYSIS {

    take:
    ch_rmats_results  // [comp_id, rmats_dir]
    ch_salmon_tpm     // path
    ch_gmt            // path
    ch_sample_ids     // [comp_id, g1_ids, g2_ids]
    ch_group_info     // path (global)

    main:

    // === SHARED: pathway activity scores computed once for all samples ===
    PEGASAS_PREPARE(ch_salmon_tpm, ch_group_info, file("${projectDir}/bin/prepare_pegasas_gene_matrix.py"))
    PEGASAS_PATHWAY(
        PEGASAS_PREPARE.out.gene_exp,
        PEGASAS_PREPARE.out.group_info,
        ch_gmt
    )

    // === PER CONTRAST: PSI matrix + 2-group order ===
    ch_contrast = ch_rmats_results
        .map { comp_id, rmats_dir ->
            def se_file = file("${rmats_dir}/SE.MATS.JC.txt")
            def b1_file = file("${rmats_dir}/b1_samples.txt")
            def b2_file = file("${rmats_dir}/b2_samples.txt")
            // If the materialized column->sample mapping is missing (e.g. rMATS
            // results generated before this fix), fail here with a clear message
            // instead of a generic staging error in PEGASAS_PREPARE_CONTRAST.
            if (!b1_file.exists() || !b2_file.exists()) {
                error "[PEGASAS_ANALYSIS] ${comp_id}: missing ${rmats_dir}/b1_samples.txt or b2_samples.txt — rMATS POST now writes these files; rerun the comparison with the current pipeline so PEGASAS can verify the column->sample order."
            }
            [comp_id, se_file, b1_file, b2_file]
        }
        .combine(ch_sample_ids, by: 0)
        .combine(ch_group_info)
        .map { comp_id, se_file, b1_samples, b2_samples, g1_ids, g2_ids, grp ->
            [comp_id, se_file, b1_samples, b2_samples, grp, g1_ids, g2_ids]
        }

    PEGASAS_PREPARE_CONTRAST(ch_contrast, file("${projectDir}/bin/prepare_pegasas_inputs.py"))

    // Correlation: shared scores + contrast PSI + contrast group_order.
    // The single shared pathway_out/ is broadcast (cross-joined) to every contrast.
    ch_corr = PEGASAS_PREPARE_CONTRAST.out.results
        .combine(PEGASAS_PATHWAY.out.results)
        .map { comp_id, psi, group_order, contrast_samples, pathway_out ->
            [comp_id, pathway_out, psi, group_order, contrast_samples]
        }
    PEGASAS_CORRELATION(ch_corr, file("${projectDir}/bin/cor_matrix_direct_perm.R"))

    // Collate: shared scores + contrast correlation_out + contrast sample ids.
    ch_collate = PEGASAS_CORRELATION.out.results
        .combine(ch_sample_ids, by: 0)
        .combine(PEGASAS_PATHWAY.out.results)
        .map { comp_id, correlation_out, g1_ids, g2_ids, pathway_out ->
            [comp_id, correlation_out, pathway_out, g1_ids, g2_ids]
        }
    PEGASAS_COLLATE(ch_collate, file("${projectDir}/bin/pegasas_collate_scores.py"))

    // === CROSS CONTRAST: upset/heatmap of significant pathways ===
    ch_cross = PEGASAS_COLLATE.out.sig_pathways
        .map { comp_id, sig_file -> sig_file }
        .collect()
    PEGASAS_CROSS_CONTRAST(ch_cross, file("${projectDir}/bin/pegasas_cross_contrast.R"))

    emit:
    results  = PEGASAS_COLLATE.out.results
    versions = PEGASAS_PREPARE.out.versions
        .mix(PEGASAS_PATHWAY.out.versions)
        .mix(PEGASAS_PREPARE_CONTRAST.out.versions)
        .mix(PEGASAS_CORRELATION.out.versions)
        .mix(PEGASAS_COLLATE.out.versions)
        .mix(PEGASAS_CROSS_CONTRAST.out.versions)
}
