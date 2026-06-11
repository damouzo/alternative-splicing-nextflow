/*
 * ========================================================================================
 *  PEGASAS_ANALYSIS subworkflow
 *  Computes pathway-splicing correlations using PEGASAS (Python 3 port).
 *
 *  Takes:
 *   - ch_rmats_results:  [comp_id, rmats_dir]
 *   - ch_salmon_tpm:     path  (single file, shared across comparisons)
 *   - ch_gmt:            path  (GMT gene signature file)
 *
 *  Emits:
 *   - results: [comp_id, correlation_dir]
 * ========================================================================================
 */

include { PEGASAS_PREPARE     } from '../../../modules/local/pegasas_prepare/main'
include { PEGASAS_PATHWAY     } from '../../../modules/local/pegasas_pathway/main'
include { PEGASAS_CORRELATION } from '../../../modules/local/pegasas_correlation/main'

workflow PEGASAS_ANALYSIS {

    take:
    ch_rmats_results  // [comp_id, rmats_dir]
    ch_salmon_tpm     // path
    ch_gmt            // path

    main:
    ch_group_info = channel.fromPath(params.pegasas_groups, checkIfExists: true)

    // Extract SE.MATS.JC.txt from each rMATS results directory
    ch_rmats_se = ch_rmats_results
        .map { comp_id, rmats_dir ->
            def se_file = file("${rmats_dir}/SE.MATS.JC.txt")
            [comp_id, se_file]
        }

    // Join: [comp_id, salmon_tpm, se_file, group_info]
    ch_prepare_input = ch_rmats_se
        .combine(ch_salmon_tpm)
        .combine(ch_group_info)
        .map { comp_id, se_file, tpm, grp ->
            [comp_id, tpm, se_file, grp]
        }

    PEGASAS_PREPARE(ch_prepare_input)

    PEGASAS_PATHWAY(
        PEGASAS_PREPARE.out.inputs,
        ch_gmt
    )

    PEGASAS_CORRELATION(PEGASAS_PATHWAY.out.results)

    emit:
    results = PEGASAS_CORRELATION.out.results
}
