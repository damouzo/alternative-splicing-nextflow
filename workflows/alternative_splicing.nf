/*
 * ========================================================================================
 *  ALTERNATIVE_SPLICING: Main workflow
 * ========================================================================================
 */

// Import subworkflows
include { INPUT_CHECK             } from '../subworkflows/local/input_check/main'
include { RMATS_ANALYSIS          } from '../subworkflows/local/rmats_analysis/main'
include { MAJIQ_ANALYSIS          } from '../subworkflows/local/majiq_analysis/main'
include { ISOFORMSWITCHR_ANALYSIS } from '../subworkflows/local/isoformswitchr_analysis/main'
include { SASHIMI_ANALYSIS        } from '../subworkflows/local/sashimi_analysis/main'
include { PEGASAS_ANALYSIS        } from '../subworkflows/local/pegasas_analysis/main'
include { LEAFCUTTER_ANALYSIS     } from '../subworkflows/local/leafcutter_analysis/main'

// Import modules
include { RENDER_REPORT } from '../modules/local/report/main'
include { PEGASAS_GROUPS } from '../modules/local/pegasas_groups/main'
include { EXPORT_METADATA    } from '../modules/local/export_metadata/main'
include { RMATS_MASTER       } from '../modules/local/rmats_master/main'
include { EXPORT_TOOL_MASTERS } from '../modules/local/export_tool_masters/main'
include { CROSS_TOOL_MASTER   } from '../modules/local/export_tool_masters/main'
include { SASHIMI_INDEX      } from '../modules/local/sashimi_index/main'
include { CONTRAST_MANIFEST  } from '../modules/local/contrast_manifest/main'
include { VALIDATE_RESULTS   } from '../modules/local/validate_results/main'

workflow ALTERNATIVE_SPLICING {

    // Input channels
    ch_samplesheet = channel.fromPath(params.input,       checkIfExists: true)
    ch_comparisons = channel.fromPath(params.comparisons, checkIfExists: true)
    ch_gtf         = channel.fromPath(params.gtf,         checkIfExists: true).first()
    ch_report_rmd  = channel.fromPath("${projectDir}/bin/render_report.Rmd", checkIfExists: true).first()

    // Placeholder directories for disabled tools — sentinel names checked by name in RENDER_REPORT
    def no_rmats_dir      = file("${workflow.projectDir}/assets/empty/NO_RMATS")
    def no_majiq_dir      = file("${workflow.projectDir}/assets/empty/NO_MAJIQ")
    def no_isar_dir       = file("${workflow.projectDir}/assets/empty/NO_ISAR")
    def no_sashimi_dir    = file("${workflow.projectDir}/assets/empty/NO_SASHIMI")
    def no_pegasas_dir    = file("${workflow.projectDir}/assets/empty/NO_PEGASAS")
    def no_leafcutter_dir = file("${workflow.projectDir}/assets/empty/NO_LEAFCUTTER")

    /*
     * SUBWORKFLOW: Input validation and channel creation
     */
    INPUT_CHECK(
        ch_samplesheet,
        ch_comparisons,
        ch_gtf
    )

    ch_samples_bam      = INPUT_CHECK.out.samples_bam
    ch_samples_salmon   = INPUT_CHECK.out.samples_salmon
    ch_comparisons_meta = INPUT_CHECK.out.comparisons

    // Fan-out comparison IDs for disabled-tool fallback branches
    ch_comparisons_meta
        .map { meta -> meta.id }
        .multiMap { comp_id ->
            rmats:      comp_id
            majiq:      comp_id
            isar:       comp_id
            sashimi:    comp_id
            pegasas:    comp_id
            leafcutter: comp_id
        }
        .set { ch_ids_split }

    /*
     * SUBWORKFLOW: rMATS-turbo (optional)
     */
    ch_rmats_for_report   = channel.empty()
    ch_rmats_for_sashimi  = channel.empty()
    ch_rmats_for_pegasas  = channel.empty()
    ch_sample_ids_for_report = channel.empty()

    if (params.run_rmats) {
        RMATS_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_rmats_for_report  = RMATS_ANALYSIS.out.results
        ch_rmats_for_sashimi = RMATS_ANALYSIS.out.results
        ch_rmats_for_pegasas = RMATS_ANALYSIS.out.results
        ch_sample_ids_for_report = RMATS_ANALYSIS.out.sample_ids
    } else {
        ch_ids_split.rmats
            .map { comp_id -> [comp_id, no_rmats_dir] }
            .set { ch_rmats_for_report }

        ch_ids_split.rmats
            .map { comp_id -> [comp_id, [], []] }
            .set { ch_sample_ids_for_report }
    }

    /*
     * SUBWORKFLOW: MAJIQ (optional)
     */
    ch_majiq_for_report = channel.empty()
    if (params.run_majiq) {
        MAJIQ_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_majiq_for_report = MAJIQ_ANALYSIS.out.results
    } else {
        ch_ids_split.majiq
            .map { comp_id -> [comp_id, no_majiq_dir] }
            .set { ch_majiq_for_report }
    }

    /*
     * SUBWORKFLOW: IsoformSwitchAnalyzeR (optional)
     */
    ch_isar_for_report = channel.empty()
    if (params.run_isar) {
        ISOFORMSWITCHR_ANALYSIS(
            ch_samples_salmon,
            ch_gtf
        )
        ch_isar_for_report = ISOFORMSWITCHR_ANALYSIS.out.results
    } else {
        ch_ids_split.isar
            .map { comp_id -> [comp_id, no_isar_dir] }
            .set { ch_isar_for_report }
    }

    /*
     * SUBWORKFLOW: Sashimi plots (optional; requires run_rmats = true)
     */
    ch_sashimi_for_report = channel.empty()
    if (params.run_sashimi && params.run_rmats) {
        SASHIMI_ANALYSIS(
            ch_rmats_for_sashimi,
            ch_samples_bam
        )
        ch_sashimi_for_report = SASHIMI_ANALYSIS.out.results
    } else {
        ch_ids_split.sashimi
            .map { comp_id -> [comp_id, no_sashimi_dir] }
            .set { ch_sashimi_for_report }
    }

    /*
     * SUBWORKFLOW: PEGASAS pathway-splicing correlation (optional; requires run_rmats = true)
     */
    ch_pegasas_for_report = channel.empty()
    if (params.run_pegasas && params.run_rmats && params.salmon_merged_tpm) {
        ch_salmon_tpm = channel.fromPath(params.salmon_merged_tpm, checkIfExists: true).first()

        def gmt_path = params.pathway_gmt ?: "${workflow.projectDir}/assets/hallmarks50.gmt.txt"
        ch_gmt = channel.fromPath(gmt_path, checkIfExists: true).first()

        if (params.pegasas_groups) {
            ch_pegasas_groups = channel.fromPath(params.pegasas_groups, checkIfExists: true)
        } else {
            ch_pegasas_groups = PEGASAS_GROUPS(channel.fromPath(params.input, checkIfExists: true))
        }

        PEGASAS_ANALYSIS(
            ch_rmats_for_pegasas,
            ch_salmon_tpm,
            ch_gmt,
            RMATS_ANALYSIS.out.sample_ids,
            ch_pegasas_groups
        )
        ch_pegasas_for_report = PEGASAS_ANALYSIS.out.results
    } else {
        ch_ids_split.pegasas
            .map { comp_id -> [comp_id, no_pegasas_dir] }
            .set { ch_pegasas_for_report }
    }

    /*
     * SUBWORKFLOW: LeafCutter intron excision differential splicing (optional)
     */
    ch_leafcutter_for_report = channel.empty()
    if (params.run_leafcutter) {
        LEAFCUTTER_ANALYSIS(
            ch_samples_bam,
            ch_gtf
        )
        ch_leafcutter_for_report = LEAFCUTTER_ANALYSIS.out.results
    } else {
        ch_ids_split.leafcutter
            .map { comp_id -> [comp_id, no_leafcutter_dir] }
            .set { ch_leafcutter_for_report }
    }

    /*
     * MODULE: Render per-comparison HTML report
     */
    ch_rmats_for_report
        .join(ch_majiq_for_report,      by: 0)
        .join(ch_isar_for_report,       by: 0)
        .join(ch_sashimi_for_report,    by: 0)
        .join(ch_pegasas_for_report,    by: 0)
        .join(ch_leafcutter_for_report, by: 0)
        .join(ch_sample_ids_for_report, by: 0)
        // Condition names (group1/group2 columns of comparisons.csv) feed the
        // report legends so plots read Healthy vs Patient, not Group1 vs Group2.
        .join(ch_comparisons_meta.map { meta -> [meta.id, meta.group1, meta.group2] }, by: 0)
        .map { comp_id, rdir, mdir, idir, sdir, pdir, ldir, g1_ids, g2_ids, g1_name, g2_name ->
            [comp_id,
             rdir.name, rdir,
             mdir.name, mdir,
             idir.name, idir,
             sdir.name, sdir,
             pdir.name, pdir,
             ldir.name, ldir,
             g1_ids, g2_ids,
             g1_name, g2_name]
        }
        .set { ch_report_inputs }

    RENDER_REPORT(ch_report_inputs, ch_report_rmd)
    ch_reports = RENDER_REPORT.out.html  // [comparison_id, html]

    /*
     * ========================================================================================
     *  DELIVERABLES: publication layer (masters, manifests, indexes, QA)
     * ========================================================================================
     */
    if (params.publish_deliverables) {

        def outdir_root = {
            params.outdir.startsWith('/') ? params.outdir :
                (workflow.launchDir.resolve().toString() + '/' + params.outdir)
        }()
        // Shippable layer: everything the recipient gets lives under deliverables/.
        // results/raw/ is audit-only. All internal paths are relative to deliverables/.
        def deliverables_root = "${outdir_root}/deliverables"

        // groovy.json package access is not resolved inside dataflow closures;
        // serialise/parse through these top-level references instead
        def json_out   = groovy.json.JsonOutput
        def json_parser = new groovy.json.JsonSlurper()
        def to_json   = { Object value -> json_out.toJson(value) }
        def from_json = { String text -> json_parser.parseText(text) }

        def pegasas_enabled = params.run_pegasas && params.run_rmats && params.salmon_merged_tpm
        def sashimi_enabled = params.run_sashimi && params.run_rmats

        def tools_record = [
            [name: 'rmats',      enabled: params.run_rmats,      container: 'docker.io/xinglab/rmats:v4.3.0'],
            [name: 'majiq',      enabled: params.run_majiq,      container: (params.majiq_sif ?: 'local/majiq:3.0')],
            [name: 'isar',       enabled: params.run_isar,       container: (params.isar_container ?: 'ghcr.io/damouzo/alternative-splicing-nextflow/isar:latest')],
            [name: 'leafcutter', enabled: params.run_leafcutter, container: (params.leafcutter_container ?: 'ghcr.io/damouzo/alternative-splicing-nextflow/leafcutter:latest')],
            [name: 'sashimi',    enabled: sashimi_enabled,       container: (params.sashimi_container ?: 'ghcr.io/damouzo/alternative-splicing-nextflow/sashimi:latest')],
            [name: 'pegasas',    enabled: pegasas_enabled,       container: (params.pegasas_container ?: 'ghcr.io/damouzo/alternative-splicing-nextflow/pegasas:latest')]
        ]

        // ---- Consolidated software versions: every process versions.yml ----
        ch_versions = INPUT_CHECK.out.versions
        if (params.run_rmats)      { ch_versions = ch_versions.mix(RMATS_ANALYSIS.out.versions) }
        if (params.run_majiq)      { ch_versions = ch_versions.mix(MAJIQ_ANALYSIS.out.versions) }
        if (params.run_isar)       { ch_versions = ch_versions.mix(ISOFORMSWITCHR_ANALYSIS.out.versions) }
        if (sashimi_enabled)       { ch_versions = ch_versions.mix(SASHIMI_ANALYSIS.out.versions) }
        if (pegasas_enabled)       { ch_versions = ch_versions.mix(PEGASAS_ANALYSIS.out.versions) }
        if (params.run_leafcutter) { ch_versions = ch_versions.mix(LEAFCUTTER_ANALYSIS.out.versions) }
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)

        // ---- Fase A: run-level metadata ----
        def pipeline_json = to_json([
            pipeline_name:   workflow.manifest.name,
            pipeline_version: workflow.manifest.version,
            nextflow_version: nextflow.version.toString(),
            session_id:       workflow.sessionId,
            execution_date:   (workflow.start ?: new Date()).format('yyyy-MM-dd HH:mm:ss'),
            outdir:           outdir_root,
            results_contract_version: params.results_contract_version,
            publish_level:    params.publish_level,
fdr_cutoff:             params.report_fdr_cutoff,
            dpsi_cutoff:            params.report_dpsi_cutoff,
            majiq_probability_threshold: params.majiq_probability_threshold,
            majiq_dpsi_cutoff:      params.majiq_delta_psi_threshold,
            strandedness:           params.strandedness,
            read_length:      params.read_length,
            organism:         params.organism,
            genome_build:     params.genome_build,
            isar_test_method: params.isar_test_method,
            rmats_cstat:      params.rmats_cstat,
            tools:            tools_record,
            tool_reliability: params.tool_reliability
        ])

        ch_metadata_payload = ch_comparisons_meta
            .map { meta -> meta.id }
            .collect()
            .map { ids ->
                def payload = from_json(pipeline_json)
                payload.comparisons = ids
                to_json([pipeline: payload])
            }

        // samples rows must flow in — gate the payload on the samples channel
        ch_samples_payload = INPUT_CHECK.out.samples_full
            .map { meta, bam, bai, sd ->
                [sample_id: meta.id, condition: meta.condition,
                 replicate: meta.replicate, comparison_id: meta.comparison_id,
                 group: meta.group, bam: bam.toString(),
                 bai: bai.toString(), salmon_dir: sd.toString()]
            }
            .collect()
            .map { rows ->
                to_json([samples: rows])
            }

        EXPORT_METADATA(
            ch_samples_payload.combine(ch_metadata_payload).map { samples_json, meta_json ->
                def payload = from_json(samples_json)
                payload.pipeline = from_json(meta_json).pipeline
                to_json(payload)
            },
            // Scripts passed as path inputs: their content is part of the task
            // hash, so edits invalidate the cache on -resume.
            file("${projectDir}/bin/build_results_manifest.py")
        )

        // ---- Fase C: rMATS master ----
        ch_rmats_master_input = ch_rmats_for_report
            .filter { comp_id, d -> d.name != 'NO_RMATS' }
            .map { comp_id, d -> [comp_id, d] }
        RMATS_MASTER(
            ch_rmats_master_input,
            params.report_fdr_cutoff,
            params.report_dpsi_cutoff,
            file("${projectDir}/bin/build_rmats_master.py")
        )

        // ---- Fase B: per-tool masters (MAJIQ / ISAR / LeafCutter / PEGASAS) ----
        def majiq_params_json = to_json([
            probability_threshold: params.majiq_probability_threshold,
            dpsi_cutoff:           params.majiq_delta_psi_threshold
        ])
        def isar_params_json = to_json([
            fdr_cutoff:  params.report_fdr_cutoff,
            dpsi_cutoff: params.report_dpsi_cutoff
        ])
        def leafcutter_params_json = to_json([
            fdr_cutoff: params.report_fdr_cutoff
        ])

        ch_export_inputs = ch_majiq_for_report
            .filter { comp_id, d -> d.name != 'NO_MAJIQ' }
            .map { comp_id, d -> [comp_id, 'majiq', d, majiq_params_json] }
        ch_export_inputs = ch_export_inputs.mix(
            ch_isar_for_report
                .filter { comp_id, d -> d.name != 'NO_ISAR' }
                .map { comp_id, d -> [comp_id, 'isar', d, isar_params_json] }
        )
        ch_export_inputs = ch_export_inputs.mix(
            ch_leafcutter_for_report
                .filter { comp_id, d -> d.name != 'NO_LEAFCUTTER' }
                .map { comp_id, d -> [comp_id, 'leafcutter', d, leafcutter_params_json] }
        )
        ch_export_inputs = ch_export_inputs.mix(
            ch_pegasas_for_report
                .filter { comp_id, d -> d.name != 'NO_PEGASAS' }
                .map { comp_id, d ->
                    def sig_path = "${outdir_root}/raw/${params.tool_ids.pegasas}/${comp_id}/${comp_id}_sig_pathways.tsv"
                    [comp_id, 'pegasas', d, to_json([
                        fdr_cutoff:   params.report_fdr_cutoff,
                        sig_pathways: sig_path
                    ])]
                }
        )

        EXPORT_TOOL_MASTERS(ch_export_inputs)

        // ---- Fase D: sashimi index ----
        ch_sashimi_index_input = ch_sashimi_for_report
            .filter { comp_id, d -> d.name != 'NO_SASHIMI' }
            .map { comp_id, d -> [comp_id, d] }
        SASHIMI_INDEX(
            ch_sashimi_index_input,
            file("${projectDir}/bin/build_sashimi_index.py")
        )

        // ---- Fase B/E: cross-tool gene master ----
        ch_master_files = channel.empty()
        if (params.build_cross_tool_master) {
            ch_master_files = RMATS_MASTER.out.tables.map { comp_id, m, _s, _su -> [comp_id, m] }
            ch_master_files = ch_master_files.mix(
                EXPORT_TOOL_MASTERS.out.tables.map { comp_id, _t, m, _s, _su -> [comp_id, m] }
            )
            CROSS_TOOL_MASTER(ch_master_files.groupTuple(by: 0))
        }

//        // ---- Fase A/B/D: per-comparison manifest ----
        // The manifest is assembled from the published deliverables on disk
        // (deterministic layout), so no channel aggregation is needed. Ordering
        // is enforced by joining per-comparison state channels that ALWAYS
        // emit one value per key — the same pattern the report section uses.
        ch_rmats_state      = params.run_rmats ?
            RMATS_MASTER.out.tables.map { comp_id, _m, _s, _su -> [comp_id, '1'] } :
            ch_ids_split.rmats.map { comp_id -> [comp_id, '0'] }
        ch_majiq_state      = params.run_majiq ?
            EXPORT_TOOL_MASTERS.out.tables.filter { _c, tool, _m, _s, _su -> tool == 'majiq' }.map { comp_id, _t, _m, _s, _su -> [comp_id, '1'] } :
            ch_ids_split.majiq.map { comp_id -> [comp_id, '0'] }
        ch_isar_state       = params.run_isar ?
            EXPORT_TOOL_MASTERS.out.tables.filter { _c, tool, _m, _s, _su -> tool == 'isar' }.map { comp_id, _t, _m, _s, _su -> [comp_id, '1'] } :
            ch_ids_split.isar.map { comp_id -> [comp_id, '0'] }
        ch_leafcutter_state = params.run_leafcutter ?
            EXPORT_TOOL_MASTERS.out.tables.filter { _c, tool, _m, _s, _su -> tool == 'leafcutter' }.map { comp_id, _t, _m, _s, _su -> [comp_id, '1'] } :
            ch_ids_split.leafcutter.map { comp_id -> [comp_id, '0'] }
        ch_pegasas_state    = pegasas_enabled ?
            EXPORT_TOOL_MASTERS.out.tables.filter { _c, tool, _m, _s, _su -> tool == 'pegasas' }.map { comp_id, _t, _m, _s, _su -> [comp_id, '1'] } :
            ch_ids_split.pegasas.map { comp_id -> [comp_id, '0'] }
        ch_sashimi_state    = sashimi_enabled ?
            SASHIMI_INDEX.out.index.map { comp_id, _idx -> [comp_id, '1'] } :
            ch_ids_split.sashimi.map { comp_id -> [comp_id, '0'] }
        ch_cross_state      = (params.build_cross_tool_master && params.run_rmats) ?
            CROSS_TOOL_MASTER.out.master.map { comp_id, _f, _gs -> [comp_id, '1'] } :
            ch_ids_split.rmats.map { comp_id -> [comp_id, '0'] }
        ch_report_state     = ch_reports.map { comp_id, _html -> [comp_id, '1'] }

        ch_manifest_gate = ch_rmats_for_report.map { comp_id, _d -> [comp_id] }
            .join(ch_rmats_state,      by: 0)
            .join(ch_majiq_state,      by: 0)
            .join(ch_isar_state,       by: 0)
            .join(ch_leafcutter_state, by: 0)
            .join(ch_pegasas_state,    by: 0)
            .join(ch_sashimi_state,    by: 0)
            .join(ch_cross_state,      by: 0)
            .join(ch_report_state,     by: 0)
            .join(ch_comparisons_meta.map { meta ->
                [meta.id, meta.group1, meta.group2]
            }, by: 0)
.map { comp_id, s_rmats, s_majiq, s_isar, s_lc, s_pegasas,
               s_sashimi, s_cross, s_report, g1, g2 ->
            [comp_id, pipeline_json, outdir_root, deliverables_root, g1, g2]
        }

        CONTRAST_MANIFEST(
            ch_manifest_gate,
            file("${projectDir}/bin/build_results_manifest.py")
        )

        // ---- Fase E: QA gate + run_info finalisation ----
        ch_versions = ch_versions
            .mix(EXPORT_METADATA.out.versions)
            .mix(RMATS_MASTER.out.versions)
            .mix(EXPORT_TOOL_MASTERS.out.versions)
            .mix(SASHIMI_INDEX.out.versions)
            .mix(CONTRAST_MANIFEST.out.versions)
        if (params.build_cross_tool_master && params.run_rmats) {
            ch_versions = ch_versions.mix(CROSS_TOOL_MASTER.out.versions)
        }
        ch_version_files = ch_versions.map { f -> f.toString() }.collect()

        ch_qa_trigger = CONTRAST_MANIFEST.out.manifest
            .map { comp_id, _mf -> comp_id }
            .collect()
            .map { comp_ids -> [comp_ids.size(), comp_ids.size()] }

        ch_validate_input = ch_qa_trigger
            .combine(ch_version_files.map { files -> files.join(':') })

        VALIDATE_RESULTS(
            ch_validate_input,
            outdir_root,
            deliverables_root,
            file("${projectDir}/bin/merge_versions.py"),
            file("${projectDir}/bin/validate_results_contract.py"),
            file("${projectDir}/assets/results_schema.yaml")
        )
    }
}
