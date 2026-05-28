/*
 * ========================================================================================
 *  ISOFORMSWITCHR_ANALYSIS: IsoformSwitchAnalyzeR workflow
 * ========================================================================================
 *  Workflow:
 *    0. GFFREAD (optional): Extract transcript FASTA from genome FASTA + GTF
 *    1. IMPORT: Import Salmon quantifications
 *    2. SWITCH_TEST: Statistical testing with satuRn
 *    3. EXTRACT_ORF: Extract open reading frames
 *    4. SWITCH_CONSEQUENCES: Analyze functional consequences
 */

include { GFFREAD_TRANSCRIPTOME    } from '../../../modules/local/gffread_transcriptome/main'
include { ISAR_WRITE_SAMPLESHEET   } from '../../../modules/local/isar_write_samplesheet/main'
include { ISAR_IMPORT              } from '../../../modules/local/isar_import/main'
include { ISAR_SWITCH_TEST         } from '../../../modules/local/isar_switch_test/main'
include { ISAR_EXTRACT_ORF         } from '../../../modules/local/isar_extract_orf/main'
include { ISAR_RUN_PFAM            } from '../../../modules/local/isar_run_pfam/main'
include { ISAR_RUN_IUPRED          } from '../../../modules/local/isar_run_iupred/main'
include { ISAR_SWITCH_CONSEQUENCES } from '../../../modules/local/isar_switch_consequences/main'

workflow ISOFORMSWITCHR_ANALYSIS {
    take:
    samples_salmon  // channel: [meta, salmon_dir] with meta.comparison_id
    gtf             // path: annotation.gtf
    
    main:
    
    /*
     * Optional: Run gffread to extract transcript FASTA
     */
    if (params.use_gffread) {
        ch_genome_fasta = channel.fromPath(params.genome_fasta, checkIfExists: true)
        
        GFFREAD_TRANSCRIPTOME(
            gtf,
            ch_genome_fasta
        )
        
        ch_transcript_fasta = GFFREAD_TRANSCRIPTOME.out.transcript_fasta.first()
        ch_gffread_versions = GFFREAD_TRANSCRIPTOME.out.versions
    } else {
        ch_transcript_fasta = channel.fromPath(params.transcript_fasta, checkIfExists: true).first()
        ch_gffread_versions = channel.empty()
    }
    
    /*
     * Group samples by comparison_id and build per-comparison CSV rows
     * salmon_dirs are kept as path objects so Nextflow can stage them in ISAR_IMPORT.
     */
    samples_salmon
        .map { meta, salmon_dir ->
            [meta.comparison_id, meta.id, meta.condition, meta.replicate, salmon_dir]
        }
        .groupTuple(by: 0)
        .map { comparison_id, sample_ids, conditions, replicates, salmon_dirs ->
            // Rows without salmon_dir — paths will be staged and appended in ISAR_IMPORT
            def rows = [sample_ids, conditions, replicates]
                .transpose()
                .collect { it.join(',') }
            [comparison_id, rows, salmon_dirs]
        }
        .set { ch_comparison_data }

    /*
     * Write per-comparison partial samplesheets (sample,condition,replicate) inside
     * proper Nextflow work directories. salmon_dir column added in ISAR_IMPORT.
     */
    ISAR_WRITE_SAMPLESHEET(
        ch_comparison_data.map { id, _rows, _dirs -> id },
        ch_comparison_data.map { _id, rows, _dirs -> rows }
    )

    // Join samplesheet output with salmon_dirs so ISAR_IMPORT gets both in sync
    ISAR_WRITE_SAMPLESHEET.out.samplesheet
        .join(ch_comparison_data.map { id, _rows, dirs -> [id, dirs] })
        .set { ch_isar_import_input }

    /*
     * Run ISAR IMPORT
     */
    ISAR_IMPORT(
        ch_isar_import_input.map { id, _f, _dirs -> id },
        ch_isar_import_input.map { _id, f, _dirs -> f },
        gtf,
        ch_transcript_fasta,
        ch_isar_import_input.map { _id, _f, dirs -> dirs }
    )
    
    /*
     * Run ISAR SWITCH TEST
     */
    ISAR_SWITCH_TEST(
        ISAR_IMPORT.out.rds
    )
    
    /*
     * Run ISAR EXTRACT ORF
     */
    ISAR_EXTRACT_ORF(
        ISAR_SWITCH_TEST.out.rds,
        gtf
    )

    /*
     * Optional: PFAM domain analysis and IUPred3 IDR prediction (Tier A full annotation)
     * Gated by --run_isar_full_annotation flag
     */
    def no_pfam_file   = file("${workflow.projectDir}/assets/empty/NO_PFAM")
    def no_iupred_file = file("${workflow.projectDir}/assets/empty/NO_IUPRED")

    ch_orf_rds = ISAR_EXTRACT_ORF.out.rds

    if (params.run_isar_full_annotation) {
        ISAR_RUN_IUPRED(ch_orf_rds)

        if (params.pfam_hmm) {
            ch_pfam_hmm = channel.fromPath(params.pfam_hmm, checkIfExists: true).first()
            ISAR_RUN_PFAM(ch_orf_rds, ch_pfam_hmm)

            ch_orf_rds
                .join(ISAR_RUN_PFAM.out.results)
                .join(ISAR_RUN_IUPRED.out.results)
                .set { ch_consequences_input }
        } else {
            ch_orf_rds
                .join(ISAR_RUN_IUPRED.out.results)
                .map { id, rds, iupred -> [id, rds, no_pfam_file, iupred] }
                .set { ch_consequences_input }
        }
    } else {
        ch_orf_rds
            .map { id, rds -> [id, rds, no_pfam_file, no_iupred_file] }
            .set { ch_consequences_input }
    }

    /*
     * Run ISAR SWITCH CONSEQUENCES
     */
    ISAR_SWITCH_CONSEQUENCES(ch_consequences_input)

    emit:
    results  = ISAR_SWITCH_CONSEQUENCES.out.results  // [comparison_id, results_dir]
    versions = ch_gffread_versions
        .mix(ISAR_IMPORT.out.versions)
        .mix(ISAR_SWITCH_TEST.out.versions)
        .mix(ISAR_EXTRACT_ORF.out.versions)
        .mix(ISAR_SWITCH_CONSEQUENCES.out.versions)
}
