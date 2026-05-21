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
include { ISAR_SWITCH_CONSEQUENCES } from '../../../modules/local/isar_switch_consequences/main'

workflow ISOFORMSWITCHR_ANALYSIS {
    take:
    samples_salmon  // channel: [meta, salmon_dir] with meta.comparison_id
    _comparisons    // channel: [comparison_meta]
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
     * Computing strings in .map{} is fine — writing files here is not.
     */
    samples_salmon
        .map { meta, salmon_dir ->
            [meta.comparison_id, meta.id, meta.condition, meta.replicate, salmon_dir.toString()]
        }
        .groupTuple(by: 0)
        .map { comparison_id, sample_ids, conditions, replicates, salmon_dirs ->
            // Build list of CSV row strings; actual file writing happens in ISAR_WRITE_SAMPLESHEET
            def rows = [sample_ids, conditions, replicates, salmon_dirs]
                .transpose()
                .collect { it.join(',') }
            [comparison_id, rows]
        }
        .set { ch_comparison_data }

    /*
     * Write per-comparison samplesheets inside proper Nextflow work directories
     */
    ISAR_WRITE_SAMPLESHEET(
        ch_comparison_data.map { id, _rows -> id },
        ch_comparison_data.map { _id, rows -> rows }
    )

    /*
     * Run ISAR IMPORT
     */
    ISAR_IMPORT(
        ISAR_WRITE_SAMPLESHEET.out.samplesheet.map { id, _f -> id },
        ISAR_WRITE_SAMPLESHEET.out.samplesheet.map { _id, f -> f },
        gtf,
        ch_transcript_fasta
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
     * Run ISAR SWITCH CONSEQUENCES
     */
    ISAR_SWITCH_CONSEQUENCES(
        ISAR_EXTRACT_ORF.out.rds
    )
    
    emit:
    results  = ISAR_SWITCH_CONSEQUENCES.out.results  // [comparison_id, results_dir]
    versions = ch_gffread_versions
        .mix(ISAR_IMPORT.out.versions)
        .mix(ISAR_SWITCH_TEST.out.versions)
        .mix(ISAR_EXTRACT_ORF.out.versions)
        .mix(ISAR_SWITCH_CONSEQUENCES.out.versions)
}
