#!/usr/bin/env nextflow
/*
 * ========================================================================================
 *  alternative-splicing-nextflow: Multi-tool differential alternative splicing pipeline
 * ========================================================================================
 *  Github: https://github.com/damouzo/alternative-splicing-nextflow
 * ----------------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

/*
 * ========================================================================================
 *  HELP MESSAGE
 * ========================================================================================
 */

def helpMessage() {
    log.info"""
    ============================================================
     alternative-splicing-nextflow v${workflow.manifest.version}
    ============================================================
    
    Usage:
      nextflow run main.nf --input samplesheet.csv --comparisons comparisons.csv --gtf ref.gtf [options]
    
    Required Arguments:
      --input              Path to samplesheet CSV
      --comparisons        Path to comparisons CSV (group1,group2)
      --gtf                Path to GTF annotation (same as used in nf-core/rnaseq)
      --strandedness       Library strandedness: unstranded, forward, or reverse
      --read_length        Read length after trimming (from nf-core/rnaseq)
    
    Reference Files:
      --genome_fasta       Genome FASTA (required when --use_gffread true)
      --transcript_fasta   Transcript FASTA (optional when --use_gffread true)
      --use_gffread        Extract transcript FASTA via gffread (recommended) [default: false]
      --genome_build       Genome build string for MAJIQ [default: hg38]
    
    Tool Switches:
      --run_rmats          Run rMATS-turbo [default: true]
      --run_majiq          Run MAJIQ [default: true]
      --run_isar           Run IsoformSwitchAnalyzeR [default: true]
    
    MAJIQ Options:
      --majiq_license      Path to MAJIQ license file (required if --run_majiq true)

    rMATS Options:
      --rmats_novel_ss     Detect unannotated (de novo) splice sites [default: true]
    
    ISAR Options:
      --run_isar_full_annotation  Pfam domains + IUPred3 IDRs (Tier A) [default: true]
      --pfam_hmm           Path to a hmmpressed Pfam-A.hmm (required when the above is true)
      --pfam_clans         Path to Pfam-A.clans.tsv from the same release (required with the above)
      --run_isar_iupred    Enable IUPred3 IDR prediction (needs iupred3 in the isar image) [default: false]

    Sashimi Options:
      --sashimi_gene_must  Gene list (one symbol per line) to always plot in every contrast

    Output:
      --outdir             Output directory [default: ./results]
      --publish_dir_mode   Mode for publishDir: copy, symlink, link [default: copy]
    
    Other:
      --help               Show this help message
    
    For full documentation, see: docs/usage.md
    ============================================================
    """.stripIndent()
}

def validateParams() {
    if (!params.input) {
        exit 1, "ERROR: --input is required"
    }
    if (!params.comparisons) {
        exit 1, "ERROR: --comparisons is required"
    }
    if (!params.gtf) {
        exit 1, "ERROR: --gtf is required"
    }
    if (!params.strandedness) {
        exit 1, "ERROR: --strandedness is required (unstranded, forward, or reverse)"
    }
    if (!params.read_length) {
        exit 1, "ERROR: --read_length is required"
    }

    if (!(params.strandedness in ['unstranded', 'forward', 'reverse'])) {
        exit 1, "ERROR: --strandedness must be one of: unstranded, forward, reverse"
    }

    if (params.use_gffread && !params.genome_fasta) {
        exit 1, "ERROR: --genome_fasta is required when --use_gffread is true"
    }

    if (params.run_majiq && !params.majiq_license) {
        log.warn "WARNING: --majiq_license not provided. MAJIQ analysis will fail without a valid license file."
    }

    if (params.run_isar && params.run_isar_full_annotation && !params.pfam_hmm) {
        exit 1, "ERROR: --pfam_hmm is required when --run_isar_full_annotation is true. Provide a hmmpressed Pfam-A.hmm, or run with --run_isar_full_annotation false to skip protein-domain consequences."
    }
    if (params.run_isar && params.run_isar_full_annotation && params.pfam_hmm &&
        !file("${params.pfam_hmm}.h3f").exists()) {
        exit 1, "ERROR: --pfam_hmm must point to an hmmpressed Pfam-A.hmm (missing the .h3f/.h3i/.h3m/.h3p index siblings next to it). Run: hmmpress ${params.pfam_hmm}"
    }
    if (params.run_isar && params.run_isar_full_annotation && params.pfam_hmm && !params.pfam_clans) {
        exit 1, "ERROR: --pfam_clans is required when --run_isar_full_annotation is true. Provide the Pfam-A.clans.tsv map (accession -> clan) from the same Pfam release as --pfam_hmm."
    }
    if (params.run_isar && params.run_isar_full_annotation && params.pfam_clans &&
        !file(params.pfam_clans).exists()) {
        exit 1, "ERROR: --pfam_clans does not exist: ${params.pfam_clans}"
    }
}

def logRunSummary() {
    log.info """
============================================================
 alternative-splicing-nextflow v${workflow.manifest.version}
============================================================
Input/Output:
  input             : ${params.input}
  comparisons       : ${params.comparisons}
  outdir            : ${params.outdir}

References:
  gtf               : ${params.gtf}
  genome_fasta      : ${params.genome_fasta ?: 'not provided'}
  transcript_fasta  : ${params.transcript_fasta ?: 'not provided'}
  use_gffread       : ${params.use_gffread}
  genome_build      : ${params.genome_build}

Library Properties:
  strandedness      : ${params.strandedness}
  read_length       : ${params.read_length}

Tool Selection:
  run_rmats         : ${params.run_rmats}
  run_majiq         : ${params.run_majiq}
  run_isar          : ${params.run_isar}

MAJIQ:
  majiq_license     : ${params.majiq_license ?: 'not provided'}

Deliverables (shippable layer):
  publish_deliverables            : ${params.publish_deliverables}
  publish_level                   : ${params.publish_level}
  publish_rmats_jcec              : ${params.publish_rmats_jcec}
  publish_rmats_individual_counts : ${params.publish_rmats_individual_counts}
  results_contract_version        : ${params.results_contract_version}
============================================================
"""
}

/*
 * ========================================================================================
 *  IMPORT WORKFLOWS
 * ========================================================================================
 */

include { ALTERNATIVE_SPLICING } from './workflows/alternative_splicing'

/*
 * ========================================================================================
 *  MAIN WORKFLOW
 * ========================================================================================
 */

workflow {
  if (params.help) {
    helpMessage()
    exit 0
  }

  validateParams()
  logRunSummary()

    ALTERNATIVE_SPLICING()
}
