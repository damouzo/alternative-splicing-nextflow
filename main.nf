#!/usr/bin/env nextflow
/*
 * ========================================================================================
 *  alternative-splicing-nextflow: Multi-tool differential alternative splicing pipeline
 * ========================================================================================
 *  Github: https://github.com/BCI-KRP/alternative-splicing-nextflow
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

Resources:
  max_cpus          : ${params.max_cpus}
  max_memory        : ${params.max_memory}
  max_time          : ${params.max_time}
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
