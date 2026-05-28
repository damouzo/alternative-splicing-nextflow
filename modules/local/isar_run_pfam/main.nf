process ISAR_RUN_PFAM {
    tag "$comparison_id"
    label 'process_high'

    // Container resolved from modules.config (same isar container that has hmmer installed)

    input:
    tuple val(comparison_id), path(rds_input)
    path pfam_hmm

    output:
    tuple val(comparison_id), path("${comparison_id}_pfam.domtblout"), emit: results
    path "versions.yml"                                               , emit: versions

    script:
    """
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    # Export AA FASTA from ISAR RDS
    Rscript - <<'REOF'
    suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
    sal <- readRDS("${rds_input}")
    if (exists("exportSequences")) {
        exportSequences(sal, pathToOutput = ".",
                        writeTranscriptSequences = FALSE,
                        writeORFSequences        = FALSE,
                        writePeptideSequences    = TRUE)
    } else {
        extractSequence(sal, writeToFile = TRUE, pathToOutput = ".",
                        extractNTseq = FALSE, extractAAseq = TRUE,
                        outputPrefix = "isoformSwitchAnalyzeR",
                        onlySwitchingGenes = FALSE, quiet = TRUE)
    }
    REOF

    AA_FASTA=isoformSwitchAnalyzeR_AA.fasta

    if [ ! -s "\$AA_FASTA" ]; then
        # No sequences available — write empty domtblout
        echo "# hmmscan :: no sequences available" > ${comparison_id}_pfam.domtblout
    else
        hmmscan \\
            --domtblout ${comparison_id}_pfam.domtblout \\
            --cpu ${task.cpus} \\
            -E 1e-5 \\
            --domE 1e-5 \\
            ${pfam_hmm} \\
            \$AA_FASTA > /dev/null
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan --version 2>&1 | head -1 | sed 's/# HMMER //')
    END_VERSIONS
    """
}
