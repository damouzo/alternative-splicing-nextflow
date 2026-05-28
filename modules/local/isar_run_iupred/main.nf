process ISAR_RUN_IUPRED {
    tag "$comparison_id"
    label 'process_medium'

    // Container resolved from modules.config (same isar container that has iupred3 installed)

    input:
    tuple val(comparison_id), path(rds_input)

    output:
    tuple val(comparison_id), path("${comparison_id}_iupred.txt"), emit: results
    path "versions.yml"                                           , emit: versions

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
        # No sequences — write empty output so analyzeIUPred2A can be skipped downstream
        touch ${comparison_id}_iupred.txt
    else
        run_iupred3.py \$AA_FASTA --output ${comparison_id}_iupred.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iupred3: \$(python3 -c "import iupred3; print(getattr(iupred3, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
