process ISAR_RUN_PFAM {
    tag "$comparison_id"
    label 'process_high'

    // Container resolved from modules.config (same isar container that has hmmer installed)

    input:
    tuple val(comparison_id), path(rds_input)
    // Pfam-A.hmm plus its hmmpressed index siblings (.h3f/.h3i/.h3m/.h3p).
    // Staged as symlinks (stageInMode in modules.config) so the ~7 GB DB is not
    // copied per task; the container bind for the source dir is added in the
    // singularity/apptainer profiles.
    path pfam_db
    // Pfam-A.clans.tsv (accession -> clan). analyzePFAM() needs a clan column.
    path pfam_clans
    // domtblout -> ISAR pfam table converter. Input so edits invalidate the cache.
    path convert_script

    output:
    tuple val(comparison_id), path("${comparison_id}_pfam.txt"), emit: results
    path "versions.yml"                                        , emit: versions

    script:
    def pfam_main = (pfam_db instanceof List || pfam_db instanceof Object[])
        ? (pfam_db.find { p -> p.name.endsWith('.hmm') && !p.name.contains('.h3') } ?: pfam_db[0])
        : pfam_db
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

    if [ -s "\$AA_FASTA" ]; then
        hmmscan \\
            --domtblout ${comparison_id}_pfam.domtblout \\
            --cpu ${task.cpus} \\
            --cut_ga \\
            ${pfam_main} \\
            \$AA_FASTA > /dev/null
    else
        : > ${comparison_id}_pfam.domtblout
    fi

    # analyzePFAM() reads the pfam_scan.pl-style table, not raw hmmscan domtblout
    python3 ${convert_script} \\
        --domtblout ${comparison_id}_pfam.domtblout \\
        --clans ${pfam_clans} \\
        --output ${comparison_id}_pfam.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h 2>&1 | head -2 | tail -1 | sed 's/# HMMER //')
    END_VERSIONS
    """
}
