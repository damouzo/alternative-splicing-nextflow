process LEAFCUTTER_DS {
    tag "$comparison_id"
    label 'process_high'

    publishDir "${params.outdir}/leafcutter/${comparison_id}", mode: params.publish_dir_mode

    input:
    tuple val(comparison_id), path(counts_gz), val(sample_ids), val(conditions)
    path gtf

    output:
    tuple val(comparison_id), path("${comparison_id}"), emit: results
    path "versions.yml"                               , emit: versions

    script:
    """
    mkdir -p ${comparison_id}

    # Generate groups file for leafcutter_ds.R (sample_id TAB condition)
    python3 - <<'PYEOF'
import sys
sample_ids = "${sample_ids.join(',')}".split(',')
conditions = "${conditions.join(',')}".split(',')
with open('groups.txt', 'w') as fh:
    for sid, cond in zip(sample_ids, conditions):
        fh.write(f'{sid}\\t{cond}\\n')
PYEOF

    # Generate per-comparison exon annotation from GTF for gene labeling
    Rscript - <<'REOF'
    suppressPackageStartupMessages({
        library(dplyr)
    })
    # Parse GTF to build exon table expected by leafcutter_ds.R
    gtf_lines <- readLines("${gtf}")
    gtf_data  <- gtf_lines[!startsWith(gtf_lines, '#') & nchar(gtf_lines) > 0]
    gtf_df    <- read.delim(text = paste(gtf_data, collapse = '\\n'),
                             header = FALSE, stringsAsFactors = FALSE,
                             col.names = c('seqname','source','feature','start','end',
                                           'score','strand','frame','attribute'))
    exons <- gtf_df[gtf_df\$feature == 'exon', ]
    if (nrow(exons) > 0) {
        # Extract gene_id and transcript_id from attribute column
        exons\$gene_id <- sub('.*gene_id "([^"]+)".*', '\\\\1', exons\$attribute)
        exons\$gene_name <- ifelse(
            grepl('gene_name', exons\$attribute),
            sub('.*gene_name "([^"]+)".*', '\\\\1', exons\$attribute),
            exons\$gene_id
        )
        exon_tbl <- exons[, c('seqname','start','end','strand','gene_id','gene_name')]
        colnames(exon_tbl) <- c('Chr','Start','End','Strand','gene_id','gene_name')
        write.table(exon_tbl, 'exons.txt', sep = '\\t', row.names = FALSE, quote = FALSE)
    }
    REOF

    EXON_ARG=""
    [ -f "exons.txt" ] && EXON_ARG="--exons exons.txt"

    Rscript \$(Rscript -e "cat(system.file('scripts', 'leafcutter_ds.R', package='leafcutter'))") \\
        --num_threads ${task.cpus} \\
        --output_prefix ${comparison_id}/${comparison_id} \\
        \$EXON_ARG \\
        ${counts_gz} \\
        groups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(Rscript -e "cat(as.character(packageVersion('leafcutter')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
