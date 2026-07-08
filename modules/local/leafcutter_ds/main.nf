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
    // Adapt LeafCutter thresholds to available replicates per group.
    def group_sizes = conditions.countBy { it }.values() as List
    def min_group_size = group_sizes ? group_sizes.min() as int : 1
    def min_samples_per_intron = Math.max(1, Math.min(5, min_group_size))
    def min_samples_per_group  = Math.max(1, Math.min(3, min_group_size))
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

    # Build exon table expected by leafcutter_ds.R without loading full GTF in memory.
    awk 'BEGIN{FS="\t"; OFS="\t"; print "chr","start","end","strand","gene_name"}
         !/^#/ && \$3=="exon" {
             attr=\$9; gene_name="";
             if (match(attr, /gene_name "[^"]+"/)) {
                 gene_name=substr(attr, RSTART+11, RLENGTH-12)
             } else if (match(attr, /gene_id "[^"]+"/)) {
                 gene_name=substr(attr, RSTART+9, RLENGTH-10)
             }
             if (gene_name != "") print \$1, \$4, \$5, \$7, gene_name
         }' ${gtf} > exons.txt

    EXON_ARG=""
    [ -s exons.txt ] && [ "\$(wc -l < exons.txt)" -gt 1 ] && EXON_ARG="--exon_file exons.txt"

    # leafcutter_ds.R lives in the repo's scripts/ dir, not installed in the R package
    Rscript /opt/leafcutter-src/scripts/leafcutter_ds.R \
        --num_threads ${task.cpus} \
        --output_prefix ${comparison_id}/${comparison_id} \
        --min_samples_per_intron ${min_samples_per_intron} \
        --min_samples_per_group ${min_samples_per_group} \
        \$EXON_ARG \
        ${counts_gz} \
        groups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(Rscript -e "cat(as.character(packageVersion('leafcutter')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
