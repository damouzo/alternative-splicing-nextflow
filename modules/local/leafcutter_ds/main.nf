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
    // Serialize as base64-encoded JSON so sample ids/conditions containing commas,
    // quotes or shell metacharacters can't break the embedded Python parsing.
    def sample_ids_b64 = groovy.json.JsonOutput.toJson(sample_ids).bytes.encodeBase64().toString()
    def conditions_b64 = groovy.json.JsonOutput.toJson(conditions).bytes.encodeBase64().toString()
    """
    mkdir -p ${comparison_id}

    SAMPLE_IDS_JSON="\$(printf '%s' '${sample_ids_b64}' | base64 -d)"
    CONDITIONS_JSON="\$(printf '%s' '${conditions_b64}' | base64 -d)"
    export SAMPLE_IDS_JSON CONDITIONS_JSON

    # Generate groups file for leafcutter_ds.R (sample_id TAB condition).
    # The order is taken from the real column header of the counts file (which
    # LEAFCUTTER_CLUSTER fixed), not from the channel list order — mapping each
    # column name to its condition via the id->condition dictionary.
    python3 - <<'PYEOF'
import gzip
import json
import os
import sys
sample_ids = json.loads(os.environ["SAMPLE_IDS_JSON"])
conditions = json.loads(os.environ["CONDITIONS_JSON"])
cond = dict(zip(sample_ids, conditions))
with gzip.open("${counts_gz}", "rt") as fh:
    header = fh.readline().split()
if not header:
    sys.stderr.write("[LEAFCUTTER_DS] counts file has no column header — aborting\\n")
    sys.exit(1)
missing = [s for s in header if s not in cond]
if missing:
    sys.stderr.write(
        "[LEAFCUTTER_DS] count columns without condition mapping: %s\\n"
        % ", ".join(missing)
    )
    sys.exit(1)
with open('groups.txt', 'w') as fh:
    for sid in header:
        fh.write(f'{sid}\\t{cond[sid]}\\n')
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
