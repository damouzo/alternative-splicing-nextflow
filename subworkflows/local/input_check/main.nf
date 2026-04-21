/*
 * ========================================================================================
 *  INPUT_CHECK: Validate samplesheet and create channels
 * ========================================================================================
 */

include { VALIDATE_INPUT } from '../../../modules/local/validate_input/main'

workflow INPUT_CHECK {
    take:
    samplesheet  // path: samplesheet.csv
    comparisons  // path: comparisons.csv
    gtf          // path: annotation.gtf
    
    main:
    
    // Validate samplesheet with critical checks
    VALIDATE_INPUT(
        samplesheet,
        gtf,
        params.read_length
    )
    
    // Parse validated samplesheet and create channels
    VALIDATE_INPUT.out.csv
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id:        row.sample,
                condition: row.condition,
                replicate: row.replicate.toInteger()
            ]
            [meta, file(row.bam), file(row.bai), file(row.salmon_dir)]
        }
        .set { ch_samples_all }
    
    // Separate into BAM and Salmon channels
    ch_samples_all
        .map { meta, bam, bai, _salmon_dir ->
            [meta, bam, bai]
        }
        .set { ch_samples_bam }
    
    ch_samples_all
        .map { meta, _bam, _bai, salmon_dir ->
            [meta, salmon_dir]
        }
        .set { ch_samples_salmon }
    
    // Parse comparisons file
    comparisons
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def comparison_meta = [
                id:     "${row.group1}_vs_${row.group2}",
                group1: row.group1,
                group2: row.group2
            ]
            comparison_meta
        }
        .set { ch_comparisons }
    
    // Add comparison metadata to sample channels
    // For each sample, determine which comparison it belongs to and which group
    ch_samples_bam
        .combine(ch_comparisons)
        .filter { meta, _bam, _bai, comp ->
            meta.condition == comp.group1 || meta.condition == comp.group2
        }
        .map { meta, bam, bai, comp ->
            def group_number = (meta.condition == comp.group1) ? 1 : 2
            def meta_updated = meta + [
                comparison_id: comp.id,
                group: group_number
            ]
            [meta_updated, bam, bai]
        }
        .set { ch_samples_bam_with_comparison }
    
    ch_samples_salmon
        .combine(ch_comparisons)
        .filter { meta, _salmon_dir, comp ->
            meta.condition == comp.group1 || meta.condition == comp.group2
        }
        .map { meta, salmon_dir, comp ->
            def group_number = (meta.condition == comp.group1) ? 1 : 2
            def meta_updated = meta + [
                comparison_id: comp.id,
                group: group_number
            ]
            [meta_updated, salmon_dir]
        }
        .set { ch_samples_salmon_with_comparison }
    
    emit:
    samples_bam    = ch_samples_bam_with_comparison  // [meta, bam, bai]
    samples_salmon = ch_samples_salmon_with_comparison // [meta, salmon_dir]
    comparisons    = ch_comparisons                    // [comparison_meta]
    versions       = VALIDATE_INPUT.out.versions
}
