/*
 * When given basecalled output directories, find reads, and get QC stats to run through pycoQC. 
 */

include { PYCOQC } from '../modules/nf-core/pycoqc/main'
include { FETCH_FASTQ } from '../modules/local/fetch_fastq/main'
include { FETCH_SEQSUM } from '../modules/local/fetch_seqsum/main'

workflow RUN_QC {

    take:
    scheme    // channel: [ val(meta) ]

    main:
    // If Guppy outdir(s) is provided as input, find FASTQ files and fetch Sequencing summary file

    ch_versions = Channel.empty()

    scheme.map{meta -> meta.run_id}.unique().set{ schemex }

    // Gather Pass FASTQ files
    FETCH_FASTQ(scheme)
    ch_versions = ch_versions.mix(FETCH_FASTQ.out.versions)

    // Get Sequencing summary file
    FETCH_SEQSUM(schemex)
    ch_versions = ch_versions.mix(FETCH_SEQSUM.out.versions)

    // Run PycoQC on Sequencing summary file
    PYCOQC(FETCH_SEQSUM.out.sequencing_summary)
    ch_versions = ch_versions.mix(PYCOQC.out.versions)

    emit:
    // Emit combined FASTQ file
    aggregated_reads = FETCH_FASTQ.out.reads

    versions = ch_versions

}
