/*
 * Perform read quality control (short reads removed, quality and adapter trimming)
 */

include { NANOPLOT as NANOPLOT_PREQC } from '../modules/nf-core/nanoplot/main'
include { PORECHOP_PORECHOP } from '../modules/nf-core/porechop/porechop/main'
include { FILTLONG } from '../modules/nf-core/filtlong/main'
include { SEQTK_FQCHK as SEQTK_FQCHK_PREQC } from '../modules/local/seqtk/fqchk/main'

workflow READ_QC {

    take:
    reads    // channel: [ val(meta), path(reads) ]

    main:
    // Takes input reads, runs two read quality control tools and then trims adapters and filters out low quality reads

    ch_versions = Channel.empty()

    // Run Nanoplot of pre-quality controlled reads
    NANOPLOT_PREQC(reads)
    ch_versions = ch_versions.mix(NANOPLOT_PREQC.out.versions)

    // Get nucloetide composition of reads
    SEQTK_FQCHK_PREQC(reads, params.SEQTK_FQCHK.endseq_len)
    ch_versions = ch_versions.mix(SEQTK_FQCHK_PREQC.out.versions)

    // Trims adapters
    PORECHOP_PORECHOP(reads)
    ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions)

    // Filters reads based on quality
    FILTLONG(PORECHOP_PORECHOP.out.reads)
    ch_versions = ch_versions.mix(FILTLONG.out.versions)

    emit:
    // Emit QC pass reads
    qc_pass_reads = FILTLONG.out.reads

    // Emit merged Seqth and Nanoplot results
    preqc_results  = SEQTK_FQCHK_PREQC.out.pbq_se.join(NANOPLOT_PREQC.out.qc_input, by: [0])

    versions = ch_versions
}
