/*
 * Identify prokaryotic and archael contigs and assign them to individual genomes (bins)
 * Rename binned and unbinned contigs and output metagenome assembled genomes (MAGs) and complete assemblies (MAGs and unbinned contigs)
 */

include { TIARA_TIARA } from '../modules/local/tiara/tiara/main'
include { FILTER_CONTIGS } from '../modules/local/filter_contigs/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_1 } from '../modules/nf-core/samtools/view/main'
include { SEMIBIN_SINGLEEASYBIN } from '../modules/local/semibin/singleeasybin/main'
include { METADECODER_COVERAGE } from '../modules/local/metadecoder/coverage/main'
include { METADECODER_SEED } from '../modules/local/metadecoder/seed/main'
include { METADECODER_CLUSTER } from '../modules/local/metadecoder/cluster/main'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_METADECODER } from '../modules/nf-core/dastool/fastatocontig2bin/main'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_SEMIBIN2 } from '../modules/nf-core/dastool/fastatocontig2bin/main'
include { DASTOOL_DASTOOL } from '../modules/nf-core/dastool/dastool/main'
include { FETCH_UNBINNED } from '../modules/local/fetch_unbinned/main'
include { RENAME_CONTIGS } from '../modules/local/rename_contigs/main'

workflow BIN_ASSIGNMENT {

    take:
    assembly_output

    main:
    ch_versions = Channel.empty()

    assembly_output.map { meta -> meta = [meta[0], meta[2]]}.set{ ch_fasta }
    assembly_output.map { meta -> meta = [meta[0], meta[2], meta[3], meta[4]]}.set{ ch_fasta_bam_index }
    assembly_output.map { meta -> meta = [meta[0], meta[1]]}.set{ ch_reads }

    TIARA_TIARA(ch_fasta)
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions)

    ch_tiara_fasta_bam_index = TIARA_TIARA.out.prokaryote_tiara.join(ch_fasta_bam_index, by: [0])

    FILTER_CONTIGS(ch_tiara_fasta_bam_index)
    FILTER_CONTIGS.out.filt_ref_bam.map { meta -> meta = [meta[0], meta[2], meta[3]]}.set{ ch_filter_bam_index }
    FILTER_CONTIGS.out.filt_ref_bam.map { meta -> meta = [meta[0], meta[2]]}.set{ ch_filter_bam }
    FILTER_CONTIGS.out.filt_ref_bam.map { meta -> meta = [meta[0], meta[1], meta[2]]}.set{ ch_filter_fasta_bam }
    FILTER_CONTIGS.out.filt_ref_bam.map { meta -> meta = [meta[0], meta[1]]}.set{ ch_filter_fasta }

    SEMIBIN_SINGLEEASYBIN(ch_filter_fasta_bam)
    ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)

    SAMTOOLS_VIEW_1(ch_filter_bam_index, [[],[]], [] )

    METADECODER_COVERAGE(SAMTOOLS_VIEW_1.out.sam)
    ch_versions = ch_versions.mix(METADECODER_COVERAGE.out.versions)

    METADECODER_SEED(ch_filter_fasta)
    ch_versions = ch_versions.mix(METADECODER_SEED.out.versions)

    ch_metadecoder_1 = ch_filter_fasta.join(METADECODER_COVERAGE.out.cov, by: [0])
    ch_metadecoder_2 = ch_metadecoder_1.join(METADECODER_SEED.out.seed, by: [0])

    METADECODER_CLUSTER(ch_metadecoder_2)
    ch_versions = ch_versions.mix(METADECODER_CLUSTER.out.versions)

    DASTOOL_FASTATOCONTIG2BIN_METADECODER(METADECODER_CLUSTER.out.binned_fastas, "fasta")
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_METADECODER.out.versions)

    DASTOOL_FASTATOCONTIG2BIN_SEMIBIN2(SEMIBIN_SINGLEEASYBIN.out.binned_fastas, "fa")
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_SEMIBIN2.out.versions)

    ch_dastool_input_0 = DASTOOL_FASTATOCONTIG2BIN_SEMIBIN2.out.fastatocontig2bin.concat( DASTOOL_FASTATOCONTIG2BIN_METADECODER.out.fastatocontig2bin ).groupTuple(by: [0])
    ch_dastool_input_1 = ch_fasta.join(ch_dastool_input_0, by: [0])

    DASTOOL_DASTOOL(ch_dastool_input_1, [], [])
    ch_versions = ch_versions.mix(DASTOOL_DASTOOL.out.versions)

    ch_fasta_bins = ch_fasta.join(DASTOOL_DASTOOL.out.bins, by: [0])

    FETCH_UNBINNED(ch_fasta_bins)
//    ch_versions = ch_versions.mix(GET_UNBINNED.out.versions)

    ch_split_assembly = DASTOOL_DASTOOL.out.bins.join( FETCH_UNBINNED.out.unbinned, by: [0] )

    RENAME_CONTIGS(ch_split_assembly)
    ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)

//    if (params.run_complete) {
//       assembled_bins = RENAME_CONTIGS.out.assembled_bins.concat( RENAME_CONTIGS.out.complete_assembly ).groupTuple(by: [0]) 
//    }
//    else {
    assembled_bins = RENAME_CONTIGS.out.assembled_bins

    emit:
    assembled_bins
    complete_assembly = RENAME_CONTIGS.out.complete_assembly
    reads_complete_assembly = ch_reads.join(RENAME_CONTIGS.out.complete_assembly, by: [0])
    bin_bam = ch_filter_bam
    versions = ch_versions
}

