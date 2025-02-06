/*
 * Perform per metagenome assembled genome (MAG) quality control and typing
 */

include { CHECKM_LINEAGEWF } from '../modules/nf-core/checkm/lineagewf/main'                                                                                                                                      
include { QUAST } from '../modules/local/quast/main'                                                                                                                                      
include { ASSEMBLY_STATS } from '../modules/local/assembly_stats/main'
include { GTDBTK_CLASSIFYWF } from '../modules/nf-core/gtdbtk/classifywf/main'
include { PLOT_BINS } from '../modules/local/plot_bins/main'
include { ASSIGN_TAXONOMY } from '../modules/local/assign_taxonomy/main'
include { FILTER_BAM } from '../modules/local/filter_bam/main'
include { MEDAKA as MEDAKA_MAG } from '../modules/local/medaka/main'   

workflow BIN_QC {

    take: 
    assembled_bins
    contig_qc_results

    main:
    ch_versions = Channel.empty()
    ch_int_3 = Channel.empty()
    ch_bin_taxonomy = Channel.empty()
    ch_binreports = Channel.empty()

    if ( params.GTDBTK_CLASSIFYWF.gtdb_db) {
        if ( params.GTDBTK_CLASSIFYWF.mash_db) {
            GTDBTK_CLASSIFYWF(assembled_bins, [params.GTDBTK_CLASSIFYWF.release, params.GTDBTK_CLASSIFYWF.gtdb_db], params.GTDBTK_CLASSIFYWF.mash_db)
            ch_versions = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions)

            ch_bin_taxonomy = assembled_bins.join(GTDBTK_CLASSIFYWF.out.summary, by: [0])}}

    if ( params.CHECKM_LINEAGEWF.db) {
        CHECKM_LINEAGEWF(assembled_bins, "fasta", params.CHECKM_LINEAGEWF.db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
    }

//    QUAST(assembled_bins)
//    ch_versions = ch_versions.mix(QUAST.out.versions)

    ASSEMBLY_STATS(assembled_bins)
    ch_versions = ch_versions.mix(ASSEMBLY_STATS.out.versions)

    if ( params.GTDBTK_CLASSIFYWF.gtdb_db) {
        if ( params.GTDBTK_CLASSIFYWF.mash_db) {
            if ( params.CHECKM_LINEAGEWF.db) {
                if ( params.TAXONOMIC_PROFILING.gtdb_metadata) {
                    ch_int_1 = CHECKM_LINEAGEWF.out.checkm_tsv.join(GTDBTK_CLASSIFYWF.out.summary, by: [0])
                    ch_int_2 = ch_int_1.join(ASSEMBLY_STATS.out.asm_stats, by: [0])
                    ch_int_3 = contig_qc_results.join(ch_int_2, by: [0])

                    PLOT_BINS(ch_int_3, params.TAXONOMIC_PROFILING.gtdb_metadata, params.PLOT_BINS.template) 
                    ch_versions = ch_versions.mix(PLOT_BINS.out.versions)

                    ch_binreports = ch_binreports.mix(PLOT_BINS.out.bin_summary)
                    ch_binreports = ch_binreports.mix(PLOT_BINS.out.contig_summary) }}}}

    emit:
    ch_int_3
    ch_bin_taxonomy
    ch_binreports = ch_binreports.groupTuple(by: [0])
    versions = ch_versions
}
