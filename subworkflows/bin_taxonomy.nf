/*
 * Use GTDB-TK results to attach a taonomic assignment to each metagenome assembled genome (MAG), and then filter BAM files to get reads specific to each MAG.
 */

include { ASSIGN_TAXONOMY } from '../modules/local/assign_taxonomy/main'
include { FILTER_BAM } from '../modules/local/filter_bam/main'
include { MEDAKA as MEDAKA_MAG } from '../modules/local/medaka/main'                                                                                                                                                          

workflow BIN_TAXONOMY {

    take:
    bin_taxonomy
    bin_bam

    main:
    ch_versions = Channel.empty()
    ch_ecoli = Channel.empty() 
    ch_salmonella = Channel.empty()
    ch_lmonocytogenes = Channel.empty()
    ch_prokarya_reads = Channel.empty()

    if (params.ASSIGN_TAXONOMY.definitiontable) {
        ASSIGN_TAXONOMY(bin_taxonomy, [params.ASSIGN_TAXONOMY.ani_cutoff, params.ASSIGN_TAXONOMY.aln_frac], params.ASSIGN_TAXONOMY.definitiontable)

        ch_string_tt = ASSIGN_TAXONOMY.out.all_fasta.map{meta -> meta = [meta[0], meta[1]]}.filter(meta -> meta[1].toString().contains(','))
        ch_string_st = ASSIGN_TAXONOMY.out.all_fasta.map{meta -> meta = [meta[0], meta[1]]}.filter(meta -> !meta[1].toString().contains(','))
        ch_filter_tt = ASSIGN_TAXONOMY.out.all_fasta.join(ch_string_tt, by: [0]).map{ meta -> meta = [meta[0], meta[1]]}.transpose()
        ch_filter_st = ASSIGN_TAXONOMY.out.all_fasta.join(ch_string_st, by: [0]).map{ meta -> meta = [meta[0], meta[1]]}
        ch_test2 = ch_filter_tt.concat(ch_filter_st)

        ch_assigned_AF_25 = ch_test2.combine(bin_bam, by: [0]).unique().map{meta -> meta = [[id: meta[0].id, run_id: meta[0].run_id, barcode: meta[0].barcode, target: meta[0].target, single_end: true, is_ont: true, bin: (meta[1].getName().split(/\./)[1]), clean_id: (meta[1].getName().split(/\./)[2].replace("_", " ")), genus: (meta[1].getName().split(/\./)[2].split("_")[0]), amfr: (meta[1].getName().split(/\./)[3]), resfinder: (meta[1].getName().split(/\./)[4].replace("~", " ")), mlst_scheme: (meta[1].getName().split(/\./)[5]), krocus_scheme: (meta[1].getName().split(/\./)[6].replace("~", " ")), gene_DB: (meta[1].getName().split(/\./)[7])], meta[1], meta[4]]}.unique()

        FILTER_BAM(ch_assigned_AF_25)
        ch_versions = ch_versions.mix(FILTER_BAM.out.versions)

        ch_assigned_AF_26 = FILTER_BAM.out.sp_sub_fastq

        if (params.BIN_TAXONOMY.medaka_mag ) {
            MEDAKA_MAG(ch_assigned_AF_26)
            ch_prokarya_reads = MEDAKA_MAG.out.assembly_reads}
        else { ch_prokarya_reads = FILTER_BAM.out.sp_sub_fastq}

        ch_salmonella = ch_prokarya_reads.filter({meta, fasta, bam -> meta.genus.toLowerCase().matches("^salmonella")})
        ch_ecoli = ch_prokarya_reads.filter({meta, fasta, bam -> meta.clean_id.toLowerCase().matches("^escherichia coli")})
        ch_lmonocytogenes = ch_prokarya_reads.filter({meta, fasta, bam -> meta.clean_id.toLowerCase().matches("^listeria monocytogenes")})

    }

    emit:
    ch_ecoli
    ch_salmonella
    ch_lmonocytogenes
    ch_prokarya_reads 
    versions = ch_versions
}

