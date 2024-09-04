/*
 * Assemble and polish metagenomes
 */

include { FLYE } from '../modules/nf-core/flye/main'                                                                                                                                                              
include { GUNZIP_FASTA } from '../modules/local/gunzip_fasta/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BIN_1 } from '../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BIN_2 } from '../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BIN_3 } from '../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BIN_4 } from '../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_BIN_5 } from '../modules/local/minimap2/align/main'
include { RACON as RACON_1 } from '../modules/local/racon/main'
include { RACON as RACON_2} from '../modules/local/racon/main'                                                                                                                                                   
include { RACON as RACON_3} from '../modules/local/racon/main'                                                                                                                                                   
include { RACON as RACON_4} from '../modules/local/racon/main'                                                                                                                                                   
include { MEDAKA } from '../modules/local/medaka/main'                                                                                                                                                          
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

workflow ASSEMBLY {

    take:
    qc_pass_reads

    main:
    ch_versions = Channel.empty()

    FLYE(qc_pass_reads, params.FLYE.read_type)
    ch_versions = ch_versions.mix(FLYE.out.versions)

    ch_assembly_polish_0 = qc_pass_reads.join(FLYE.out.fasta, by: [0])

    if (params.ASSEMBLY.racon_rounds == 0) {
        GUNZIP_FASTA(ch_assembly_polish_0)
        ch_polished_asm = GUNZIP_FASTA.out.gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }

    if (params.ASSEMBLY.racon_rounds > 0) {
        MINIMAP2_ALIGN_BIN_1(ch_assembly_polish_0,false,false,false)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_BIN_1.out.versions)

        if (params.ASSEMBLY.racon_rounds < 2) {
            RACON_1(MINIMAP2_ALIGN_BIN_1.out.paf, true)
            ch_versions = ch_versions.mix(RACON_1.out.versions)

            ch_polished_asm = RACON_1.out.improved_assembly
        }
        else {
            RACON_1(MINIMAP2_ALIGN_BIN_1.out.paf, false)
            ch_versions = ch_versions.mix(RACON_1.out.versions)
        }
    }

    if (params.ASSEMBLY.racon_rounds > 1) {
        MINIMAP2_ALIGN_BIN_2(RACON_1.out.improved_assembly,false,false,false)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_BIN_2.out.versions)

        if (params.ASSEMBLY.racon_rounds < 3) {
            RACON_2(MINIMAP2_ALIGN_BIN_2.out.paf, true)
            ch_versions = ch_versions.mix(RACON_2.out.versions)

            ch_polished_asm = RACON_2.out.improved_assembly
        }

        else {
            RACON_2(MINIMAP2_ALIGN_BIN_2.out.paf, false)
            ch_versions = ch_versions.mix(RACON_2.out.versions)
        }
    }

    if (params.ASSEMBLY.racon_rounds > 2) {
        MINIMAP2_ALIGN_BIN_3(RACON_2.out.improved_assembly,false,false,false)

        if (params.ASSEMBLY.racon_rounds < 4) {
            RACON_3(MINIMAP2_ALIGN_BIN_3.out.paf, true)
            ch_versions = ch_versions.mix(RACON_3.out.versions)

            ch_polished_asm = RACON_3.out.improved_assembly
        }

        else {
            RACON_3(MINIMAP2_ALIGN_BIN_3.out.paf, false)
            ch_versions = ch_versions.mix(RACON_3.out.versions)
        }
    }

    if (params.ASSEMBLY.racon_rounds == 4) {
        MINIMAP2_ALIGN_BIN_4(RACON_3.out.improved_assembly,false,false,false)

        RACON_4(MINIMAP2_ALIGN_BIN_4.out.paf, true)
        ch_versions = ch_versions.mix(RACON_4.out.versions)

        ch_polished_asm = RACON_4.out.improved_assembly
    }

    if (params.ASSEMBLY.medaka == true ) {
        MEDAKA(RACON_4.out.improved_assembly)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
        
        ch_polished_asm = qc_pass_reads.join(MEDAKA.out.assembly, by:[0])
    }

    MINIMAP2_ALIGN_BIN_5(ch_polished_asm, true, false, false)

    ch_polished_bam = MINIMAP2_ALIGN_BIN_5.out.bam.map{ meta -> meta = [meta[0], meta[3]]}

    SAMTOOLS_INDEX(ch_polished_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_polished_all = MINIMAP2_ALIGN_BIN_5.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])

    emit:
    final_assembly = ch_polished_all.map{meta -> meta = [meta[0], meta[2]]}
    final_assembly_mapped = ch_polished_all
    versions = ch_versions
}

