/*
 * QC and type whole-metagenome assemblies
 */


include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_1 } from '../modules/local/minimap2/align/main'
include { GENOMAD_ENDTOEND } from '../modules/nf-core/genomad/endtoend/main'
include { SKANI_SEARCH } from '../modules/local/skani/search/main.nf'
include { SEQKIT_FX2TAB } from '../modules/nf-core/seqkit/fx2tab/main'
include { SAMTOOLS_COVERAGE } from '../modules/local/samtools/coverage/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_1 } from '../modules/nf-core/samtools/index/main'

workflow CONTIG_QC {

    take:
    reads_complete_assembly

    main:
    ch_versions = Channel.empty()

    reads_complete_assembly.map{meta -> meta = [meta[0], meta[2]]}.set { ch_complete_assembly }

    MINIMAP2_ALIGN_1(reads_complete_assembly,true,false,false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_1.out.versions)

    ch_ca_bam = MINIMAP2_ALIGN_1.out.bam.map {meta -> meta = [meta[0], meta[3]]}

    if (params.GENOMAD_ENDTOEND.db) {
        GENOMAD_ENDTOEND(ch_complete_assembly, params.GENOMAD_ENDTOEND.db)
        ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions)
    }

    if (params.SKANI_SEARCH.db) {    
        SKANI_SEARCH(ch_complete_assembly, params.SKANI_SEARCH.db)
        ch_versions = ch_versions.mix(SKANI_SEARCH.out.versions)
    }

    SEQKIT_FX2TAB(ch_complete_assembly)
    ch_versions = ch_versions.mix(SEQKIT_FX2TAB.out.versions)

    SAMTOOLS_INDEX_1(ch_ca_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_1.out.versions)

    ch_merged_index_1 = ch_complete_assembly.join(ch_ca_bam, by: [0])
    ch_merged_index_2 = ch_merged_index_1.join(SAMTOOLS_INDEX_1.out.bai, by: [0])

    SAMTOOLS_COVERAGE(ch_merged_index_2)
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    if (params.GENOMAD_ENDTOEND.db) {
        if (params.SKANI_SEARCH.db) {
            ch_p1_1 = SAMTOOLS_COVERAGE.out.coverage.join(SEQKIT_FX2TAB.out.text, by:[0])
            ch_p1_2 = ch_p1_1.join(SKANI_SEARCH.out.summary, by:[0])
            ch_plt_p1 = ch_p1_2.join(GENOMAD_ENDTOEND.out.plasmid_summary, by:[0])
        }
    }
    else (ch_plt_p1 = Channel.empty())

    emit:
    complete_assembly_bam = MINIMAP2_ALIGN_1.out.bam
    ch_plt_p1
    versions = ch_versions
}

