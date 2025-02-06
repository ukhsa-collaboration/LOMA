/*
 * Remove host reads from assembly using Kraken2 and/or Minimap2 against supplied database(s) of host(s).
 */

include { NANOPLOT as NANOPLOT_POSTQC} from '../modules/nf-core/nanoplot/main'
include { MINIMAP2_ALIGN } from '../modules/nf-core/minimap2/align/main'
include { FILTER_MAPPED } from '../modules/local/filter_mapped/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HOST } from '../modules/nf-core/kraken2/kraken2/main'
include { SEQKIT_GREP } from '../modules/nf-core/seqkit/grep/main'
include { SEQTK_FQCHK as SEQTK_FQCHK_POSTQC } from '../modules/local/seqtk/fqchk/main'
include { FILTER_READLIST } from '../modules/local/filter_readlist/main'

workflow READ_DECONTAMINATION {

    take:
    qc_pass_reads      // channel: [ val(meta), path(reads) ]

    main:
    ch_versions = Channel.empty()
    ch_candidate_reads = Channel.empty()

    if (params.READ_DECONTAMINATION.host_assembly) {
       MINIMAP2_ALIGN(qc_pass_reads, [[],params.READ_DECONTAMINATION.host_assembly], true, false, false)
       ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

       FILTER_MAPPED(MINIMAP2_ALIGN.out.bam)
       ch_versions = ch_versions.mix(FILTER_MAPPED.out.versions)

       ch_candidate_reads =  qc_pass_reads.join(FILTER_MAPPED.out.read_list.map{meta -> meta = [meta[0],meta[1],meta[1]]})
    }

    if (params.READ_DECONTAMINATION.host_krakendb) {
       KRAKEN2_HOST(qc_pass_reads, params.READ_DECONTAMINATION.host_krakendb,false,true)
       ch_versions = ch_versions.mix(KRAKEN2_HOST.out.versions)

       ch_candidate_reads =  qc_pass_reads.join(KRAKEN2_HOST.out.classified_reads_assignment.map{meta -> meta = [meta[0],meta[1], meta[1]]})
    }

    if (params.READ_DECONTAMINATION.host_assembly) {
       if (params.READ_DECONTAMINATION.host_krakendb) {
          ch_decom_reads = FILTER_MAPPED.out.read_list.join(KRAKEN2_HOST.out.classified_reads_assignment, by: [0])
          ch_candidate_reads = qc_pass_reads.join(ch_decom_reads)
       }
    }

    FILTER_READLIST(ch_candidate_reads)
    //versions

    SEQKIT_GREP(FILTER_READLIST.out.host_readlist)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    NANOPLOT_POSTQC(SEQKIT_GREP.out.filter)
    ch_versions = ch_versions.mix(NANOPLOT_POSTQC.out.versions)

    SEQTK_FQCHK_POSTQC(SEQKIT_GREP.out.filter, params.SEQTK_FQCHK.endseq_len)
    ch_versions = ch_versions.mix(SEQTK_FQCHK_POSTQC.out.versions)

    emit:
    host_readlist = FILTER_READLIST.out.host_readlist.map{meta -> meta = [meta[0], meta[2]]}
    postqc_reads  = SEQKIT_GREP.out.filter.map{meta -> meta = [meta[0], meta[1]] }
    postqc_results = SEQTK_FQCHK_POSTQC.out.pbq_se.join(NANOPLOT_POSTQC.out.qc_input, by: [0])
    versions = ch_versions
}

