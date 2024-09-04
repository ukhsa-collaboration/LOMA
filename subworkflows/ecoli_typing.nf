/*
 * Type Escherichia coli genomes
 */

include { ECTYPER } from '../modules/nf-core/ectyper/main'                                                                                                                                                        
include { STECFINDER } from '../modules/nf-core/stecfinder/main'                                                                                                                                                  
include { SHIGEIFINDER } from '../modules/nf-core/shigeifinder/main'                                                                                                                                              
include { SHIGATYPER } from '../modules/nf-core/shigatyper/main'                                                                                                                                                  
include { MYKROBE_PREDICT as MYKROBE_PREDICT_SSONNEI } from '../modules/nf-core/mykrobe/predict/main'                                                                                                                                        
include { FILTER_ECOLI } from '../modules/local/filter_ecoli/main'
include { FILTER_TYPING as FILTER_TYPING_ECOLI } from '../modules/local/filter_typing/main'

workflow ECOLI_TYPING {

    take:
    ch_merged    // channel: [ val(meta), path(assembly), path(reads) ]
    ch_mlst    // channel: [ val(meta), path(mlst_results), path(krocus_results) ]

    main:
    // Runs typing method specific to Escherichia coli/Shigella

    ch_versions = Channel.empty()
    ch_mykrobe_out = Channel.empty()

    // Select only assembly or reads depending on requirements
    ch_assembly = ch_merged.map{meta -> meta = [meta[0], meta[1]]}
    ch_reads = ch_merged.map{meta -> meta = [meta[0], meta[2]]}

    // Serotype E. coli
    ECTYPER(ch_assembly)
    ch_versions = ch_versions.mix(ECTYPER.out.versions)

    // Identify serotype, cluster specific genes, O/H antigen genes
    STECFINDER(ch_assembly)
    ch_versions = ch_versions.mix(STECFINDER.out.versions)

    // Differentiate Shigella and Enteroinvasive E. coli
    SHIGEIFINDER(ch_assembly)
    ch_versions = ch_versions.mix(SHIGEIFINDER.out.versions)

    // Determine Shigella serotype
    SHIGATYPER(ch_reads)
    ch_versions = ch_versions.mix(SHIGATYPER.out.versions)

    ch_prefilt = ch_merged.join(ch_mlst, by: [0]).map{meta -> meta = [meta[0], meta[1], meta[3], meta[4]]}

    FILTER_TYPING_ECOLI(ch_prefilt, "ecoli")
    ch_versions = ch_versions.mix(FILTER_TYPING_ECOLI.out.versions)

    ch_mykrobe_in_ecoli = FILTER_TYPING_ECOLI.out.typing.map{meta -> meta = [meta[0], [type: (meta[2].getName().split(/\./)[2])], meta[1]]}
       .filter({meta, type, seq -> type.type.toLowerCase().matches("cc152")})
       .map{meta -> meta = [meta[0], meta[2]]}

    // Shigella sonnei specific genotyping and AMR prediction
    MYKROBE_PREDICT_SSONNEI(ch_mykrobe_in_ecoli, "sonnei")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT_SSONNEI.out.versions)

    // Merge outputs together
    ch_in_ecoli_4 = SHIGEIFINDER.out.tsv
       .join(SHIGATYPER.out.tsv, by: [0])
       .join(ECTYPER.out.tsv, by: [0])
       .join(STECFINDER.out.tsv, by: [0])
       .join(ch_mlst, by: [0])

    ch_in_ecoli_4 = ch_in_ecoli_4.mix(MYKROBE_PREDICT_SSONNEI.out.csv)

    ch_in_ecoli_4.filter(meta -> meta.toString().contains('mykrobe')).set{ch_in_ecoli_5}
    ch_in_ecoli_4.filter(meta -> !meta.toString().contains('mykrobe')).map{meta -> meta = [meta[0], meta[1], meta[2], meta[3], meta[4], meta[5],meta[6], []]}.set{ch_in_ecoli_5}

    // Summarize typing results, try and infer likely pathotype
    FILTER_ECOLI(ch_in_ecoli_5)
    ch_versions = ch_versions.mix(FILTER_ECOLI.out.versions)


    emit:
    ch_ecoli_report = FILTER_ECOLI.out.typing
    versions = ch_versions
}
