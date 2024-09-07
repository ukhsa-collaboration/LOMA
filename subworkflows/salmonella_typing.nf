/*
 * Type Salmonella genomes
 */

include { SEQSERO2 } from '../modules/nf-core/seqsero2/main'                                                                                                                                                      
include { SISTR } from '../modules/nf-core/sistr/main'                                                                                                                                                            
include { MYKROBE_PREDICT as MYKROBE_PREDICT_SPARATYPHIB } from '../modules/nf-core/mykrobe/predict/main'
include { MYKROBE_PREDICT as MYKROBE_PREDICT_STYPHI } from '../modules/nf-core/mykrobe/predict/main'
include { FILTER_TYPING as FILTER_TYPING_SALMONELLA } from '../modules/local/filter_typing/main'
include { FILTER_SALMONELLA } from '../modules/local/filter_salmonella/main'

workflow SALMONELLA_TYPING {

    take:
    ch_salmonella    // channel: [ val(meta), path(assembly) ]
    ch_mlst    // channel: [ val(meta), path(mlst), path(krocus) ]

    main:
    // Runs typing method specific to Salmonella

    ch_versions = Channel.empty()
    ch_mykrobe_out = Channel.empty()

    // Serotype prediction
    SEQSERO2(ch_salmonella)
    ch_versions = ch_versions.mix(SEQSERO2.out.versions)

    // Serovar prediction, identification of antigen and cgMLST genes 
    SISTR(ch_salmonella)
    ch_versions = ch_versions.mix(SISTR.out.versions)

    ch_prefilt_salmonella1 = ch_salmonella.join(SEQSERO2.out.tsv, by: [0])
    ch_prefilt_salmonella2 = ch_prefilt_salmonella1.join(SISTR.out.tsv, by: [0])

    // Filter Seqsero2 and SISTR results to identify Paratyphi B genomes
    FILTER_TYPING_SALMONELLA(ch_prefilt_salmonella2, "salmonella")
    ch_versions = ch_versions.mix(FILTER_TYPING_SALMONELLA.out.versions)

    // Separate Paratyphi B and all other assemblies for Mykrobe
    ch_mykrobe_in_typhi = FILTER_TYPING_SALMONELLA.out.typing.map{meta -> meta = [meta[0], [type: (meta[2].getName().split(/\./)[2])], meta[1]]}.filter({meta, type, seq -> !(type.type.toLowerCase().matches("paratyphi_b"))}).map{meta -> meta = [meta[0], meta[2]]}
    ch_mykrobe_in_paratyphib = FILTER_TYPING_SALMONELLA.out.typing.map{meta -> meta = [meta[0], [type: (meta[2].getName().split(/\./)[2])], meta[1]]}.filter({meta, type, seq -> type.type.toLowerCase().matches("paratyphi_b")}).map{meta -> meta = [meta[0], meta[2]]}

    // Runs Mykrobe on all non-Paratyphi B genomes
    MYKROBE_PREDICT_STYPHI(ch_mykrobe_in_typhi,"typhi")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT_STYPHI.out.versions)

    // Runs Mykrobe on all Paratyphi B genomes
    MYKROBE_PREDICT_SPARATYPHIB(ch_mykrobe_in_paratyphib,"paratyphiB")
    ch_versions = ch_versions.mix(MYKROBE_PREDICT_SPARATYPHIB.out.versions)

//    ch_mykrobe_out = MYKROBE_PREDICT_SPARATYPHIB.out.csv
//    ch_mykrobe_in_paratyphib.ifEmpty(MYKROBE_PREDICT_STYPHI.out.csv)

    ch_mykrobe_out = ch_mykrobe_out.mix(MYKROBE_PREDICT_STYPHI.out.csv)
    ch_mykrobe_out = ch_mykrobe_out.mix(MYKROBE_PREDICT_SPARATYPHIB.out.csv)

    // Merge outputs together    
    ch_in_salmonella_1 = SEQSERO2.out.tsv.join(SISTR.out.tsv, by: [0])
    ch_in_salmonella_2 = ch_in_salmonella_1.join(ch_mlst, by: [0])
    ch_in_salmonella_3 = ch_in_salmonella_2.join(ch_mykrobe_out, by: [0])

    // Summarize typing results
    FILTER_SALMONELLA(ch_in_salmonella_3)
    ch_versions = ch_versions.mix(FILTER_SALMONELLA.out.versions)

    emit:

    ch_salmonella_report = FILTER_SALMONELLA.out.typing
    versions = ch_versions
}
