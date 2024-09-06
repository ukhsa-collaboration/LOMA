/*
 * Perform antimicrobial resistance typing and summarize with hAMRonization
 */

include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/run/main'
include { ABRICATE_RUN } from '../modules/nf-core/abricate/run/main'
include { RESFINDER } from '../modules/local/resfinder/run/main'
//include { ABRITAMR_RUN } from '../modules/nf-core/abritamr/run/main'                                                                                                                                              
include { RGI_MAIN } from '../modules/local/rgi/main/main'                                                                                                                                                      
include { HAMRONIZATION_ABRICATE } from '../modules/local/hamronization/abricate/main'                                                                                                                          
include { HAMRONIZATION_RGI } from '../modules/local/hamronization/rgi/main'                                                                                                                                    
include { HAMRONIZATION_AMRFINDERPLUS } from '../modules/local/hamronization/amrfinderplus/main'                                                                                                                
include { HAMRONIZATION_RESFINDER } from '../modules/local/hamronization/resfinder/main'
include { HAMRONIZATION_POINTFINDER } from '../modules/local/hamronization/pointfinder/main'
include { HAMRONIZATION_SUMMARIZE } from '../modules/local/hamronization/summarize/main'                                                                                                                        


workflow AMR_TYPING {

    take:
    ch_assembly

    main:
    ch_versions = Channel.empty()
    ch_amr_summarize = Channel.empty()

    if (params.AMRFINDERPLUS_RUN.db) {
        AMRFINDERPLUS_RUN(ch_assembly, params.AMRFINDERPLUS_RUN.db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions) 

        HAMRONIZATION_AMRFINDERPLUS(AMRFINDERPLUS_RUN.out.report, "tsv", params.AMRFINDERPLUS_RUN.software_version, params.AMRFINDERPLUS_RUN.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_AMRFINDERPLUS.out.versions)

        ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_AMRFINDERPLUS.out.tsv)
    }


    ABRICATE_RUN(ch_assembly)
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)

    HAMRONIZATION_ABRICATE(ABRICATE_RUN.out.report,"tsv", params.ABRICATE_RUN.software_version, params.ABRICATE_RUN.db_version)
    ch_versions = ch_versions.mix(HAMRONIZATION_ABRICATE.out.versions)

    ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_ABRICATE.out.tsv)

    RGI_MAIN(ch_assembly)
    ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

    HAMRONIZATION_RGI(RGI_MAIN.out.tsv, "tsv", params.RGI_MAIN.software_version, params.RGI_MAIN.db_version)
    ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions)

    ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_RGI.out.tsv)

    RESFINDER(ch_assembly)
    ch_versions = ch_versions.mix(RESFINDER.out.versions)

    HAMRONIZATION_RESFINDER(RESFINDER.out.results_resfinder,"tsv", params.RESFINDER.db_version, params.RESFINDER.db_version)
    ch_versions = ch_versions.mix(HAMRONIZATION_RESFINDER.out.versions)

    ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_RESFINDER.out.tsv)

    HAMRONIZATION_POINTFINDER(RESFINDER.out.results_pointfinder, "tsv", params.POINTFINDER.db_version, params.POINTFINDER.db_version)
    ch_versions = ch_versions.mix(HAMRONIZATION_POINTFINDER.out.versions)

    ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_POINTFINDER.out.tsv)
    ch_amr_summarize = ch_amr_summarize.groupTuple(by: [0])

    HAMRONIZATION_SUMMARIZE(ch_amr_summarize, ['tsv', 'json', 'interactive'])
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE.out.versions)

    ch_amr_summarize = ch_amr_summarize.mix(HAMRONIZATION_SUMMARIZE.out.tsv)


    emit:
    ch_amr_reports = ch_amr_summarize
    versions = ch_versions
}
