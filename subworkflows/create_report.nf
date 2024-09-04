/*
 * 
 */

include { SUMMARIZE_RESULTS } from '../modules/local/summarize_results/main'
include { SUMMARIZE_AMR } from '../modules/local/summarize_amr/main'

workflow CREATE_REPORT {

    take:
    ch_final_report    // channel: [ val(meta), path(files) ]

    main:
    // 
    ch_versions = Channel.empty()

    ch_final_report
       .map{meta -> meta = [[id: meta[0].id, run_id: meta[0].run_id, barcode:meta[0].barcode, target: meta[0].target, single_end: meta[0].single_end], meta[1]]}
       .transpose()
       .groupTuple(by: [0])
       .set{ch_final_report_merged}


    SUMMARIZE_RESULTS(ch_final_report_merged, params.CREATE_REPORT.summary_template)
    ch_versions = ch_versions.mix(SUMMARIZE_RESULTS.out.versions)

    SUMMARIZE_AMR(ch_final_report_merged, params.CREATE_REPORT.amr_template)
    ch_versions = ch_versions.mix(SUMMARIZE_AMR.out.versions)


    emit:
    
    versions = ch_versions
}
