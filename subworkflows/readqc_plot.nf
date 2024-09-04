/*
 * Plot read quality control metrics
 */

include { PLOT_READQCMETRICS } from '../modules/local/plot_readqcmetrics/main'

workflow READQC_PLOT {

    take:
    ch_metrics

    main:
    ch_versions = Channel.empty()

    PLOT_READQCMETRICS(ch_metrics, params.PLOT_READQCMETRICS.template)
    ch_versions = ch_versions.mix(PLOT_READQCMETRICS.out.versions)

    emit:
    ch_readqc_stats = PLOT_READQCMETRICS.out.nanostats
    versions = ch_versions
}

