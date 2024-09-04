/*
 * Perform per metagenome assembled genome (MAG) quality control and typing
 */

include { KRAKEN2_KRAKEN2 as KRAKEN2_TAXPROFILING } from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN_KRAKEN2 } from '../modules/nf-core/bracken/bracken/main'
include { BRACKEN_BRACKEN as BRACKEN_CENTRIFUGER } from '../modules/nf-core/bracken/bracken/main'
include { SYLPH_SKETCH } from '../modules/local/sylph/sketch/main'
include { SYLPH_PROFILE } from '../modules/local/sylph/profile/main'
include { TAXPASTA_STANDARDISE as TAXPASTA_STANDARDISE_KRAKEN2 } from '../modules/nf-core/taxpasta/standardise/main'
include { TAXPASTA_STANDARDISE as TAXPASTA_STANDARDISE_CENTRIFUGER } from '../modules/nf-core/taxpasta/standardise/main'
include { CENTRIFUGER_CENTRIFUGER } from '../modules/local/centrifuger/centrifuger/main'
include { CENTRIFUGER_KREPORT } from '../modules/local/centrifuger/kreport/main'
include { PLOT_TAXHITS } from '../modules/local/plot_taxhits/main' 
include { PLOT_INDVTAXHITS as PLOT_KRAKEN2 } from '../modules/local/plot_indvtaxhits/main'
include { PLOT_INDVTAXHITS as PLOT_CENTRIFUGER } from '../modules/local/plot_indvtaxhits/main'
include { PLOT_INDVTAXHITS as PLOT_KRAKEN2BRACKEN } from '../modules/local/plot_indvtaxhits/main'
include { PLOT_INDVTAXHITS as PLOT_CENTRIFUGERBRACKEN } from '../modules/local/plot_indvtaxhits/main'
include { PLOT_INDVTAXHITS as PLOT_SYLPH } from '../modules/local/plot_indvtaxhits/main'
include { PARSE_TAXHITS as PARSE_KRAKEN2HITS } from '../modules/local/parse_taxhits/main'
include { PARSE_TAXHITS as PARSE_CENTRIFUGERHITS } from '../modules/local/parse_taxhits/main'
include { PARSE_TAXHITS as PARSE_SYLPHHITS } from '../modules/local/parse_taxhits/main'

workflow TAXONOMIC_PROFILING {

    take:
    postqc_reads

    main:
    ch_versions = Channel.empty()
    ch_parsedreports = Channel.empty()
    ch_parsedreports_gt = Channel.empty()


    if (params.TAXONOMIC_PROFILING.krakendb ) {
        KRAKEN2_TAXPROFILING(postqc_reads, params.TAXONOMIC_PROFILING.krakendb, false, false)
        ch_versions = ch_versions.mix(KRAKEN2_TAXPROFILING.out.versions)

        BRACKEN_KRAKEN2(KRAKEN2_TAXPROFILING.out.report, params.TAXONOMIC_PROFILING.krakendb)
        ch_versions = ch_versions.mix(BRACKEN_KRAKEN2.out.versions)

        PLOT_KRAKEN2BRACKEN(BRACKEN_KRAKEN2.out.reports, "Kraken2, Bracken", [], params.TAXONOMIC_PROFILING.template)
        ch_versions = ch_versions.mix(PLOT_KRAKEN2BRACKEN.out.versions)

        if (params.TAXONOMIC_PROFILING.dbdir) {
            TAXPASTA_STANDARDISE_KRAKEN2(KRAKEN2_TAXPROFILING.out.report, params.TAXONOMIC_PROFILING.dbdir)
            ch_versions = ch_versions.mix(TAXPASTA_STANDARDISE_KRAKEN2.out.versions)

            PLOT_KRAKEN2(TAXPASTA_STANDARDISE_KRAKEN2.out.standardised_profile, "Kraken2, Taxpasta", [], params.TAXONOMIC_PROFILING.template)
            ch_versions = ch_versions.mix(PLOT_KRAKEN2.out.versions)

            ch_parsedreports = ch_parsedreports.mix(BRACKEN_KRAKEN2.out.reports)

            if (params.TAXONOMIC_PROFILING.target_species) {
                PARSE_KRAKEN2HITS(BRACKEN_KRAKEN2.out.reports, params.TAXONOMIC_PROFILING.target_species, [params.PARSE_KRAKEN2HITS.min_target_reads, params.PARSE_KRAKEN2HITS.min_target_fraction], "Bracken")
                ch_versions = ch_versions.mix(PARSE_KRAKEN2HITS.out.versions)

                ch_parsedreports = ch_parsedreports.mix(PARSE_KRAKEN2HITS.out.targets_filtered)
            }
        }
    }

    if (params.TAXONOMIC_PROFILING.sylphdb ) {
        SYLPH_SKETCH(postqc_reads)
        ch_versions = ch_versions.mix(SYLPH_SKETCH.out.versions)

        SYLPH_PROFILE(SYLPH_SKETCH.out.sketch, params.TAXONOMIC_PROFILING.sylphdb)
        ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions)

        if (params.TAXONOMIC_PROFILING.gtdb_metadata) {
            PLOT_SYLPH(SYLPH_PROFILE.out.results, "Sylph", params.TAXONOMIC_PROFILING.gtdb_metadata, params.TAXONOMIC_PROFILING.template)
            ch_versions = ch_versions.mix(PLOT_SYLPH.out.versions)

            PARSE_SYLPHHITS(SYLPH_PROFILE.out.results, params.TAXONOMIC_PROFILING.target_species, [params.PARSE_SYLPHHITS.min_target_reads, params.PARSE_SYLPHHITS.min_target_fraction], "Sylph")
            ch_versions = ch_versions.mix(PARSE_SYLPHHITS.out.versions)

            ch_parsedreports = ch_parsedreports.mix(SYLPH_PROFILE.out.results)
            ch_parsedreports = ch_parsedreports.mix(PARSE_SYLPHHITS.out.targets_filtered)

        }
    }

    if (params.TAXONOMIC_PROFILING.centrifugerdb ) {
        CENTRIFUGER_CENTRIFUGER(postqc_reads,params.TAXONOMIC_PROFILING.centrifugerdb)
        ch_versions = ch_versions.mix(CENTRIFUGER_CENTRIFUGER.out.versions)

        CENTRIFUGER_KREPORT(CENTRIFUGER_CENTRIFUGER.out.results, params.TAXONOMIC_PROFILING.centrifugerdb)
        ch_versions = ch_versions.mix(CENTRIFUGER_KREPORT.out.versions)

        ch_parsedreports = ch_parsedreports.mix(CENTRIFUGER_KREPORT.out.kreport)

        if (params.TAXONOMIC_PROFILING.krakendb ) {
            BRACKEN_CENTRIFUGER(CENTRIFUGER_KREPORT.out.kreport, params.TAXONOMIC_PROFILING.krakendb)
            ch_versions = ch_versions.mix(BRACKEN_CENTRIFUGER.out.versions)

            PLOT_CENTRIFUGERBRACKEN(BRACKEN_CENTRIFUGER.out.reports, "Centrifuger, Bracken", [], params.TAXONOMIC_PROFILING.template)
            ch_versions = ch_versions.mix(PLOT_CENTRIFUGERBRACKEN.out.versions)

            ch_parsedreports = ch_parsedreports.mix(BRACKEN_CENTRIFUGER.out.reports)

            if (params.TAXONOMIC_PROFILING.target_species) {
                PARSE_CENTRIFUGERHITS(BRACKEN_CENTRIFUGER.out.reports, params.TAXONOMIC_PROFILING.target_species, [params.PARSE_CENTRIFUGERHITS.min_target_reads, params.PARSE_CENTRIFUGERHITS.min_target_fraction], "Bracken")
                ch_versions = ch_versions.mix(PARSE_CENTRIFUGERHITS.out.versions)

                ch_parsedreports = ch_parsedreports.mix(PARSE_CENTRIFUGERHITS.out.targets_filtered)

            }
       }

       if (params.TAXONOMIC_PROFILING.dbdir) {
            TAXPASTA_STANDARDISE_CENTRIFUGER(CENTRIFUGER_KREPORT.out.kreport, params.TAXONOMIC_PROFILING.dbdir)
            ch_versions = ch_versions.mix(TAXPASTA_STANDARDISE_CENTRIFUGER.out.versions)

            PLOT_CENTRIFUGER(TAXPASTA_STANDARDISE_CENTRIFUGER.out.standardised_profile, "Centrifuger, Taxpasta", [], params.TAXONOMIC_PROFILING.template)
            ch_versions = ch_versions.mix(PLOT_CENTRIFUGER.out.versions)
        }
    }

    if (params.TAXONOMIC_PROFILING.krakendb ) {
        if (params.TAXONOMIC_PROFILING.sylphdb) {
            if (params.TAXONOMIC_PROFILING.centrifugerdb) {
                ch_in_1 = TAXPASTA_STANDARDISE_KRAKEN2.out.standardised_profile.join(TAXPASTA_STANDARDISE_CENTRIFUGER.out.standardised_profile, by: [0])
                ch_in_2 = ch_in_1.join(SYLPH_PROFILE.out.results, by: [0])
                ch_in_3 = ch_in_2.join(BRACKEN_KRAKEN2.out.reports, by: [0])
                ch_in_4 = ch_in_3.join(BRACKEN_CENTRIFUGER.out.reports, by: [0])

                PLOT_TAXHITS(ch_in_4, params.TAXONOMIC_PROFILING.gtdb_metadata, params.TAXONOMIC_PROFILING.template)
                ch_versions = ch_versions.mix(PLOT_TAXHITS.out.versions)
          }
       }
    } 


    emit:
    ch_taxreports = ch_parsedreports.groupTuple(by: [0])
    versions = ch_versions
}

