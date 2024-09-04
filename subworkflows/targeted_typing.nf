/*
 * Use various tools that only type specific subsets of species (for which databases exist)
 */

include { PLASMIDFINDER } from '../modules/nf-core/plasmidfinder/main'                                                                                                                                            
include { VIRULENCEFINDER } from '../modules/local/virulencefinder/main'
include { BLAST_BLASTN } from '../modules/local/blast/blastn/main'
include { FILTER_BLAST } from '../modules/local/filter_blast/main'
include { GENEFINDER } from '../modules/local/genefinder/main'
include { SUMMARIZE_GENEFINDING } from '../modules/local/summarize_genefinding/main'

workflow TARGETED_TYPING {

    take:
    ch_merged    // channel: [ val(meta), path(assembly), path(reads) ]

    main:
    // Perform typing on specific combinations of species (those with relevant databases) and generally type plasmids

    ch_versions = Channel.empty()
    ch_targetedtyping_reports = Channel.empty()

    ch_assembly = ch_merged.map{meta -> meta = [meta[0],meta[1]]}

    // Subset target species for VirulenceFinder
    ch_secondary_remainder = ch_assembly.filter({meta, fasta -> (meta.clean_id.matches("Listeria [a-z]*|Enterococcus [a-z]*|Escherichia coli|Staphylococcus aureus"))})

    // Identify and type plasmid replicons
    PLASMIDFINDER(ch_assembly)
    ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions)

    ch_targetedtyping_reports = ch_targetedtyping_reports.mix(PLASMIDFINDER.out.tsv)

    ch_vf_out = ch_assembly.map{meta -> meta = [meta[0], []]}

    // If VirulenceFinder database is provided, then run on target species
    if (params.VIRULENCEFINDER.db) {
        VIRULENCEFINDER(ch_secondary_remainder, params.VIRULENCEFINDER.db) 
        ch_vf_out = VIRULENCEFINDER.out.tsv
    } 

    // Subset species with gene databases (if provided), select assemblies or assemblies and reads
//    ch_merged.filter({meta, fasta, fastq -> !(meta.gene_DB.matches(""))}).map{meta -> meta = [meta[0], meta[1]]}.set{ch_secondary_remainder_1B}
//    ch_merged.filter({meta, fasta, fastq -> !(meta.gene_DB.matches(""))}).set{ch_secondary_remainder_1C}

//    // If gene databases provided, run BLAST (on assemblies) and Genefinder (on reads)
//    if (params.TARGETED_TYPING.genedbdir) {
//        BLAST_BLASTN(ch_secondary_remainder_1B, params.TARGETED_TYPING.genedbdir)
//        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

//        FILTER_BLAST(BLAST_BLASTN.out.txt)
//        ch_versions = ch_versions.mix(FILTER_BLAST.out.versions)

//        GENEFINDER(ch_secondary_remainder_1C, params.TARGETED_TYPING.genedbdir)
//        ch_versions = ch_versions.mix(GENEFINDER.out.versions)

//        ch_summarize_in_1 = FILTER_BLAST.out.filtered_hits.join(GENEFINDER.out.results, by: [0])
//        ch_summarize_in_2 = ch_summarize_in_1.join(ch_vf_out, by: [0])

//        SUMMARIZE_GENEFINDING(ch_summarize_in_2)
//        ch_versions = ch_versions.mix(SUMMARIZE_GENEFINDING.out.versions)

//        ch_targetedtyping_reports = ch_targetedtyping_reports.mix(SUMMARIZE_GENEFINDING.out.merged_results)
//        ch_targetedtyping_reports = ch_targetedtyping_reports.mix(SUMMARIZE_GENEFINDING.out.summary)

//    }

    emit:
    ch_targetedtyping_reports
    versions = ch_versions
}
