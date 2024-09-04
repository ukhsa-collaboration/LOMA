/*
 * Perform general typing and identify metagenome assembled genomes that should undergo species/genus specific typing
 */

include { AMR_TYPING } from '../subworkflows/amr_typing.nf'
include { SEQUENCE_TYPING } from '../subworkflows/sequence_typing.nf'
include { TARGETED_TYPING } from '../subworkflows/targeted_typing.nf'
include { SALMONELLA_TYPING } from '../subworkflows/salmonella_typing.nf'
include { LMONOCYTOGENES_TYPING } from '../subworkflows/lmonocytogenes_typing.nf'
include { ECOLI_TYPING } from '../subworkflows/ecoli_typing.nf'

workflow PROKARYA_TYPING {

    take:
    ch_merged    // channel: [ val(meta), path(assembly), path(reads) ]

    main:
    // Using taxonomic assignments made during the BIN_TAXONOMY step, pass relevant reads/assemblies to the appropriate subworkflows
    ch_versions = Channel.empty()

    ch_assembly = ch_merged.map{meta -> meta = [meta[0],meta[1]]}

    // Create channels for specific genera/species
    ch_salmonella = ch_merged.filter({meta, fasta, fastq -> meta.genus.toLowerCase().matches("^salmonella")}).map{meta -> meta = [meta[0],meta[1]]}
    ch_ecoli = ch_merged.filter({meta, fasta, fastq -> meta.clean_id.toLowerCase().matches("^escherichia coli")})
    ch_lmonocytogenes = ch_merged.filter({meta, fasta, fastq -> meta.clean_id.toLowerCase().matches("^listeria monocytogenes")}).map{meta -> meta = [meta[0],meta[1]]}

    // AMR typing steps
    AMR_TYPING(ch_assembly)
    ch_versions = ch_versions.mix(AMR_TYPING.out.versions)

    // Multi-locus sequence typing
    SEQUENCE_TYPING(ch_merged)
    ch_versions = ch_versions.mix(SEQUENCE_TYPING.out.versions)

    // Screen for specific genes/mobile genetic elements
    TARGETED_TYPING(ch_merged)
    ch_versions = ch_versions.mix(TARGETED_TYPING.out.versions)

    // E. coli/Shigella specific typing
    ECOLI_TYPING(ch_ecoli, SEQUENCE_TYPING.out.ch_mlst)
    ch_versions = ch_versions.mix(ECOLI_TYPING.out.versions)

    // Salmonella specific typing
    SALMONELLA_TYPING(ch_salmonella, SEQUENCE_TYPING.out.ch_mlst)
    ch_versions = ch_versions.mix(SALMONELLA_TYPING.out.versions)

    // Listeria monocytogenes specific typing
    LMONOCYTOGENES_TYPING(ch_lmonocytogenes)
    ch_versions = ch_versions.mix(LMONOCYTOGENES_TYPING.out.versions)

    emit:
    ch_amr_reports = AMR_TYPING.out.ch_amr_reports
    ch_targetedtyping_reports = TARGETED_TYPING.out.ch_targetedtyping_reports
    ch_sequencetyping_reports = SEQUENCE_TYPING.out.ch_mlst_reports
    ch_ecoli_report = ECOLI_TYPING.out.ch_ecoli_report
    ch_salmonella_report = SALMONELLA_TYPING.out.ch_salmonella_report
    ch_lmonocytogenes_report = LMONOCYTOGENES_TYPING.out.ch_lmonocytogenes_report
    versions = ch_versions
}
