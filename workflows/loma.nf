/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RUN_QC } from '../subworkflows/run_qc.nf'
include { READ_QC } from '../subworkflows/read_qc.nf'
include { READ_DECONTAMINATION } from '../subworkflows/read_decontamination.nf'
include { ASSEMBLY } from '../subworkflows/assembly.nf'
include { BIN_ASSIGNMENT } from '../subworkflows/bin_assignment.nf'
include { TAXONOMIC_PROFILING } from '../subworkflows/taxonomic_profiling.nf'
include { BIN_QC } from '../subworkflows/bin_qc.nf'
include { CONTIG_QC } from '../subworkflows/contig_qc.nf'
include { BIN_TAXONOMY } from '../subworkflows/bin_taxonomy.nf'
include { PROKARYA_TYPING } from '../subworkflows/prokarya_typing.nf'
include { READQC_PLOT } from '../subworkflows/readqc_plot.nf'

include { CREATE_REPORT } from '../subworkflows/create_report.nf'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOMA {

    ch_versions = Channel.empty()
    ch_final_report = Channel.empty()

if (params.input_type =="guppy") {
    Channel
     .fromPath(params.input)
     .ifEmpty {exit 1, log.info "Cannot find input file"}
     .splitCsv(sep: "\t", header: ['run_id', 'barcode', 'id', 'target'])
     .filter { row -> "$row.target" != 'NONE'}
     .filter { row -> if ((row.id.contains('.'))) {throw new IllegalArgumentException("Value '$row.id' contains illegal character: '.'")} else {return true} }
     .map { row -> meta = [id: row.id, run_id: row.run_id, barcode: row.barcode, target: row.target, single_end: true]}
     .set {scheme_all}

    RUN_QC(scheme_all)
    ch_versions = ch_versions.mix(RUN_QC.out.versions)
    ch_reads_pre = RUN_QC.out.aggregated_reads
}

else if (params.input_type =="fastq") {
    Channel
     .fromPath(params.input)
     .ifEmpty{ throw new IllegalArgumentException("No input rows in file: $params.input") }
     .splitCsv(sep: "\t", header: ['run_id', 'barcode', 'id', 'target', 'fastq'])
     .filter { row -> "$row.target" != 'NONE'}
     .ifEmpty{ throw new IllegalArgumentException("No input rows in file: $params.input") }
     .filter { row -> if (row.id == null || row.id == "") {throw new IllegalArgumentException("Missing or null value of 'id' for input: \n\nid:\t\t$row.id\nrow_id:\t\t\t$row.run_id\nsample_type:\t\t$row.target\nfastq:\t\t\t$row.fastq\n")} else {return true} }
     .filter { row -> if (row.run_id == null || row.run_id == "") {throw new IllegalArgumentException("Missing or null value of 'run_id' for input: \n\nid:\t\t$row.id\nrow_id:\t\t\t$row.run_id\ntarget:\t\t$row.target\nfastq:\t\t\t$row.fastq\n")} else {return true} }
     .filter { row -> if (row.fastq == null || row.fastq == "") {throw new IllegalArgumentException("Missing or null value of 'fastq' for input: \n\nid:\t\t$row.id\nrow_id:\t\t\t$row.run_id\ntarget:\t\t$row.target\nfastq:\t\t\t$row.fastq\n")} else {return true} }
     .map { row -> meta = [[id: row.id.replaceAll('\\.','_'), run_id: row.run_id, barcode: (row.barcode ?: 'NA'), target: (row.target ?: 'NA'), single_end: true], file(row.fastq, checkIfExists: true)]}
     .set {ch_reads_pre}
}

    READ_QC(ch_reads_pre)
    ch_versions = ch_versions.mix(READ_QC.out.versions)

    READ_DECONTAMINATION(READ_QC.out.qc_pass_reads)
    ch_versions = ch_versions.mix(READ_DECONTAMINATION.out.versions)

    ch_final_report = ch_final_report.mix(READ_DECONTAMINATION.out.host_readlist)

    ch_metrics = READ_QC.out.preqc_results.join(READ_DECONTAMINATION.out.postqc_results, by: [0])

    READQC_PLOT(ch_metrics)
    ch_versions = ch_versions.mix(READQC_PLOT.out.versions)

    ch_final_report = ch_final_report.mix(READQC_PLOT.out.ch_readqc_stats)

    if (!params.skip_taxonomic_profiling ) {
        if (params.TAXONOMIC_PROFILING.use_decontaminated_reads) {
            ch_taxonomic_profiling_reads = READ_DECONTAMINATION.out.postqc_reads
        } else {
            ch_taxonomic_profiling_reads = READ_QC.out.qc_pass_reads
        }

        TAXONOMIC_PROFILING(ch_taxonomic_profiling_reads)
        ch_versions = ch_versions.mix(TAXONOMIC_PROFILING.out.versions)

        ch_final_report = ch_final_report.mix(TAXONOMIC_PROFILING.out.ch_taxreports)
    }
 
    if (!params.skip_assembly ) {
        ASSEMBLY(READ_DECONTAMINATION.out.postqc_reads)
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        BIN_ASSIGNMENT(ASSEMBLY.out.final_assembly_mapped)
        ch_versions = ch_versions.mix(BIN_ASSIGNMENT.out.versions)

        CONTIG_QC(BIN_ASSIGNMENT.out.reads_complete_assembly)
        ch_versions = ch_versions.mix(CONTIG_QC.out.versions)

        BIN_QC(BIN_ASSIGNMENT.out.assembled_bins, CONTIG_QC.out.ch_plt_p1)
        ch_versions = ch_versions.mix(BIN_QC.out.versions)

        ch_final_report = ch_final_report.mix(BIN_QC.out.ch_binreports)

        BIN_TAXONOMY(BIN_QC.out.ch_bin_taxonomy, CONTIG_QC.out.complete_assembly_bam)
        ch_versions = ch_versions.mix(BIN_TAXONOMY.out.versions)

        if (!params.skip_bacterial_typing ) {
            PROKARYA_TYPING(BIN_TAXONOMY.out.ch_prokarya_reads)
            ch_versions = ch_versions.mix(PROKARYA_TYPING.out.versions)

            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_amr_reports)
            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_ecoli_report)
            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_salmonella_report)
            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_lmonocytogenes_report)
            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_targetedtyping_reports)
            ch_final_report = ch_final_report.mix(PROKARYA_TYPING.out.ch_sequencetyping_reports)
        }
    }

    CREATE_REPORT(ch_final_report)
    ch_versions = ch_versions.mix(CREATE_REPORT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

}
