#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nxf/loma
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://gitlab.phe.gov.uk/Duncan.Berger2/loma
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def String command = "./run_loma --input <input_file> \n\n  ./run_loma --input <input_file> --input_type guppy --run_dir <guppy output directory>"
//    helpMessage()
//    log.info paramsHelp(command) + NfcoreTemplate.dashedLine(params.monochrome_logs)
    log.info paramsHelp(command)
    exit 0
}

if (!params.input) {
    def String command = "\nMissing parameter: --input\n\nUsage:\n\tnextflow run main.nf --input <input_file> --input_type fastq\n\n\tnextflow run main.nf --input <input_file> --input_type guppy --run_dir <guppy output directory>\n"
    log.info paramsHelp(command)
    exit 0
}
//if (!params.input_type) {
//    def String command = "\nMissing parameter: --input_type\n\nUsage:\n\tnextflow run main.nf --input <input_file> --input_type fastq\n\n\tnextflow run main.nf --input <input_file> --input_type guppy --run_dir <guppy output directory>\n"
//    log.info paramsHelp(command)
//    exit 0
//}
if (params.input_type == "guppy") {
    if (!params.run_dir) {
        def String command = "\nMissing parameter: --run_dir\n\nUsage:\n\tnextflow run main.nf --input <input_file> --input_type fastq\n\n\tnextflow run main.nf --input <input_file> --input_type guppy --run_dir <guppy output directory>\n"
        log.info paramsHelp(command)
        exit 0
    }
}




// Validate input parameters
//if (params.validate_params) {
//    validateParameters()
//}

//validateParameters()

//WorkflowMain.initialise(workflow, params, log)

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LOMA } from './workflows/loma'

//
// WORKFLOW: Run main nxf/loma analysis pipeline
//
workflow NXF_LOMA {
    LOMA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    NXF_LOMA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
