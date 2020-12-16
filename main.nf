#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/nanoseq
========================================================================================
 nf-core/nanoseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/nanoseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$baseDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/nanoseq --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /*
     * SUBWORKFLOW: Run main nf-core/nanoseq analysis pipeline
     */
    include { NANOSEQ } from './nanoseq' addParams( summary_params: summary_params )
    NANOSEQ ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
