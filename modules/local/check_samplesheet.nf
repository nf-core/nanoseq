// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str ]
def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {

    // Resolve fasta and gtf file if using iGenomes
    def fasta = false
    def gtf   = false
    if (sample.genome) {
        if (genomeMap && genomeMap.containsKey(sample.genome)) {
            fasta = file(genomeMap[sample.genome].fasta, checkIfExists: true)
            gtf   = file(genomeMap[sample.genome].gtf, checkIfExists: true)
        } else {
            fasta = file(sample.genome, checkIfExists: true)
        }
    }

    // Check if input file and gtf file exists
    input_file = sample.input_file ? file(sample.input_file, checkIfExists: true) : null
    gtf        = sample.gtf        ? file(sample.gtf, checkIfExists: true)        : gtf

    return [ sample.sample, input_file, sample.barcode, fasta, gtf, sample.is_transcripts.toBoolean(), fasta.toString()+';'+gtf.toString() ]
}
