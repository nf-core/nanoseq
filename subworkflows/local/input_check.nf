/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    input_path

    main:
    /*
     * Check samplesheet is valid
     */
    SAMPLESHEET_CHECK ( samplesheet, input_path )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { ch_sample }

    emit:
    ch_sample // [ meta, input_file, barcode, nanopolish_fast5 ]
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str, nanopolish_fast5 ]
def get_sample_info(LinkedHashMap sample) {
    def meta = [:]
    meta.id  = sample.sample

    // Resolve fasta and gtf file if using iGenomes
//    def fasta = false
//    def gtf   = false
//    if (sample.fasta) {
//        if (genomeMap && genomeMap.containsKey(sample.fasta)) {
//            fasta = file(genomeMap[sample.fasta].fasta, checkIfExists: true)
//            gtf   = file(genomeMap[sample.fasta].gtf, checkIfExists: true)
//        } else {
//            fasta = file(sample.fasta, checkIfExists: true)
//        }
//    }

    // Check if input file and gtf file exists
    input_file = sample.input_file ? file(sample.input_file, checkIfExists: true) : null
//    gtf        = sample.gtf        ? file(sample.gtf, checkIfExists: true)        : gtf

    return [ meta, input_file, sample.barcode, sample.nanopolish_fast5 ]
}

