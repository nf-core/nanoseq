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
    ch_sample // [[id:, barcode:, nanopolish_fast5:], [input_file]]
}

// Function to get list of [ meta, fastq ]
def get_sample_info(LinkedHashMap row) {
    def meta              = [:]
    meta.id               = row.sample
    meta.barcode          = row.barcode
    meta.nanopolish_fast5 = row.nanopolish_fast5
    input_file            = row.reads //? file(row.reads, checkIfExists: true) : null

    fastq_meta = [ meta, input_file ]

    return fastq_meta
}
