<<<<<<< HEAD
//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )
=======
/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { CHECK_SAMPLESHEET } from '../../modules/local/check_samplesheet' addParams( options: params.options )
>>>>>>> origin/dev

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
<<<<<<< HEAD
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
=======
    CHECK_SAMPLESHEET ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it, params.genomes) }
        .map { it -> [ it[0], it[2], it[3], it[4], it[5], it[6], it[1] , it[7] ] }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str ]
def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {
    def meta = [:]
    meta.id  = sample.sample

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

    return [ meta, input_file, sample.barcode, fasta, gtf, sample.is_transcripts.toBoolean(), fasta.toString()+';'+gtf.toString() ]
>>>>>>> origin/dev
}
