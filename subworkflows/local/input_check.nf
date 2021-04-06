/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { CHECK_SAMPLESHEET;
          get_sample_info } from '../../modules/local/check_samplesheet' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    CHECK_SAMPLESHEET ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it, params.genomes) }
        .map { it -> [ it[0], it[2], it[3], it[4], it[5], it[6], it[1] , it[7] ] }
        .set { ch_sample }

    emit:
    ch_sample // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
}
