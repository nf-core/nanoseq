//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    input_path  // path: directory used to resolve relative 'input_file' entries, or [] if not needed

    main:
    SAMPLESHEET_CHECK ( samplesheet, input_path )
    def ch_sample = SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header: true, sep: ',' )
        .map { row -> getSampleInfo(row) }

    emit:
    sample = ch_sample // channel: [ val(meta), reads ] where meta = [ id, barcode, fast5 ]
}

//
// Convert a row of the validated samplesheet into [ meta, reads ]
//
def getSampleInfo(row) {
    def meta = [
        id     : row.sample,
        barcode: row.barcode,
        fast5  : row.fast5
    ]
    def reads = row.reads ? file(row.reads) : null
    return [ meta, reads ]
}
