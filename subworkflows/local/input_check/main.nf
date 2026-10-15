//
// Read the input samplesheet and resolve input files for each sample
//
// Each row is validated against assets/schema_input.json with nf-schema. Checks that span
// rows, and the resolution of raw data directories, are done here.
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet // channel (value): path to samplesheet
    input_path  // channel (value): directory used to resolve 'input_file' entries by name, or [] if not needed

    main:
    def ch_base_dir = input_path instanceof List ? channel.value('') : input_path

    def ch_sample = samplesheet
        .combine(ch_base_dir)
        .flatMap { sheet, base_dir ->
            def rows = samplesheetToList(sheet, "${projectDir}/assets/schema_input.json")
                .collect { meta, input_file -> getSampleInfo(meta, input_file, base_dir) }
            validateSamples(rows)
            return rows
        }

    emit:
    sample = ch_sample // channel: [ val(meta), reads ] where meta = [ id, group, replicate, barcode, fast5 ]
}

//
// Convert a samplesheet row into [ meta, reads ]
//
def getSampleInfo(row_meta, input_file, base_dir) {
    def meta = [
        id       : "${row_meta.group}_R${row_meta.replicate}".toString(),
        group    : row_meta.group,
        replicate: row_meta.replicate,
        barcode  : row_meta.barcode ? "barcode${row_meta.barcode.toString().padLeft(2, '0')}".toString() : '',
        fast5    : ''
    ]
    if (!input_file) {
        return [ meta, null ]
    }

    // Entries are resolved by name against a staged directory (e.g. cloned test data)
    def input = base_dir ? "${base_dir}/${input_file.toString().tokenize('/')[-1]}".toString() : input_file.toString()

    if (input ==~ /.*\.(f(ast)?q\.gz|bam|pod5)$/) {
        return [ meta, file(input, checkIfExists: true) ]
    }
    if (input.endsWith('.fast5')) {
        return [ meta + [ fast5: input ], file(input, checkIfExists: true) ]
    }

    // Directory of raw data: either only FAST5/POD5 files, or 'fast5/' and 'fastq/' subdirectories
    def dir = file(input, type: 'dir', checkIfExists: true)
    def names = files("${dir}/*").collect { f -> f.name }
    if (names && (names.every { n -> n.endsWith('.fast5') } || names.every { n -> n.endsWith('.pod5') })) {
        return [ meta + [ fast5: dir.toString() ], dir ]
    }
    if (names.containsAll([ 'fast5', 'fastq' ])) {
        def fast5_files = files("${dir}/fast5/*")
        if (!fast5_files || !fast5_files.every { f -> f.name.endsWith('.fast5') }) {
            error("Please check input samplesheet -> '${dir}/fast5' must only contain FAST5 files: ${meta.id}")
        }
        def fastq = files("${dir}/fastq/*.{fastq,fq}.gz")
        if (fastq.size() != 1) {
            error("Please check input samplesheet -> '${dir}/fastq' must contain exactly one '.fastq.gz' or '.fq.gz' file: ${meta.id}")
        }
        return [ meta + [ fast5: "${dir}/fast5".toString() ], fastq[0] ]
    }
    error("Please check input samplesheet -> '${input_file}' is not a '.fastq.gz', '.fq.gz', '.bam', '.pod5' or '.fast5' file, or a directory of FAST5/POD5 files or with 'fast5' and 'fastq' subdirectories: ${meta.id}")
}

//
// Checks that span samplesheet rows
//
def validateSamples(samples) {
    // Replicates of each group must be numbered 1..<number of replicates>
    samples.groupBy { meta, _reads -> meta.group }.each { group, group_samples ->
        def replicates = group_samples.collect { meta, _reads -> meta.replicate as Integer }.sort()
        if (replicates != (1..replicates.size()).toList()) {
            error("Please check input samplesheet -> Replicate ids must start with 1..<number of replicates>: group '${group}'")
        }
    }

    // All input files must be of the same type
    def extensions = samples
        .findAll { _meta, reads -> reads }
        .collect { _meta, reads -> reads.name.replaceFirst(/.*?\.((f(ast)?q\.gz)|bam|pod5|fast5)$/, '$1').replace('fq.gz', 'fastq.gz') }
        .collect { ext -> ext in [ 'fastq.gz', 'bam', 'pod5', 'fast5' ] ? ext : 'directory' }
        .unique()
    if (extensions.size() > 1) {
        error("Please check input samplesheet -> All input files must be of the same type. Found: ${extensions.join(', ')}")
    }
}
