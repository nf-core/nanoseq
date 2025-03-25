include { DORADO_DEMUX } from '../modules/local/dorado/demux'
include { DORADO_ALIGN } from '../modules/local/dorado/align'

workflow DORADO {
    take:
    path samplesheet
    path input_files //path to dir containing fastqs/ubams
    path reference_fasta
    // path reference_fai
    val barcode_kit
    val barcode_both_ends

    main:
    //TODO make demux optional
    //Demultiplexing
    DORADO_DEMUX(samplesheet, fastq, barcode_kit, barcode_both_ends)
    sample_fastq = DORADO_DEMUX.out.fastq //list of all sample fastqs

    /*
    * convert sample_fastq into [meta, fastq]
    */

    //create meta maps per sample from the samplesheet
    meta = samplesheet.splitCsv(header: true)
                .map { row -> 
                    def meta = [
                        id: row.sample_id, 
                        position: row.position_id, 
                        flow_cell: row.flow_cell_id, 
                        kit: row.kit, 
                        experiment: row.experiment_id, 
                        barcode: row.barcode, 
                        alias: row.alias
                    ]
                    return [row.alias, meta]  // Return metadata and alias for matching
                }
    //restructure fastq channel into alias, fastq 
    sample_fastq_alias = sample_fastq.map { file -> 
        def alias = file.name.replaceAll('.fastq.gz','')
        return [alias, file]
    }
    //join sample fastqs and meta maps using alias as matching key
    sample_fastq_ch = meta.join(sample_fastq_alias)

    //chopper fastq trimming
    CHOPPER(sample_fastq_ch)
    sample_clean_fastq = CHOPPER.out.clean_fastq

    DORADO_ALIGN(sample_fastq, reference_fasta)
    sample_bams = DORADO_ALIGN.out.bam //{[meta1, bam1], [meta2, bam2]}

    emit:
    sample_bams
}
