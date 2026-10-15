/*
 * FastQ QC with NanoPlot, ToulligQC and FastQC
 */

include { NANOPLOT  } from '../../modules/nf-core/nanoplot/main'
include { TOULLIGQC } from '../../modules/nf-core/toulligqc/main'
include { FASTQC    } from '../../modules/nf-core/fastqc/main'

workflow QCFASTQ_NANOPLOT_FASTQC {
    take:
    ch_fastq       // channel: [ val(meta), path(fastq) ]
    skip_nanoplot  // boolean
    skip_toulligqc // boolean
    skip_fastqc    // boolean

    main:
    ch_multiqc_files = channel.empty()

    /*
     * FastQ QC using NanoPlot
     */
    ch_nanoplot_html = channel.empty()
    ch_nanoplot_txt  = channel.empty()
    if (!skip_nanoplot) {
        NANOPLOT(ch_fastq)
        ch_nanoplot_html = NANOPLOT.out.html
        ch_nanoplot_txt  = NANOPLOT.out.txt
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.map { _meta, txt -> txt })
    }

    /*
     * FastQ QC using ToulligQC
     */
    ch_toulligqc_report_html = channel.empty()
    ch_toulligqc_report_data = channel.empty()
    if (!skip_toulligqc) {
        TOULLIGQC(ch_fastq)
        ch_toulligqc_report_html = TOULLIGQC.out.report_html
        ch_toulligqc_report_data = TOULLIGQC.out.report_data
    }

    /*
     * FastQ QC using FastQC
     */
    ch_fastqc_zip  = channel.empty()
    ch_fastqc_html = channel.empty()
    if (!skip_fastqc) {
        FASTQC(ch_fastq)
        ch_fastqc_zip    = FASTQC.out.zip
        ch_fastqc_html   = FASTQC.out.html
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })
    }

    emit:
    nanoplot_html         = ch_nanoplot_html          // channel: [ val(meta), path(html) ]
    nanoplot_txt          = ch_nanoplot_txt           // channel: [ val(meta), path(txt) ]
    toulligqc_report_html = ch_toulligqc_report_html  // channel: [ val(meta), path(html) ]
    toulligqc_report_data = ch_toulligqc_report_data  // channel: [ val(meta), path(data) ]
    fastqc_zip            = ch_fastqc_zip             // channel: [ val(meta), path(zip) ]
    fastqc_html           = ch_fastqc_html            // channel: [ val(meta), path(html) ]
    multiqc_files         = ch_multiqc_files          // channel: path(file)
}
