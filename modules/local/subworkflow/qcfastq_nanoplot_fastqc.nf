/*
 * FastQ QC with NanoPlot and fastqc
 */

params.nanoplot_fastq_options = [:]
params.fastqc_options         = [:]

include { NANOPLOT_FASTQ     } from '../process/nanoplot_fastq'            addParams( options: params.nanoplot_fastq_options )
include { FASTQC             } from '../process/fastqc'                    addParams( options: params.fastqc_options )

workflow QCFASTQ_NANOPLOT_FASTQC {
    take:
    ch_fastq
    skip_nanoplot
    skip_fastqc    

    main:
    ch_fastq.map{ ch -> [ ch[0], ch[1] ] }

    /*
     * FastQ QC using NanoPlot
     */
    nanoplot_png     = Channel.empty()
    nanoplot_html    = Channel.empty()
    nanoplot_txt     = Channel.empty()
    nanoplot_log     = Channel.empty()
    nanoplot_version = Channel.empty()
    if (!skip_nanoplot){
       NANOPLOT_FASTQ ( ch_fastq )
       nanoplot_png     = NANOPLOT_FASTQ.out.png
       nanoplot_html    = NANOPLOT_FASTQ.out.html
       nanoplot_txt     = NANOPLOT_FASTQ.out.txt
       nanoplot_log     = NANOPLOT_FASTQ.out.log
       nanoplot_version = NANOPLOT_FASTQ.out.version
    }

    /*
     * FastQ QC using FASTQC
     */
    fastqc_zip     = Channel.empty()
    fastqc_html    = Channel.empty()
    fastqc_multiqc = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc){
       FASTQC ( ch_fastq )
       fastqc_zip     = FASTQC.out.zip
       fastqc_html    = FASTQC.out.html  
       fastqc_multiqc = fastqc_multiqc.mix( fastqc_zip, fastqc_html )
       fastqc_version = FASTQC.out.version
    }
    
    emit:
    nanoplot_png
    nanoplot_html
    nanoplot_txt
    nanoplot_log
    nanoplot_version

    fastqc_zip
    fastqc_html
    fastqc_version
    fastqc_multiqc
}
