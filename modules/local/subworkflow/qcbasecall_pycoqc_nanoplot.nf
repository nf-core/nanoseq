/*
 * Basecalling QC with PycoQC and NanoPlot
 */

params.pycoqc_options             = [:]
params.nanoplot_summary_options   = [:]

include { PYCOQC              } from '../process/pycoqc'             addParams( options: params.pycoqc_options )
include { NANOPLOT_SUMMARY    } from '../process/nanoplot_summary'   addParams( options: params.nanoplot_summary_options )

workflow QCBASECALL_PYCOQC_NANOPLOT {
    take:
    ch_guppy_summary_txt // channel: [ val(meta), [ reads ] ]
    skip_pycoqc
    skip_nanoplot
    
    main:
    /*
     * QC using PycoQC
     */
    pycoqc_html    = Channel.empty()
    pycoqc_json    = Channel.empty()
    pycoqc_version = Channel.empty()
    if (!skip_pycoqc){
       PYCOQC ( ch_guppy_summary_txt )
       pycoqc_html    = PYCOQC.out.html
       pycoqc_json    = PYCOQC.out.pycoqc_multiqc
       pycoqc_version = PYCOQC.out.version
    }

    /*
     * QC using NanoPlot
     */
    nanoplot_png    = Channel.empty()
    nanoplot_html   = Channel.empty()
    nanoplot_txt    = Channel.empty()
    nanoplot_log    = Channel.empty()
    if (!skip_nanoplot){
       NANOPLOT_SUMMARY ( ch_guppy_summary_txt )
       nanoplot_png     = NANOPLOT_SUMMARY.out.png
       nanoplot_html    = NANOPLOT_SUMMARY.out.html
       nanoplot_txt     = NANOPLOT_SUMMARY.out.txt
       nanoplot_log     = NANOPLOT_SUMMARY.out.log
    }

    emit:
    pycoqc_html
    pycoqc_json
    pycoqc_version

    nanoplot_png
    nanoplot_html
    nanoplot_txt
    nanoplot_log
}
