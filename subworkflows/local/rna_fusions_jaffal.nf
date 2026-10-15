/*
 * RNA FUSION DETECTION WITH JAFFAL
 */

include { WGET as WGET_JAFFAL_REF } from '../../modules/nf-core/wget/main'
include { UNTAR          } from '../../modules/nf-core/untar/main'
include { JAFFAL         } from '../../modules/local/jaffal'

workflow RNA_FUSIONS_JAFFAL {
    take:
    ch_fastq       // channel: [ val(meta), path(fastq) ]
    jaffal_ref_dir // val: path to a JAFFAL reference directory, or null to download it

    main:
    if (jaffal_ref_dir) {
        ch_jaffal_ref_dir = channel.value(file(jaffal_ref_dir, checkIfExists: true))
    } else {
        /*
         * Download the JAFFAL reference
         */
        WGET_JAFFAL_REF( channel.of([ [ id: 'jaffal_ref' ], 'https://ndownloader.figshare.com/files/28168755', 'tar.gz' ]) )

        /*
         * Untar the JAFFAL reference
         */
        UNTAR( WGET_JAFFAL_REF.out.outfile )
        ch_jaffal_ref_dir = UNTAR.out.untar
            .map { _meta, ref_dir -> ref_dir }
            .first()
    }

    /*
     * Detect fusion transcripts with JAFFAL
     */
    JAFFAL( ch_fastq, ch_jaffal_ref_dir )

    emit:
    fusions = JAFFAL.out.jaffal_results // channel: path(csv)
    fasta   = JAFFAL.out.jaffal_fastq   // channel: [ val(meta), path(fasta) ]
}
