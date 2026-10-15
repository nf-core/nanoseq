/*
 * RNA FUSION DETECTION WITH JAFFAL
 */

include { GET_JAFFAL_REF } from '../../modules/local/get_jaffal_ref'
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
        GET_JAFFAL_REF()

        /*
         * Untar the JAFFAL reference
         */
        UNTAR( GET_JAFFAL_REF.out.ref )
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
