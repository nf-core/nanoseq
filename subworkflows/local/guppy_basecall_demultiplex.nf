/*
 * Basecalling and demultiplexing fast5 using Guppy and ont_fast5_api
 */

params.guppy_options          = [:]
params.demux_fast5_options    = [:]

include { GUPPY         } from '../../modules/local/guppy'       	addParams( options: params.guppy_options    )
include { DEMUX_FAST5   } from '../../modules/local/demux_fast5'    addParams( options: params.demux_fast5_options    )

workflow GUPPY_BASECALL_DEMULTIPLEX {
    take:
	ch_input_path
	ch_sample_name
	ch_guppy_config
	ch_guppy_model
	ch_sample
	
    main:
    /*
     * Basecalling with guppy
     */
    GUPPY ( ch_input_path, ch_sample_name, ch_guppy_config.ifEmpty([]), ch_guppy_model.ifEmpty([]) )
    ch_guppy_summary = GUPPY.out.summary
    guppy_version = GUPPY.out.versions
	
	if (params.skip_demultiplexing){
		ch_sample
			.map { it -> [ it[0], it[0].id, it[2], it[3], it[4], it[5] ] }
			.set { ch_sample }
		}

	GUPPY.out.fastq
		.flatten()
		.map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.')) ] }
		.join(ch_sample, by: 1) // join on barcode
		.map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
		.set { ch_fastq }
	
	GUPPY.out.summary.view()
	demultiplexed_fast5 = GUPPY.out.fast5 
	input_summary = GUPPY.out.summary
    /*
     * Demultiplex fast5 files
     */	
    if (params.output_fast5 && !params.skip_demultiplexing_fast5) {
        DEMUX_FAST5 ( demultiplexed_fast5, input_summary )
    }            

	dumex_fast5_version = DEMUX_FAST5.out.versions 
	demultiplexed_fast5 = DEMUX_FAST5.out.fast5
	
    emit:
    guppy_version
    ch_guppy_summary
    ch_fastq
	//demultiplexed_fast5
    dumex_fast5_version
	
}
