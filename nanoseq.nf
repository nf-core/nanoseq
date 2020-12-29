////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

def ch_guppy_model  = Channel.empty()
def ch_guppy_config = Channel.empty()
include { GET_TEST_DATA         } from './modules/local/process/get_test_data'         addParams( options: [:] )
include { GET_NANOLYSE_FASTA    } from './modules/local/process/get_nanolyse_fasta'    addParams( options: [:] )
if (!params.skip_basecalling) { 

    // Pre-download test-dataset to get files for '--input_path' parameter
    // Nextflow is unable to recursively download directories via HTTPS
    if (workflow.profile.contains('test')) {
        if (!isOffline()) {
            GET_TEST_DATA ().set { ch_input_path }
        } else {
            exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download and run any test dataset!"
        }
    } else {
        if (params.input_path) { 
            ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true) 
        } 
    }

    // Need to stage guppy_config properly depending on whether its a file or string
    if (!params.guppy_config) {
        if (!params.flowcell) { exit 1, "Please specify a valid flowcell identifier for basecalling!" }
        if (!params.kit)      { exit 1, "Please specify a valid kit identifier for basecalling!"      }
    } else if (file(params.guppy_config).exists()) {
        ch_guppy_config = Channel.fromPath(params.guppy_config)
    }

    // Need to stage guppy_model properly depending on whether its a file or string
    if (params.guppy_model) {
        if (file(params.guppy_model).exists()) {
            ch_guppy_model = Channel.fromPath(params.guppy_model)
        }
    }

} else {
    if (!params.skip_demultiplexing) {
        if (!params.barcode_kit) {
            params.barcode_kit = 'Auto'
        }

        def qcatBarcodeKitList = ['Auto', 'RBK001', 'RBK004', 'NBD103/NBD104',
                                  'NBD114', 'NBD104/NBD114', 'PBC001', 'PBC096',
                                  'RPB004/RLB001', 'PBK004/LWB001', 'RAB204', 'VMK001', 'DUAL']

        if (params.barcode_kit && qcatBarcodeKitList.contains(params.barcode_kit)) {
            if (params.input_path) { 
                ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true) 
            } else { 
                exit 1, "Please specify a valid input fastq file to perform demultiplexing!" 
            }
        } else {
            exit 1, "Please provide a barcode kit to demultiplex with qcat. Valid options: ${qcatBarcodeKitList}"
        }
    }
}

if (params.run_nanolyse){
    if (!params.nanolyse_fasta){
        if (!isOffline()){
            GET_NANOLYSE_FASTA ().set { ch_nanolyse_fasta }
        } else {
            exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download lambda.fasta.gz file for running NanoLyse! Please explicitly specify --nanolyse_fasta."
        }
    } else {
        ch_nanolyse_fasta = file(params.nanolyse_fasta, checkIfExists: true)
    }
}

if (!params.skip_alignment) {
    if (params.aligner != 'minimap2' && params.aligner != 'graphmap2') {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'minimap2', 'graphmap2'"
    }
    if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
      exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
    }
}

if (!params.skip_quantification) {
    if (params.quantification_method != 'bambu' && params.quantification_method != 'stringtie2') {
        exit 1, "Invalid transcript quantification option: ${params.quantification_method}. Valid options: 'bambu', 'stringtie2'"
    }
    if (params.protocol != 'cDNA' && params.protocol != 'directRNA') {
        exit 1, "Invalid protocol option if performing quantification: ${params.protocol}. Valid options: 'cDNA', 'directRNA'"
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''
if (params.skip_alignment)  { multiqc_options['publish_dir'] = '' }

// TO DO -- define options for the processes below
def guppy_options    = modules['guppy']
def qcat_options     = modules['qcat']
def nanolyse_options = modules['nanolyse']
def bambu_options    = modules['bambu']

include { GUPPY                 } from './modules/local/process/guppy'                 addParams( options: guppy_options                )
include { QCAT                  } from './modules/local/process/qcat'                  addParams( options: qcat_options                 ) 
include { NANOLYSE              } from './modules/local/process/nanolyse'              addParams( options: nanolyse_options             )
include { BAM_RENAME            } from './modules/local/process/bam_rename'            addParams( options: [:]                          )
include { BAMBU                 } from './modules/local/process/bambu'                 addParams( options: bambu_options                )
include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions' addParams( options: [publish_files : ['csv':'']] )
include { MULTIQC               } from './modules/local/process/multiqc'               addParams( options: multiqc_options              )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
// TO DO -- define options for the subworkflows below
def pycoqc_options              = modules['pycoqc']
def nanoplot_summary_options    = modules['nanoplot_summary']
def nanoplot_fastq_options      = modules['nanoplot_fastq']
def fastqc_options              = modules['fastqc']
def genome_options              = modules['get_chrom_size']
def graphmap2_index_options     = modules['graphmap2_index']
def graphmap2_align_options     = modules['graphmap2_align']
def minimap2_index_options      = modules['minimap2_index']
def minimap2_align_options      = modules['minimap2_align']
def samtools_sort_options       = modules['samtools_sort']
def bigwig_options              = modules['ucsc_bed12tobigwig'] 
def bigbed_options              = modules['ucsc_bed12tobigbed']
def stringtie2_options          = modules['stringtie2']
def featurecounts_options       = modules['subread_featurecounts']
def deseq2_options              = modules['deseq2']
def dexseq_options              = modules['dexseq']

include { INPUT_CHECK                      } from './modules/local/subworkflow/input_check'                       addParams( options: [:] )
include { QCBASECALL_PYCOQC_NANOPLOT       } from './modules/local/subworkflow/qcbasecall_pycoqc_nanoplot'        addParams( pycoqc_options: pycoqc_options, nanoplot_summary_options: nanoplot_summary_options )
include { QCFASTQ_NANOPLOT_FASTQC          } from './modules/local/subworkflow/qcfastq_nanoplot_fastqc'           addParams( nanoplot_fastq_options: nanoplot_fastq_options, fastqc_options: fastqc_options )
include { PREPARE_GENOME                   } from './modules/local/subworkflow/prepare_genome'                    addParams( genome_options: genome_options )
include { ALIGN_GRAPHMAP2                  } from './modules/local/subworkflow/align_graphmap2'                   addParams( index_options: graphmap2_index_options, align_options: graphmap2_align_options, samtools_options: samtools_sort_options )
include { ALIGN_MINIMAP2                   } from './modules/local/subworkflow/align_minimap2'                    addParams( index_options: minimap2_index_options, align_options: minimap2_align_options, samtools_options: samtools_sort_options )
include { BEDTOOLS_UCSC_BIGWIG             } from './modules/local/subworkflow/bedtools_ucsc_bigwig'              addParams( bigwig_options: bigwig_options )
include { BEDTOOLS_UCSC_BIGBED             } from './modules/local/subworkflow/bedtools_ucsc_bigbed'              addParams( bigbed_options: bigbed_options )
include { QUANTIFY_STRINGTIE_FEATURECOUNTS } from './modules/local/subworkflow/quantify_stringtie_featurecounts'  addParams( stringtie2_options: stringtie2_options, featurecounts_options: featurecounts_options )
include { DIFFERENTIAL_DESEQ2_DEXSEQ       } from './modules/local/subworkflow/differential_deseq2_dexseq'        addParams( deseq2_options: deseq2_options, dexseq_options: dexseq_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
// TO DO -- define options for the processes below

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules (BAM_SORT_SAMTOOLS & BAM_STAT_SAMTOOLS are within two local subworkflows)
 */

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

workflow NANOSEQ{
    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    if (!params.skip_basecalling){
        ch_sample
            .first()
            .map { it[0] }
            .set { ch_sample_name }

        GUPPY ( ch_input_path, ch_sample_name, ch_guppy_config.ifEmpty([]), ch_guppy_model.ifEmpty([]) )
        ch_guppy_summary = GUPPY.out.summary

        if (params.skip_demultiplexing){
            ch_sample
                .map { it -> [ it[0], it[0], it[2], it[3], it[4], it[5] ] }
                .set { ch_sample }
        }

        GUPPY.out.fastq
           .flatten()
           .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.')) ] }
           .join(ch_sample, by: 1) // join on barcode
           .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
           .set { ch_fastq }
    } else {
        ch_guppy_summary = Channel.empty()
        
        if (!params.skip_demultiplexing){
            QCAT ( ch_input_path )
            ch_fastq = Channel.empty()
            QCAT.out.fastq
               .flatten()
               .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.'))] }
               .join(ch_sample, by: 1) // join on barcode
               .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
               .set { ch_fastq } 
        } else {
            if (!params.skip_alignment){
               ch_sample
                   .map { it -> if (it[6].toString().endsWith('.gz')) [ it[0], it[6], it[2], it[1], it[4], it[5] ] }
                   .set { ch_fastq }
            } else {
               ch_fastq = Channel.empty()
            }
        }
    }

    if (params.run_nanolyse) {
       ch_fastq
           .map { it -> [ it[0], it[1] ] }
           .set { ch_fastq_nanolyse }
       NANOLYSE ( ch_fastq_nanolyse, ch_nanolyse_fasta )
       NANOLYSE.out.nanolyse_fastq
          .join( ch_sample )
          .map { it -> [ it[0], it[1], it[3], it[4], it[5], it[6] ]}
          .set { ch_fastq }
    }

    if (!params.skip_qc){
        if (!params.skip_basecalling){
            QCBASECALL_PYCOQC_NANOPLOT ( ch_guppy_summary , params.skip_pycoqc, params.skip_nanoplot )
        }
        QCFASTQ_NANOPLOT_FASTQC ( ch_fastq, params.skip_nanoplot, params.skip_fastqc)        
    }

    if (!params.skip_alignment){
       PREPARE_GENOME ( ch_fastq )
       ch_fasta_index = PREPARE_GENOME.out.ch_fasta_index
       ch_gtf_bed     = PREPARE_GENOME.out.ch_gtf_bed
       if (params.aligner == 'minimap2') {
          ALIGN_MINIMAP2 ( ch_fasta_index, ch_fastq )
          ch_view_sortbam = ALIGN_MINIMAP2.out.ch_sortbam
       } else {
          ALIGN_GRAPHMAP2 ( ch_fasta_index, ch_fastq )
          ch_view_sortbam = ALIGN_GRAPHMAP2.out.ch_sortbam
       }
       if (!params.skip_bigwig){
          BEDTOOLS_UCSC_BIGWIG ( ch_view_sortbam )
       }
       if (!params.skip_bigbed){
          BEDTOOLS_UCSC_BIGBED ( ch_view_sortbam )
       }
       ch_view_sortbam
          .map { it -> [ it[0], it[3] ] }
          .set { ch_sortbam }
    } else {
       ch_sample
           .map { it -> if (it[6].toString().endsWith('.bam')) [ it[0], it[6] ] }
           .set { ch_sample_bam }
       BAM_RENAME ( ch_sample_bam )
       ch_sortbam = BAM_RENAME.out.sortbam_quant
    }

    if (!params.skip_quantification && (params.protocol == 'cDNA' || params.protocol == 'directRNA')) {
       if (params.quantification_method == 'bambu') {
          ch_sample
              .map { it -> [ it[2], it[3] ]}
              .unique()
              .set { ch_sample_annotation }          
          BAMBU ( ch_sample_annotation, ch_sortbam.collect{ it [1] } )
          ch_gene_counts       = BAMBU.out.ch_gene_counts
          ch_transcript_counts = BAMBU.out.ch_transcript_counts
       } else {
          QUANTIFY_STRINGTIE_FEATURECOUNTS( ch_sample, ch_sortbam )
          ch_gene_counts       = QUANTIFY_STRINGTIE_FEATURECOUNTS.out.ch_gene_counts
          ch_transcript_counts = QUANTIFY_STRINGTIE_FEATURECOUNTS.out.ch_transcript_counts
       }
       if (!params.skip_differential_analysis) {
          DIFFERENTIAL_DESEQ2_DEXSEQ( ch_gene_counts, ch_transcript_counts )
       }
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
//    Completion.email(workflow, params, params.summary_params, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
