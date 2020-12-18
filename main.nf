#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/nanoseq
========================================================================================
 nf-core/nanoseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/nanoseq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow run nf-core/nanoseq \
        --input samplesheet.csv \
        --protocol DNA \
        --input_path ./fast5/ \
        --flowcell FLO-MIN106 \
        --kit SQK-LSK109 \
        --barcode_kit SQK-PBK004 \
        -profile docker

    Mandatory arguments
      --input [file]                      Comma-separated file containing information about the samples in the experiment (see docs/usage.md)
      --protocol [str]                    Specifies the type of sequencing i.e. "DNA", "cDNA" or "directRNA"
      -profile [str]                      Configuration profile to use. Can use multiple (comma separated)
                                          Available: docker, singularity, test, awsbatch, <institute> and more.

    Basecalling/Demultiplexing
      --input_path [file]                 Path to Nanopore run directory (e.g. fastq_pass/) or a basecalled fastq file that requires demultiplexing
      --flowcell [str]                    Flowcell used to perform the sequencing e.g. FLO-MIN106. Not required if '--guppy_config' is specified
      --kit [str]                         Kit used to perform the sequencing e.g. SQK-LSK109. Not required if '--guppy_config' is specified
      --barcode_kit [str]                 Barcode kit used to perform the sequencing e.g. SQK-PBK004
      --barcode_both_ends [bool]          Require barcode on both ends for Guppy basecaller (Default: false)
      --guppy_config [file/str]           Guppy config file used for basecalling. Cannot be used in conjunction with '--flowcell' and '--kit'
      --guppy_model [file/str]            Custom basecalling model file (JSON) to use for Guppy basecalling, such as the output from Taiyaki (Default: false)
      --guppy_gpu [bool]                  Whether to perform basecalling with Guppy in GPU mode (Default: false)
      --guppy_gpu_runners [int]           Number of '--gpu_runners_per_device' used for Guppy when using '--guppy_gpu' (Default: 6)
      --guppy_cpu_threads [int]           Number of '--cpu_threads_per_caller' used for Guppy when using '--guppy_gpu' (Default: 1)
      --gpu_device [str]                  Basecalling device specified to Guppy in GPU mode using '--device' (Default: 'auto')
      --gpu_cluster_options [str]         Cluster options required to use GPU resources (e.g. '--part=gpu --gres=gpu:1')
      --qcat_min_score [int]              Minimum scores of '--min-score' used for qcat (Default: 60)
      --qcat_detect_middle [bool]         Search for adapters in the whole read '--detect-middle' used for qcat (Default: false)
      --skip_basecalling [bool]           Skip basecalling with Guppy (Default: false)
      --skip_demultiplexing [bool]        Skip demultiplexing with Guppy (Default: false)
    
    Raw read cleaning
      --run_nanolyse [bool]               Filter reads from FastQ files mapping to the lambda phage genome using NanoLyse (Default: false)
      --nanolyse_fasta [file]             Provide a fasta file for NanoLyse to use for filtering FastQ files (Default: false)

    Alignment
      --aligner [str]                     Specifies the aligner to use (available are: minimap2 or graphmap2). Both are capable of performing unspliced/spliced alignment (Default: 'minimap2')
      --stranded [bool]                   Specifies if the data is strand-specific. Automatically activated when using '--protocol directRNA' (Default: false)
      --save_align_intermeds [bool]       Save the '.sam' files from the alignment step (Default: false)
      --skip_alignment [bool]             Skip alignment and subsequent process (Default: false)

    Coverage tracks
      --skip_bigwig [bool]                Skip BigWig file generation (Default: false)
      --skip_bigbed [bool]                Skip BigBed file generation (Default: false)

    Quantification and differential analysis
      --quantification_method [str]       Specifies the transcript quantification method to use (available are: bambu or stringtie2). Only available when protocol is cDNA or directRNA.
      --skip_quantification [bool]        Skip transcript quantification and differential analysis (Default: false)
      --skip_differential_analysis [bool] Skip differential analysis with DESeq2 and DEXSeq (Default: false)

    QC
      --skip_qc [bool]                    Skip all QC steps apart from MultiQC (Default: false)
      --skip_pycoqc [bool]                Skip pycoQC (Default: false)
      --skip_nanoplot [bool]              Skip NanoPlot (Default: false)
      --skip_fastqc [bool]                Skip FastQC (Default: false)
      --skip_multiqc [bool]               Skip MultiQC (Default: false)
           
    Other
      --outdir [file]                     The output directory where the results will be saved (Default: '/results')
      --publish_dir_mode [str]            Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]             Same as --email, except only send mail if the workflow is not successful (Default: false)
      --max_multiqc_email_size [str]      Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch
      --awsqueue [str]                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                   The AWS Region for your AWS Batch job to run on (Default: 'eu-west-1')
      --awscli [str]                      Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

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
if (!params.skip_basecalling) {

    // Pre-download test-dataset to get files for '--input_path' parameter
    // Nextflow is unable to recursively download directories via HTTPS
    if (workflow.profile.contains('test')) {
        if (!isOffline()) {
            process GET_TEST_DATA {
                output:
                path "test-datasets/fast5/$barcoded/" into ch_input_path

                script:
                barcoded = workflow.profile.contains('test_nonbc') ? "nonbarcoded" : "barcoded"
                """
                git clone https://github.com/nf-core/test-datasets.git --branch nanoseq --single-branch
                """
            }
        } else {
            exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download and run any test dataset!"
        }
    } else {
        if (params.input_path) { ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true) } else { exit 1, "Please specify a valid run directory to perform basecalling!" }
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
            if (params.input_path) { ch_input_path = Channel.fromPath(params.input_path, checkIfExists: true) } else { exit 1, "Please specify a valid input fastq file to perform demultiplexing!" }
        } else {
            exit 1, "Please provide a barcode kit to demultiplex with qcat. Valid options: ${qcatBarcodeKitList}"
        }
    }
}

if (params.run_nanolyse){
    if (!params.nanolyse_fasta){
        if (!isOffline()){
            process GET_NANOLYSE_FASTA {
                output:
                path "*fasta.gz" into ch_nanolyse_fasta

                script:
                """
                wget https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz
                """
            }
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

// Stage config files
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs           = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images    = file("$projectDir/docs/images/", checkIfExists: true)

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']               = custom_runName ?: workflow.runName
summary['Samplesheet']            = params.input
summary['Protocol']               = params.protocol
summary['Stranded']               = (params.stranded || params.protocol == 'directRNA') ? 'Yes' : 'No'
summary['Skip Basecalling']       = params.skip_basecalling ? 'Yes' : 'No'
summary['Skip Demultiplexing']    = params.skip_demultiplexing ? 'Yes' : 'No'
if (!params.skip_basecalling) {
    summary['Run Dir']            = params.input_path
    summary['Flowcell ID']        = params.flowcell ?: 'Not required'
    summary['Kit ID']             = params.kit ?: 'Not required'
    summary['Barcode Kit ID']     = params.barcode_kit ?: 'Unspecified'
    summary['Barcode Both Ends']  = params.barcode_both_ends ? 'Yes' : 'No'
    summary['Guppy Config File']  = params.guppy_config ?: 'Unspecified'
    summary['Guppy Model File']   = params.guppy_model ?:'Unspecified'
    summary['Guppy GPU Mode']     = params.guppy_gpu ? 'Yes' : 'No'
    summary['Guppy GPU Runners']  = params.guppy_gpu_runners
    summary['Guppy CPU Threads']  = params.guppy_cpu_threads
    summary['Guppy GPU Device']   = params.gpu_device ?: 'Unspecified'
    summary['Guppy GPU Options']  = params.gpu_cluster_options ?: 'Unspecified'
}
if (params.skip_basecalling && !params.skip_demultiplexing) {
    summary['Input FastQ File']   = params.input_path
    summary['Barcode Kit ID']     = params.barcode_kit ?: 'Unspecified'
    summary['Qcat Min Score']     = params.qcat_min_score
    summary['Qcat Detect Middle'] = params.qcat_detect_middle ? 'Yes': 'No'
}
if (params.run_nanolyse) {
    summary['Run NanoLyse']       = params.run_nanolyse
    summary['NanoLyse Fasta']     = params.nanolyse_fasta
}
if (!params.skip_alignment) {
    summary['Aligner']            = params.aligner
    summary['Save Intermeds']     = params.save_align_intermeds ? 'Yes' : 'No'
}
if (!params.skip_quantification && params.protocol != 'DNA') {
    summary['Quantification Method'] = params.quantification_method
    summary['Skip Diff Analysis']    = params.skip_differential_analysis
}
if (params.skip_alignment) summary['Skip Alignment'] = 'Yes'
if (params.skip_bigbed)    summary['Skip BigBed']    = 'Yes'
if (params.skip_bigwig)    summary['Skip BigWig']    = 'Yes'
if (params.skip_qc)        summary['Skip QC']        = 'Yes'
if (params.skip_pycoqc)    summary['Skip PycoQC']    = 'Yes'
if (params.skip_nanoplot)  summary['Skip NanoPlot']  = 'Yes'
if (params.skip_fastqc)    summary['Skip FastQC']    = 'Yes'
if (params.skip_multiqc)   summary['Skip MultiQC']   = 'Yes'
summary['Max Resources']          = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']             = params.outdir
summary['Launch dir']             = workflow.launchDir
summary['Working dir']            = workflow.workDir
summary['Script dir']             = workflow.projectDir
summary['User']                   = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']         = params.awsregion
    summary['AWS Queue']          = params.awsqueue
    summary['AWS CLI']            = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']     = params.email
    summary['E-mail on failure']  = params.email_on_fail
    summary['MultiQC maxsize']    = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

/*
 * Check input samplesheet
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet from ch_input
    
    output:
    path "*reformat.csv" into ch_samplesheet_reformat
    
    script:  // This script is bundled with the pipeline, in nf-core/nanoseq/bin/
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}

// Function to resolve fasta and gtf file if using iGenomes
// Returns [ sample, input_file, barcode, fasta, gtf, is_transcripts, annotation_str ]
def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {

    // Resolve fasta and gtf file if using iGenomes
    def fasta = false
    def gtf   = false
    if (sample.genome) {
        if (genomeMap && genomeMap.containsKey(sample.genome)) {
            fasta = file(genomeMap[sample.genome].fasta, checkIfExists: true)
            gtf   = file(genomeMap[sample.genome].gtf, checkIfExists: true)
        } else {
            fasta = file(sample.genome, checkIfExists: true)
        }
    }

    // Check if input file and gtf file exists
    input_file = sample.input_file ? file(sample.input_file, checkIfExists: true) : null
    gtf        = sample.gtf        ? file(sample.gtf, checkIfExists: true)        : gtf

    return [ sample.sample, input_file, sample.barcode, fasta, gtf, sample.is_transcripts.toBoolean(), fasta.toString()+';'+gtf.toString() ]
}

// Create channels = [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it, params.genomes) }
    .map { it -> [ it[0], it[2], it[3], it[4], it[5], it[6] ] }
    .into {
        ch_sample_fasta
        ch_sample_gtf
        ch_sample_replicates
        ch_sample_groups
        ch_sample_info
        ch_sample_name
        ch_sample_annotation
        ch_sample_nanolyse
    }

// Check if reference genome and annotation are the same for all samples if perfoming quantification
// Check if we have replicates and multiple conditions in the input samplesheet
def REPLICATES_EXIST    = false
def MULTIPLE_CONDITIONS = false
if (!params.skip_quantification) {

    def fastas = ch_sample_fasta.map { it[2] }.unique().toList().val
    def gtfs   = ch_sample_gtf.map   { it[3] }.unique().toList().val
    if (fastas.size() != 1 || gtfs.size() != 1 || fastas[0] == false || gtfs[0] == false) {
        exit 1, """Quantification can only be performed if all samples in the samplesheet have the same reference fasta and GTF file."
                   Please specify the '--skip_quantification' parameter if you wish to skip these steps."""
    }

    REPLICATES_EXIST    = ch_sample_replicates.map { it -> it[0].split('_')[-1].replaceAll('R','').toInteger() }.max().val > 1
    MULTIPLE_CONDITIONS = ch_sample_groups.map { it -> it[0].split('_')[0..-2].join('_') }.unique().count().val > 1
}

if (!params.skip_basecalling) {

    // Get sample name for single sample when --skip_demultiplexing
    ch_sample_name
        .first()
        .map { it[0] }
        .set { ch_sample_name }

    /*
     * Basecalling and demultipexing using Guppy
     */
    process GUPPY {
        tag "$input_path"
        label 'process_high'
        publishDir path: "${params.outdir}/guppy", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (!filename.endsWith("guppy.txt")) filename
                    }

        input:
        path input_path from ch_input_path
        val name from ch_sample_name
        path guppy_config from ch_guppy_config.ifEmpty([])
        path guppy_model from ch_guppy_model.ifEmpty([])

        output:
        path "fastq/*.fastq.gz" into ch_fastq
        path "basecalling/*.txt" into ch_guppy_pycoqc_summary,
                                      ch_guppy_nanoplot_summary
        path "basecalling/*"
        path "v_guppy.txt" into ch_guppy_version

        script:
        barcode_kit  = params.barcode_kit ? "--barcode_kits $params.barcode_kit" : ""
        barcode_ends = params.barcode_both_ends ? "--require_barcodes_both_ends" : ""
        proc_options = params.guppy_gpu ? "--device $params.gpu_device --num_callers $task.cpus --cpu_threads_per_caller $params.guppy_cpu_threads --gpu_runners_per_device $params.guppy_gpu_runners" : "--num_callers 2 --cpu_threads_per_caller ${task.cpus/2}"
        def config   = "--flowcell $params.flowcell --kit $params.kit"
        if (params.guppy_config) config = file(params.guppy_config).exists() ? "--config ./$guppy_config" : "--config $params.guppy_config"
        def model    = ""
        if (params.guppy_model)  model  = file(params.guppy_model).exists() ? "--model ./$guppy_model" : "--model $params.guppy_model"
        """
        guppy_basecaller \\
            --input_path $input_path \\
            --save_path ./basecalling \\
            --records_per_fastq 0 \\
            --compress_fastq \\
            $barcode_kit \\
            $proc_options \\
            $barcode_ends \\
            $config \\
            $model
        guppy_basecaller --version &> v_guppy.txt

        ## Concatenate fastq files
        mkdir fastq
        cd basecalling
        if [ "\$(find . -type d -name "barcode*" )" != "" ]
        then
            for dir in barcode*/
            do
                dir=\${dir%*/}
                cat \$dir/*.fastq.gz > ../fastq/\$dir.fastq.gz
            done
        else
            cat *.fastq.gz > ../fastq/${name}.fastq.gz
        fi
        """
    }

    // Create channels = [ sample, fastq, fasta, gtf, is_transcripts, annotation_str ]
    ch_fastq
        .flatten()
        .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.')) ] } // [barcode001.fastq, barcode001]
        .join(ch_sample_info, by: 1) // join on barcode
        .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
        .into { 
            ch_fastq_nanoplot
            ch_fastq_fastqc
            ch_fastq_sizes
            ch_fastq_gtf
            ch_fastq_index
            ch_fastq_align 
            ch_fastq_nanolyse
        }

} else {
    if (!params.skip_demultiplexing) {

        /*
         * Demultipexing using qcat
         */
        process QCAT {
            tag "$input_path"
            label 'process_medium'
            publishDir path: "${params.outdir}/qcat", mode: params.publish_dir_mode

            input:
            path input_path from ch_input_path

            output:
            path "fastq/*.fastq.gz" into ch_fastq

            script:
            detect_middle = params.qcat_detect_middle ? "--detect-middle $params.qcat_detect_middle" : ""
            """
            ## Unzip fastq file
            ## qcat doesnt support zipped files yet
            FILE=$input_path
            if [[ \$FILE == *.gz ]]
            then
                zcat $input_path > unzipped.fastq
                FILE=unzipped.fastq
            fi

            qcat  \\
                -f \$FILE \\
                -b ./fastq \\
                --kit $params.barcode_kit \\
                --min-score $params.qcat_min_score \\
                $detect_middle
            
            ## Zip fastq files
            pigz -p $task.cpus fastq/*
            """
        }

        // Create channels = [ sample, fastq, fasta, gtf, is_transcripts, annotation_str ]
        ch_fastq
            .flatten()
            .map { it -> [ it, it.baseName.substring(0,it.baseName.lastIndexOf('.'))] } // [barcode001.fastq, barcode001]
            .join(ch_sample_info, by: 1) // join on barcode
            .map { it -> [ it[2], it[1], it[3], it[4], it[5], it[6] ] }
            .into { 
                ch_fastq_nanoplot
                ch_fastq_fastqc
                ch_fastq_sizes
                ch_fastq_gtf
                ch_fastq_index
                ch_fastq_align
                ch_fastq_nanolyse 
            }
    } else {
        if (!params.skip_alignment) {
            // Create channels = [ sample, fastq, fasta, gtf, is_transcripts, annotation_str ]
            ch_samplesheet_reformat
                .splitCsv(header:true, sep:',')
                .map { get_sample_info(it, params.genomes) }
                .map { it -> if (it[1].toString().endsWith('.gz')) [ it[0], it[1], it[3], it[4], it[5], it[6] ] }
                .into { 
                    ch_fastq_nanoplot
                    ch_fastq_fastqc
                    ch_fastq_sizes
                    ch_fastq_gtf
                    ch_fastq_index
                    ch_fastq_align 
                    ch_fastq_nanolyse
                }
        } else {
            ch_fastq_nanoplot = Channel.empty()
            ch_fastq_fastqc   = Channel.empty()
        }
    }
    ch_guppy_version          = Channel.empty()
    ch_guppy_pycoqc_summary   = Channel.empty()
    ch_guppy_nanoplot_summary = Channel.empty()
}


if (params.run_nanolyse){
    /*
     * DNA contaminant removal using NanoLyse
     */
    process NANOLYSE {
        tag "$sample"
        label 'process_medium'
        publishDir "${params.outdir}/nanolyse", mode: params.publish_dir_mode,
            saveAs: { filename ->
                            if (filename.endsWith(".log")) "log/$filename"
                            else filename
                    }

        input:
        tuple val(sample), path(fastq) from ch_fastq_nanolyse.map{ ch -> [ ch[0], ch[1] ] }
        path fasta from ch_nanolyse_fasta
    
        output:
        tuple val(sample), path("*.fastq.gz") into ch_nanolyse_fastq
        path "*.log"

        script:
        """
        gunzip -c $fastq | NanoLyse -r $fasta | gzip > ${sample}.clean.fastq.gz
        cp NanoLyse.log ${sample}.nanolyse.log
        """
    }

    ch_nanolyse_fastq
        .join(ch_sample_nanolyse)
        .map { it -> [ it[0], it[1], it[3], it[4], it[5], it[6] ] }
        .into {
            ch_fastq_nanoplot
            ch_fastq_fastqc
            ch_fastq_sizes
            ch_fastq_gtf
            ch_fastq_index
            ch_fastq_align
        }
}

/*
 * QC using PycoQC
 */
process PYCOQC {
    tag "$summary_txt"
    label 'process_low'
    publishDir "${params.outdir}/pycoqc", mode: params.publish_dir_mode

    when:
    !params.skip_basecalling && !params.skip_qc && !params.skip_pycoqc

    input:
    path summary_txt from ch_guppy_pycoqc_summary

    output:
    path "*.html" 
    path "*.json" into ch_pycoqc_multiqc
    path "v_pycoqc.txt" into ch_pycoqc_version

    script:
    """
    pycoQC -f $summary_txt -o pycoqc.html -j pycoqc.json
    pycoQC --version &> v_pycoqc.txt
    """
}

/*
 * QC using NanoPlot
 */
process NANOPLOT_SUMMARY {
    tag "$summary_txt"
    label 'process_low'
    publishDir "${params.outdir}/nanoplot/summary", mode: params.publish_dir_mode

    when:
    !params.skip_basecalling && !params.skip_qc && !params.skip_nanoplot

    input:
    path summary_txt from ch_guppy_nanoplot_summary

    output:
    path "*.{png,html,txt,log}"

    script:
    """
    NanoPlot -t $task.cpus --summary $summary_txt
    """
}

/*
 * FastQ QC using NanoPlot
 */
process NANOPLOT_FASTQ {
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}/nanoplot/fastq/${sample}", mode: params.publish_dir_mode

    when:
    !params.skip_qc && !params.skip_nanoplot

    input:
    tuple val(sample), path(fastq) from ch_fastq_nanoplot.map { ch -> [ ch[0], ch[1] ] }

    output:
    path "*.{png,html,txt,log}"

    script:
    """
    NanoPlot -t $task.cpus --fastq $fastq
    """
}

/*
 * FastQ QC using FastQC
 */
process FASTQC {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    tuple val(sample), path(fastq) from ch_fastq_fastqc.map { ch -> [ ch[0], ch[1] ] }

    output:
    path "*.{zip,html}" into ch_fastqc_multiqc

    script:
    """
    [ ! -f  ${sample}.fastq.gz ] && ln -s $fastq ${sample}.fastq.gz

    fastqc \\
        -q \\
        -t $task.cpus \\
        ${sample}.fastq.gz
    """
}

if (!params.skip_alignment) {

    // Get unique list of all fasta files
    ch_fastq_sizes
        .filter { it[2] }
        .map { it -> [ it[2], it[-1].toString() ] }  // [ fasta, annotation_str ]
        .unique()
        .set { ch_fastq_sizes }

    /*
     * Make chromosome sizes file
     */
    process GET_CHROM_SIZES {
        tag "$fasta"

        input:
        tuple path(fasta), val(name) from ch_fastq_sizes

        output:
        tuple path("*.sizes"), val(name) into ch_chrom_sizes

        script:
        """
        samtools faidx $fasta
        cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
        """
    }

    // Get unique list of all gtf files
    ch_fastq_gtf
        .filter { it[3] }
        .map { it -> [ it[3], it[-1] ] }  // [ gtf, annotation_str ]
        .unique()
        .set { ch_fastq_gtf }

    /*
     * Convert GTF to BED12
     */
    process GTF_TO_BED {
        tag "$gtf"
        label 'process_low'

        input:
        tuple path(gtf), val(name) from ch_fastq_gtf

        output:
        tuple path("*.bed"), val(name) into ch_gtf_bed

        script: // This script is bundled with the pipeline, in nf-core/nanoseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }

    ch_chrom_sizes
        .join(ch_gtf_bed, by: 1, remainder:true)
        .map { it -> [ it[1], it[2], it[0] ] }
        .cross(ch_fastq_index) { it -> it[-1] }
        .flatten()
        .collate(9)
        .map { it -> [ it[5], it[0], it[6], it[1], it[7], it[8] ]} // [ fasta, sizes, gtf, bed, is_transcripts, annotation_str ]
        .unique()
        .set { ch_fasta_index }

    /*
     * Create genome/transcriptome index
     */
    if (params.aligner == 'minimap2') {
        process MINIMAP2_INDEX {
            tag "$fasta"
            label 'process_medium'

            input:
            tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), val(annotation_str) from ch_fasta_index

            output:
            tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path("*.mmi"), val(annotation_str) into ch_index

            script:
            preset    = (params.protocol == 'DNA' || is_transcripts) ? "-ax map-ont" : "-ax splice"
            kmer      = (params.protocol == 'directRNA') ? "-k14" : ""
            stranded  = (params.stranded || params.protocol == 'directRNA') ? "-uf" : ""
            // TODO pipeline: Should be staging bed file properly as an input
            junctions = (params.protocol != 'DNA' && bed) ? "--junc-bed ${file(bed)}" : ""
            """
            minimap2 \\
                $preset \\
                $kmer \\
                $stranded \\
                $junctions \\
                -t $task.cpus \\
                -d ${fasta}.mmi \\
                $fasta
            """
        }
    } else {
        process GRAPHMAP2_INDEX {
            tag "$fasta"
            label 'process_high'

            input:
            tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), val(annotation_str) from ch_fasta_index

            output:
            tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path("*.gmidx"), val(annotation_str) into ch_index

            script:
            preset = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
            // TODO pipeline: Should be staging gtf file properly as an input
            junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf ${file(gtf)}" : ""
            """
            graphmap2 \\
                align \\
                $preset \\
                $junctions \\
                -t $task.cpus \\
                -I \\
                -r $fasta
            """
        }
    }

    ch_index
        .cross(ch_fastq_align) { it -> it[-1] }
        .flatten()
        .collate(13)
        .map { it -> [ it[7], it[8], it[0], it[1], it[2], it[3], it[4], it[5] ] } // [ sample, fastq, fasta, sizes, gtf, bed, is_transcripts, index ]
        .set { ch_index }

    /*
     * Align fastq files
     */
    if (params.aligner == 'minimap2') {
        process MINIMAP2_ALIGN {
            tag "$sample"
            label 'process_medium'
            if (params.save_align_intermeds) {
                publishDir path: "${params.outdir}/${params.aligner}", mode: params.publish_dir_mode,
                    saveAs: { filename ->
                                  if (filename.endsWith(".sam")) filename
                            }
            }

            input:
            tuple val(sample), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index) from ch_index


            output:
            tuple val(sample), path(sizes), val(is_transcripts), path("*.sam") into ch_align_sam

            script:
            preset    = (params.protocol == 'DNA' || is_transcripts) ? "-ax map-ont" : "-ax splice"
            kmer      = (params.protocol == 'directRNA') ? "-k14" : ""
            stranded  = (params.stranded || params.protocol == 'directRNA') ? "-uf" : ""
            // TODO pipeline: Should be staging bed file properly as an input
            junctions = (params.protocol != 'DNA' && bed) ? "--junc-bed ${file(bed)}" : ""
            """
            minimap2 \\
                $preset \\
                $kmer \\
                $stranded \\
                $junctions \\
                -t $task.cpus \\
                $index \\
                $fastq > ${sample}.sam
            """
        }
    } else {
        process GRAPHMAP2_ALIGN {
            tag "$sample"
            label 'process_medium'
            if (params.save_align_intermeds) {
                publishDir path: "${params.outdir}/${params.aligner}", mode: params.publish_dir_mode,
                    saveAs: { filename ->
                                  if (filename.endsWith(".sam")) filename
                            }
            }

            input:
            tuple val(sample), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index) from ch_index

            output:
            tuple val(sample), path(sizes), val(is_transcripts), path("*.sam") into ch_align_sam

            script:
            preset    = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
            // TODO pipeline: Should be staging gtf file properly as an input
            junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf ${file(gtf)}" : ""
            """
            graphmap2 \\
                align \\
                $preset \\
                $junctions \\
                -t $task.cpus \\
                -r $fasta \\
                -i $index \\
                -d $fastq \\
                -o ${sample}.sam \\
                --extcigar
            """
        }
    }

    /*
    * Coordinate sort BAM files
    */
    process SAMTOOLS_SORT {
        tag "$sample"
        label 'process_medium'
        publishDir path: "${params.outdir}/${params.aligner}", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                        else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                        else if (filename.endsWith(".sorted.bam")) filename
                        else if (filename.endsWith(".sorted.bam.bai")) filename
                        else null
                    }

        input:
        tuple val(sample), path(sizes), val(is_transcripts), path(sam) from ch_align_sam

        output:
        tuple val(sample), path(sizes), val(is_transcripts), path("*.sorted.bam"), path("*.sorted.bam.bai") into ch_sortbam_bedgraph,
                                                                                                                 ch_sortbam_bed12,
                                                                                                                 ch_sortbam_quant
        path "*.{flagstat,idxstats,stats}" into ch_sortbam_stats_multiqc
        
        script:
        """
        samtools view -b -h -O BAM -@ $task.cpus -o ${sample}.bam $sam
        samtools sort -@ $task.cpus -o ${sample}.sorted.bam -T $sample ${sample}.bam
        samtools index ${sample}.sorted.bam
        samtools flagstat ${sample}.sorted.bam > ${sample}.sorted.bam.flagstat
        samtools idxstats ${sample}.sorted.bam > ${sample}.sorted.bam.idxstats
        samtools stats ${sample}.sorted.bam > ${sample}.sorted.bam.stats
        """
    }

    ch_sortbam_quant
        .map { it -> [ it[0], it[3] ] }
        .set { ch_sortbam_quant }

} else {
    ch_sortbam_bedgraph      = Channel.empty()
    ch_sortbam_bed12         = Channel.empty()
    ch_sortbam_stats_multiqc = Channel.empty()

    ch_samplesheet_reformat
        .splitCsv(header:true, sep:',')
        .map { get_sample_info(it, params.genomes) }
        .map { it -> if (it[1].toString().endsWith('.bam')) [ it[0] , it[1] ] }
        .set { ch_sortbam_rename }

    /*
     * Rename BAM inputs to <GROUP_REPLICATE>
     */
    process BAM_RENAME {
        tag "$sample"
        
        input:
        tuple val(sample), path(bam) from ch_sortbam_rename
        
        output:
        tuple val(sample), path("*.bam") into ch_sortbam_quant
        
        script:
        """
        [ ! -f ${sample}.bam ] && ln -s $bam ${sample}.bam
        """
    }
}


/*
 * Convert BAM to BEDGraph
 */
process BEDTOOLS_GENOMECOV {
    tag "$sample"
    label 'process_medium'

    when:
    !params.skip_alignment && !params.skip_bigwig

    input:
    tuple val(sample), path(sizes), val(is_transcripts), path(bam), path(bai) from ch_sortbam_bedgraph

    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.bedGraph") into ch_bedgraph

    script:
    split = (params.protocol == 'DNA' || is_transcripts) ? "" : "-split"
    """
    bedtools \\
        genomecov \\
        $split \\
        -ibam ${bam[0]} \\
        -bg \\
        | bedtools sort > ${sample}.bedGraph
    """
}

/*
 * Convert BEDGraph to BigWig
 */
process UCSC_BEDGRAPHTOBIGWIG {
    tag "$sample"
    label 'process_medium'
    publishDir path: "${params.outdir}/${params.aligner}/bigwig/", mode: params.publish_dir_mode

    when:
    !params.skip_alignment && !params.skip_bigwig

    input:
    tuple val(sample), path(sizes), val(is_transcripts), path(bedgraph) from ch_bedgraph

    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.bigWig") into ch_bigwig

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${sample}.bigWig
    """
}

/*
 * Convert BAM to BED12
 */
process BEDTOOLS_BAMTOBED {
    tag "$sample"
    label 'process_medium'

    when:
    !params.skip_alignment && !params.skip_bigbed && (params.protocol == 'directRNA' || params.protocol == 'cDNA')

    input:
    tuple val(sample), path(sizes), val(is_transcripts), path(bam), path(bai) from ch_sortbam_bed12

    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.bed12") into ch_bed12

    script:
    """
    bedtools \\
        bamtobed \\
        -bed12 \\
        -cigar \\
        -i ${bam[0]} \\
        | bedtools sort > ${sample}.bed12
    """
}

/*
 * Convert BED12 to BigBED
 */
process UCSC_BED12TOBIGBED {
    tag "$sample"
    label 'process_medium'
    publishDir path: "${params.outdir}/${params.aligner}/bigbed/", mode: params.publish_dir_mode

    when:
    !params.skip_alignment && !params.skip_bigbed && (params.protocol == 'directRNA' || params.protocol == 'cDNA')

    input:
    tuple val(sample), path(sizes), val(is_transcripts), path(bed12) from ch_bed12

    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.bigBed") into ch_bigbed

    script:
    """
    bedToBigBed \\
        $bed12 \\
        $sizes \\
        ${sample}.bigBed
    """
}

if (!params.skip_quantification && (params.protocol == 'cDNA' || params.protocol == 'directRNA')) {
    if (params.quantification_method == 'bambu') {
        
        ch_sample_annotation
            .map { it -> [ it[2], it[3] ] }
            .unique()
            .set { ch_sample_annotation }

        /*
         * Quantification and novel isoform detection with bambu
         */
        process BAMBU {
            label 'process_high'
            publishDir "${params.outdir}/${params.quantification_method}", mode: params.publish_dir_mode,
                saveAs: { filename ->
                            if (!filename.endsWith("bambu.txt")) filename
                        }

            input:
            tuple path(fasta), path(gtf) from ch_sample_annotation
            path bams from ch_sortbam_quant.collect{ it[1] }
            
            output:
            path "counts_gene.txt" into ch_gene_counts
            path "counts_transcript.txt" into ch_transcript_counts
            path "extended_annotations.gtf" 
            
            script:
            """
            run_bambu.r \\
                --tag=. \\
                --ncore=$task.cpus \\
                --annotation=$gtf \\
                --fasta=$fasta $bams
            """
        }
        ch_featurecounts_transcript_multiqc = Channel.empty()
        ch_featurecounts_gene_multiqc       = Channel.empty()

    } else {
        
        ch_sample_annotation
            .map { it -> [ it[0], it[2], it[3] ] }
            .join(ch_sortbam_quant)
            .into { 
                ch_sample_annotation
                ch_sample_gtf
                ch_sample_featurecounts
            }
           
        /*
         * Novel isoform detection with StringTie
         */
        process STRINGTIE2 {
            tag "$sample"
            label 'process_medium'
            publishDir "${params.outdir}/${params.quantification_method}", mode: params.publish_dir_mode

            input:
            tuple val(sample), path(fasta), path(gtf), path(bam) from ch_sample_annotation

            output:
            path "*.stringtie.gtf" into ch_stringtie_gtf

            script:
            """
            stringtie \\
                -L \\
                -G $gtf \\
                -o ${sample}.stringtie.gtf $bam
            """
        }

        ch_sample_gtf
            .map { it -> [ it[2]] }
            .unique()
            .set { ch_sample_gtf }

        /*
         * Merge isoforms across samples caleed by StringTie
         */
        process STRINGTIE2_MERGE {
            label 'process_medium'
            publishDir "${params.outdir}/${params.quantification_method}", mode: params.publish_dir_mode

            input:
            path gtf from ch_sample_gtf
            path gtfs from ch_stringtie_gtf.collect()
            
            output:
            path "stringtie.merged.gtf" into ch_stringtie_merged_gtf

            script:
            """
            stringtie \\
                --merge $gtfs \\
                -G $gtf \\
                -o stringtie.merged.gtf
            """
        }

        /*
         * Gene and transcript quantification with featureCounts
         */
        process SUBREAD_FEATURECOUNTS {
            label 'process_medium'
            publishDir "${params.outdir}/${params.quantification_method}", mode: params.publish_dir_mode

            input:
            path gtf from ch_stringtie_merged_gtf
            path bams from ch_sample_featurecounts.collect{ it[-1] }

            output:
            path "counts_gene.txt" into ch_gene_counts
            path "counts_transcript.txt" into ch_transcript_counts
            path "counts_gene.txt.summary" into ch_featurecounts_gene_multiqc
            path "counts_transcript.txt.summary" into ch_featurecounts_transcript_multiqc
            
            script:
            """
            featureCounts \\
                -L \\
                -O \\
                -f \\
                -g gene_id \\
                -t exon \\
                -T $task.cpus \\
                -a $gtf \\
                -o counts_gene.txt \\
                $bams
            
            featureCounts \\
                -L \\
                -O \\
                -f \\
                --primary \\
                --fraction \\
                -F GTF \\
                -g transcript_id \\
                -t transcript \\
                --extraAttributes gene_id \\
                -T $task.cpus \\
                -a $gtf \\
                -o counts_transcript.txt \\
                $bams
            """
        }
    }

    if (!params.skip_differential_analysis) {

        /*
         * DESeq2 differential expression of genes
         */
        process DESEQ2 {
            label 'process_medium'
            publishDir "${params.outdir}/${params.quantification_method}/deseq2", mode: params.publish_dir_mode

            when:
            MULTIPLE_CONDITIONS && REPLICATES_EXIST

            input:
            path counts from ch_gene_counts
            
            output:
            path "*.txt"

            script:
            """
            run_deseq2.r $params.quantification_method $counts
            """
        }

        /*
         * DEXseq differential expression of transcripts
         */
        process DEXSEQ {
            label 'process_medium'
            publishDir "${params.outdir}/${params.quantification_method}/dexseq", mode: params.publish_dir_mode

            when:
            MULTIPLE_CONDITIONS && REPLICATES_EXIST

            input:
            path counts from ch_transcript_counts

            output:
            path "*.txt"
            
            script:
            """
            run_dexseq.r $params.quantification_method $counts
            """
        }
    }
} else {
    ch_featurecounts_transcript_multiqc = Channel.empty()
    ch_featurecounts_gene_multiqc       = Channel.empty()
}

/*
 * Output Description HTML
 */
process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path output_docs from ch_output_docs
    path images from ch_output_docs_images

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Parse software version numbers
 */
process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    input:
    path guppy from ch_guppy_version.collect().ifEmpty([])
    path pycoqc from ch_pycoqc_version.collect().ifEmpty([])

    output:
    path 'software_versions_mqc.yaml' into software_versions_yaml
    path "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    qcat --version &> v_qcat.txt
    NanoLyse --version &> v_nanolyse.txt
    NanoPlot --version &> v_nanoplot.txt
    fastqc --version > v_fastqc.txt
    minimap2 --version &> v_minimap2.txt
    echo \$(graphmap2 2>&1) > v_graphmap2.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt
    stringtie --version > v_stringtie.txt
    echo \$(featureCounts -v 2>&1) > v_featurecounts.txt
    echo \$(R --version 2>&1) > v_r.txt
    Rscript -e "library(bambu); write(x=as.character(packageVersion('bambu')), file='v_bambu.txt')"
    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='v_deseq2.txt')"
    Rscript -e "library(DRIMSeq); write(x=as.character(packageVersion('DRIMSeq')), file='v_drimseq.txt')"
    Rscript -e "library(DEXSeq); write(x=as.character(packageVersion('DEXSeq')), file='v_dexseq.txt')"
    Rscript -e "library(stageR); write(x=as.character(packageVersion('stageR')), file='v_stager.txt')"
    Rscript -e "library(BSgenome); write(x=as.character(packageVersion('BSgenome')), file='v_bsgenome.txt')"
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-nanoseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/nanoseq Workflow Summary'
    section_href: 'https://github.com/nf-core/nanoseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * MultiQC
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    when:
    !params.skip_multiqc

    input:
    path (multiqc_config) from ch_multiqc_config
    path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    
    path ('pycoqc/*') from ch_pycoqc_multiqc.collect().ifEmpty([])
    path ('fastqc/*') from ch_fastqc_multiqc.collect().ifEmpty([])
    path ('samtools/*') from ch_sortbam_stats_multiqc.collect().ifEmpty([])
    path ('featurecounts/gene/*') from ch_featurecounts_gene_multiqc.collect().ifEmpty([])
    path ('featurecounts/transcript/*') from ch_featurecounts_transcript_multiqc.collect().ifEmpty([])
    path ('software_versions/*') from software_versions_yaml.collect()
    path workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    path "*multiqc_report.html" into ch_multiqc_report
    path "*_data"
    path "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/nanoseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/nanoseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/nanoseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/nanoseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def max_multiqc_email_size = params.max_multiqc_email_size as nextflow.util.MemoryUnit
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/nanoseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/nanoseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/nanoseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/nanoseq]${c_red} Pipeline completed with errors${c_reset}-"
    }
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/nanoseq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
