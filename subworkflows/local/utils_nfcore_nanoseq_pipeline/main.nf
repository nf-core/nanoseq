//
// Subworkflow with functionality specific to the nf-core/nanoseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/nanoseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/nanoseq/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Samplesheet contents are validated against assets/schema_input.json by nf-schema
    // (see the `schema` key of `--input` in nextflow_schema.json). Paths in the samplesheet
    // are resolved by INPUT_CHECK because raw data directories may only exist at runtime.
    //
    def ch_samplesheet = channel.value(file(input, checkIfExists: true))

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    if (!params.fasta) {
        error("Reference FASTA not specified! Please provide --fasta or a valid --genome.")
    }

    def rna_protocols = ['cDNA', 'directRNA']
    if (!params.skip_quantification && !rna_protocols.contains(params.protocol)) {
        error("Invalid protocol option for quantification: ${params.protocol}. Valid options: ${rna_protocols.join(', ')}. Use --skip_quantification for DNA data.")
    }
    if (!params.skip_quantification && !params.gtf) {
        error("A GTF annotation is required for quantification. Please provide --gtf, a --genome with a GTF, or use --skip_quantification.")
    }

    if (!params.skip_demultiplexing) {
        def qcat_barcode_kits = [
            'Auto', 'RBK001', 'RBK004', 'NBD103/NBD104', 'NBD114', 'NBD104/NBD114', 'PBC001',
            'PBC096', 'RPB004/RLB001', 'PBK004/LWB001', 'RAB204', 'VMK001', 'DUAL'
        ]
        if (params.barcode_kit && !qcat_barcode_kits.contains(params.barcode_kit)) {
            error("Invalid barcode kit for demultiplexing with qcat: ${params.barcode_kit}. Valid options: ${qcat_barcode_kits.join(', ')}")
        }
        if (!params.input_path) {
            error("Please specify --input_path to perform demultiplexing!")
        }
    }

    if (!params.skip_basecalling && !params.dorado_model) {
        error("Please specify --dorado_model to perform basecalling, or use --skip_basecalling.")
    }

    if (params.run_nanolyse && !params.nanolyse_fasta && env('NXF_OFFLINE') == 'true') {
        error("Running offline so cannot download the lambda FASTA for NanoLyse! Please specify --nanolyse_fasta.")
    }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            !params.skip_basecalling           ? "dorado," : "",
            !params.skip_demultiplexing        ? "qcat," : "",
            params.run_nanolyse                ? "NanoLyse (De Coster et al. 2018)," : "",
            !params.skip_qc && !params.skip_nanoplot  ? "NanoPlot (De Coster et al. 2018)," : "",
            !params.skip_qc && !params.skip_toulligqc ? "ToulligQC," : "",
            !params.skip_qc && !params.skip_fastqc    ? "FastQC (Andrews 2010)," : "",
            !params.skip_alignment ? (params.aligner == 'minimap2' ? "minimap2 (Li 2018)," : "GraphMap2 (Sović et al. 2016),") : "",
            "SAMtools (Li et al. 2009),",
            !params.skip_alignment && !(params.skip_bigwig && params.skip_bigbed) ? "BEDTools (Quinlan and Hall 2010), UCSC tools (Kent et al. 2010)," : "",
            !params.skip_quantification ? (params.quantification_method == 'bambu' ? "bambu (Chen et al. 2023)," : "StringTie2 (Kovaka et al. 2019), featureCounts (Liao et al. 2014),") : "",
            !params.skip_quantification && !params.skip_differential_analysis ? "DESeq2 (Love et al. 2014), DEXSeq (Anders et al. 2012)," : "",
            !params.skip_modification_analysis && !params.skip_xpore  ? "xPore (Pratanwanich et al. 2021)," : "",
            !params.skip_modification_analysis && !params.skip_m6anet ? "m6anet (Hendra et al. 2022)," : "",
            !params.skip_fusion_analysis ? "JAFFAL (Davidson et al. 2022)," : "",
            params.bedmethyl_out ? "modkit," : "",
            "and MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim().replaceAll(/\s+/, ' ')

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            params.run_nanolyse || (!params.skip_qc && !params.skip_nanoplot) ? "<li>De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-2669. doi: 10.1093/bioinformatics/bty149</li>" : "",
            !params.skip_qc && !params.skip_fastqc ? "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>" : "",
            !params.skip_alignment && params.aligner == 'minimap2'  ? "<li>Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191</li>" : "",
            !params.skip_alignment && params.aligner == 'graphmap2' ? "<li>Sović I, Šikić M, Wilm A, Fenlon SN, Chen S, Nagarajan N. Fast and sensitive mapping of nanopore sequencing reads with GraphMap. Nat Commun. 2016;7:11307. doi: 10.1038/ncomms11307</li>" : "",
            "<li>Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-2079. doi: 10.1093/bioinformatics/btp352</li>",
            !params.skip_alignment && !(params.skip_bigwig && params.skip_bigbed) ? "<li>Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842. doi: 10.1093/bioinformatics/btq033</li><li>Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010;26(17):2204-2207. doi: 10.1093/bioinformatics/btq351</li>" : "",
            !params.skip_quantification && params.quantification_method == 'bambu' ? "<li>Chen Y, Sim A, Wan YK, et al. Context-aware transcript quantification from long-read RNA-seq data with Bambu. Nat Methods. 2023;20:1187-1195. doi: 10.1038/s41592-023-01908-w</li>" : "",
            !params.skip_quantification && params.quantification_method == 'stringtie2' ? "<li>Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biol. 2019;20(1):278. doi: 10.1186/s13059-019-1910-1</li><li>Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-930. doi: 10.1093/bioinformatics/btt656</li>" : "",
            !params.skip_quantification && !params.skip_differential_analysis ? "<li>Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi: 10.1186/s13059-014-0550-8</li><li>Anders S, Reyes A, Huber W. Detecting differential usage of exons from RNA-seq data. Genome Res. 2012;22(10):2008-2017. doi: 10.1101/gr.133744.111</li>" : "",
            !params.skip_modification_analysis && !params.skip_xpore ? "<li>Pratanwanich PN, Yao F, Chen Y, et al. Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. Nat Biotechnol. 2021;39(11):1394-1402. doi: 10.1038/s41587-021-00949-w</li>" : "",
            !params.skip_modification_analysis && !params.skip_m6anet ? "<li>Hendra C, Pratanwanich PN, Wan YK, Goh WSS, Thiery A, Göke J. Detection of m6A from direct RNA sequencing using a multiple instance learning framework. Nat Methods. 2022;19:1590-1598. doi: 10.1038/s41592-022-01666-1</li>" : "",
            !params.skip_fusion_analysis ? "<li>Davidson NM, Chen Y, Sadras T, et al. JAFFAL: detecting fusion genes with long-read transcriptome sequencing. Genome Biol. 2022;23(1):10. doi: 10.1186/s13059-021-02588-5</li>" : "",
            "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
