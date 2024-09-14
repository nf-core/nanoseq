from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Union

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)


@dataclass
class SampleSheet:
    group: str
    replicate: int
    barcode: Optional[int]
    input_file: Optional[LatchFile]
    fasta: Optional[Union[LatchFile, str]]
    gtf: Optional[LatchFile]


class Protocol(Enum):
    cDNA = "cDNA"
    DNA = "DNA"
    directRNA = "directRNA"


class BarcodeKit(Enum):
    Auto = "Auto"
    RBK001 = "RBK001"
    RBK004 = "RBK004"
    NBD103_104 = "NBD103/NBD104"
    NBD114 = "NBD114"
    NBD104_114 = "NBD104/NBD114"
    PBC001 = "PBC001"
    PBC096 = "PBC096"
    RPB004_RLB001 = "RPB004/RLB001"
    RPB004_LWB001 = "RPB004/LWB001"
    RAB204 = "RAB204"
    VMK001 = "VMK001"


flow = [
    Section(
        "Inputs",
        Params(
            "input",
            "protocol",
        ),
    ),
    Section(
        "Advanced Demultiplexing",
        Fork(
            "input_path_source",
            "",
            file=ForkBranch("File", Params("input_path_file")),
            directory=ForkBranch("Directory", Params("input_path_dir")),
        ),
        Params(
            "skip_demultiplexing",
            "barcode_kit",
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Advanced Demultiplexing",
        Params(
            "barcode_both_ends",
            "trim_barcodes",
            "gpu_device",
            "gpu_cluster_options",
            "qcat_min_score",
            "qcat_detect_middle",
            "run_nanolyse",
            "nanolyse_fasta",
        ),
    ),
    Spoiler(
        "Alignment",
        Params(
            "aligner",
            "stranded",
            "save_align_intermeds",
            "skip_alignment",
        ),
    ),
    Spoiler(
        "Variant Calling",
        Params(
            "call_variants",
            "variant_caller",
            "structural_variant_caller",
            "split_mnps",
            "deepvariant_gpu",
            "phase_vcf",
            "skip_vc",
            "skip_sv",
        ),
    ),
    Spoiler(
        "Differential Analysis",
        Params(
            "quantification_method",
            "skip_quantification",
            "skip_differential_analysis",
        ),
    ),
    Spoiler(
        "RNA Fusion Analysis",
        Params(
            "jaffal_ref_dir",
            "skip_fusion_analysis",
        ),
    ),
    Spoiler(
        "RNA Modification Analysis",
        Params(
            "skip_modification_analysis",
            "skip_xpore",
            "skip_m6anet",
        ),
    ),
    Spoiler(
        "General Options",
        Params(
            "email",
            "multiqc_title",
            "skip_bigbed",
            "skip_bigwig",
            "skip_nanoplot",
            "skip_fastqc",
            "skip_multiqc",
            "skip_qc",
        ),
    ),
]


generated_parameters = {
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "input": NextflowParameter(
        type=List[SampleSheet],
        display_name="Samplesheet",
        samplesheet=True,
        samplesheet_type="csv",
        description="Samplesheet containing information about the samples in the experiment.",
    ),
    "protocol": NextflowParameter(
        type=Protocol,
        display_name="Protocol",
        default=Protocol.cDNA,
        description="Input sample type. Valid options: 'DNA', 'cDNA', and 'directRNA'.",
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        display_name="Output Directory",
        default=None,
        description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
    ),
    "email": NextflowParameter(
        type=Optional[str],
        display_name="Email",
        default=None,
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Title",
        default=None,
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "input_path_file": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Input Basecalled File",
        default=None,
        description="Path to Nanopore basecalled fastq file that requires demultiplexing.",
    ),
    "input_path_dir": NextflowParameter(
        type=Optional[LatchDir],
        display_name="Input Run Directory",
        default=None,
        description="Path to Nanopore run directory files (e.g. 'fastq_pass/*') that requires demultiplexing.",
    ),
    "barcode_kit": NextflowParameter(
        type=Optional[BarcodeKit],
        display_name="Barcode Kit",
        default=None,
        description="Barcode kit used to perform the sequencing e.g. 'SQK-PBK004'.",
    ),
    "barcode_both_ends": NextflowParameter(
        type=bool,
        display_name="Barcode Both Ends",
        default=None,
        description="Require barcode on both ends for basecaller.",
    ),
    "trim_barcodes": NextflowParameter(
        type=bool,
        display_name="Trim Barcodes",
        default=None,
        description="Trim barcodes from the output sequences in the FastQ files from basecaller.",
    ),
    "gpu_device": NextflowParameter(
        type=Optional[str],
        display_name="GPU Device",
        default="auto",
        description="Device specified in GPU mode using '--device'.",
    ),
    "gpu_cluster_options": NextflowParameter(
        type=Optional[str],
        display_name="GPU Cluster Options",
        default=None,
        description="Cluster options required to use GPU resources (e.g. '--part=gpu --gres=gpu:1').",
    ),
    "qcat_min_score": NextflowParameter(
        type=Optional[int],
        display_name="Qcat Min Score",
        default=60,
        description="Specify the minimum quality score for qcat in the range 0-100.",
    ),
    "qcat_detect_middle": NextflowParameter(
        type=bool,
        display_name="Qcat Detect Middle",
        default=None,
        description="Search for adapters in the whole read by applying the '--detect-middle' parameter in qcat.",
    ),
    "skip_demultiplexing": NextflowParameter(
        type=bool,
        display_name="Skip Demultiplexing",
        default=None,
        description="Skip demultiplexing with qcat.",
    ),
    "run_nanolyse": NextflowParameter(
        type=bool,
        display_name="Run NanoLyse",
        default=None,
        description="Filter reads from FastQ files using NanoLyse",
    ),
    "nanolyse_fasta": NextflowParameter(
        type=Optional[str],
        display_name="NanoLyse FASTA",
        default=None,
        description="Fasta file to be filtered against using NanoLyse",
    ),
    "aligner": NextflowParameter(
        type=Optional[str],
        display_name="Aligner",
        default="minimap2",
        description="Specifies the aligner to use i.e. 'minimap2' or 'graphmap2'.",
    ),
    "stranded": NextflowParameter(
        type=bool,
        display_name="Stranded",
        default=None,
        description="Specifies if the data is strand-specific. Automatically activated when using '--protocol directRNA'.",
    ),
    "save_align_intermeds": NextflowParameter(
        type=bool,
        display_name="Save Alignment Intermediates",
        default=None,
        description="Save the '.sam' files from the alignment step - not done by default.",
    ),
    "skip_alignment": NextflowParameter(
        type=bool,
        display_name="Skip Alignment",
        default=None,
        description="Skip alignment and downstream processes.",
    ),
    "call_variants": NextflowParameter(
        type=bool,
        display_name="Call Variants",
        default=None,
        description="Specifies if variant calling will executed.",
    ),
    "variant_caller": NextflowParameter(
        type=Optional[str],
        display_name="Variant Caller",
        default="medaka",
        description="Specifies the variant caller that will used to call small variants (available are: medaka, deepvariant, and pepper_margin_deepvariant). Variant calling is only available if '--call_variants' is set and the protocol is set to `DNA`. Please note `deepvariant` and `pepper_margin_deepvariant` are only avaible if using singularity or docker.",
    ),
    "structural_variant_caller": NextflowParameter(
        type=Optional[str],
        display_name="Structural Variant Caller",
        default="sniffles",
        description="Specifies the variant caller that will be used to call structural variants (available are: sniffles and cutesv). Structural variant calling is only available if '--call_variants' is set and the protocol is set to `DNA`.",
    ),
    "split_mnps": NextflowParameter(
        type=bool,
        display_name="Split MNPs",
        default=None,
        description="Specifies if MNPs will be split into SNPs when using medaka.",
    ),
    "deepvariant_gpu": NextflowParameter(
        type=bool,
        display_name="DeepVariant GPU",
        default=None,
        description="Specifies whether to call variants with pepper_margin_deepvariant in GPU mode.",
    ),
    "phase_vcf": NextflowParameter(
        type=bool,
        display_name="Phase VCF",
        default=None,
        description="Specifies if vcf will be phased when using medaka.",
    ),
    "skip_vc": NextflowParameter(
        type=bool,
        display_name="Skip Variant Calling",
        default=None,
        description="Skip variant calling.",
    ),
    "skip_sv": NextflowParameter(
        type=bool,
        display_name="Skip Structural Variant Calling",
        default=None,
        description="Skip structural variant calling.",
    ),
    "quantification_method": NextflowParameter(
        type=Optional[str],
        display_name="Quantification Method",
        default="bambu",
        description="Specifies the transcript quantification method to use (available are: bambu or stringtie2). Only available when protocol is cDNA or directRNA.",
    ),
    "skip_quantification": NextflowParameter(
        type=bool,
        display_name="Skip Quantification",
        default=None,
        description="Skip transcript quantification and differential analysis.",
    ),
    "skip_differential_analysis": NextflowParameter(
        type=bool,
        display_name="Skip Differential Analysis",
        default=None,
        description="Skip differential analysis with DESeq2 and DEXSeq.",
    ),
    "jaffal_ref_dir": NextflowParameter(
        type=Optional[str],
        display_name="JAFFAL Reference Directory",
        default="for_jaffal",
        description="Specifies the reference directory for JAFFAL.",
    ),
    "skip_fusion_analysis": NextflowParameter(
        type=bool,
        display_name="Skip Fusion Analysis",
        default=None,
        description="Skip differential analysis with DESeq2 and DEXSeq.",
    ),
    "skip_modification_analysis": NextflowParameter(
        type=bool,
        display_name="Skip Modification Analysis",
        default=None,
        description="Skip RNA modification analysis.",
    ),
    "skip_xpore": NextflowParameter(
        type=bool,
        display_name="Skip xPore",
        default=None,
        description="Skip differential modification analysis with xpore.",
    ),
    "skip_m6anet": NextflowParameter(
        type=bool,
        display_name="Skip m6anet",
        default=None,
        description="Skip m6A detection with m6anet.",
    ),
    "skip_bigbed": NextflowParameter(
        type=bool,
        display_name="Skip BigBed",
        default=None,
        description="Skip BigBed file generation.",
    ),
    "skip_bigwig": NextflowParameter(
        type=bool,
        display_name="Skip BigWig",
        default=None,
        description="Skip BigWig file generation.",
    ),
    "skip_nanoplot": NextflowParameter(
        type=bool,
        display_name="Skip NanoPlot",
        default=None,
        description="Skip NanoPlot.",
    ),
    "skip_fastqc": NextflowParameter(
        type=bool,
        display_name="Skip FastQC",
        default=None,
        description="Skip FastQC.",
    ),
    "skip_multiqc": NextflowParameter(
        type=bool,
        display_name="Skip MultiQC",
        default=None,
        description="Skip MultiQC.",
    ),
    "skip_qc": NextflowParameter(
        type=bool,
        display_name="Skip QC",
        default=None,
        description="Skip all QC steps apart from MultiQC.",
    ),
}
