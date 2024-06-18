
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

generated_parameters = {
    'input': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title='Input/output options',
        description='Path to comma-separated file containing information about the samples in the experiment.',
    ),
    'protocol': NextflowParameter(
        type=str,
        default=None,
        section_title=None,
        description="Input sample type. Valid options: 'DNA', 'cDNA',  and 'directRNA'.",
    ),
    'outdir': NextflowParameter(
        type=typing.Optional[typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})]],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'email': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Email address for completion summary.',
    ),
    'multiqc_title': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='MultiQC report title. Printed as page header, used for filename if not otherwise specified.',
    ),
    'input_path': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title='Demultiplexing options',
        description="Path to Nanopore run directory files (e.g. 'fastq_pass/*') or a basecalled fastq file that requires demultiplexing.",
    ),
    'barcode_kit': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Barcode kit used to perform the sequencing e.g. 'SQK-PBK004'.",
    ),
    'barcode_both_ends': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Require barcode on both ends for basecaller.',
    ),
    'trim_barcodes': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Trim barcodes from the output sequences in the FastQ files from basecaller.',
    ),
    'gpu_device': NextflowParameter(
        type=typing.Optional[str],
        default='auto',
        section_title=None,
        description="Device specified in GPU mode using '--device'.",
    ),
    'gpu_cluster_options': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description="Cluster options required to use GPU resources (e.g. '--part=gpu --gres=gpu:1').",
    ),
    'qcat_min_score': NextflowParameter(
        type=typing.Optional[int],
        default=60,
        section_title=None,
        description='Specify the minimum quality score for qcat in the range 0-100.',
    ),
    'qcat_detect_middle': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Search for adapters in the whole read by applying the '--detect-middle' parameter in qcat.",
    ),
    'skip_demultiplexing': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip demultiplexing with qcat.',
    ),
    'run_nanolyse': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Filter reads from FastQ files using NanoLyse',
    ),
    'nanolyse_fasta': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Fasta file to be filtered against using NanoLyse',
    ),
    'aligner': NextflowParameter(
        type=typing.Optional[str],
        default='minimap2',
        section_title='Alignment options',
        description="Specifies the aligner to use i.e. 'minimap2' or 'graphmap2'.",
    ),
    'stranded': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Specifies if the data is strand-specific. Automatically activated when using '--protocol directRNA'.",
    ),
    'save_align_intermeds': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description="Save the '.sam' files from the alignment step - not done by default.",
    ),
    'skip_alignment': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip alignment and downstream processes.',
    ),
    'call_variants': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Variant calling options',
        description='Specifies if variant calling will executed.',
    ),
    'variant_caller': NextflowParameter(
        type=typing.Optional[str],
        default='medaka',
        section_title=None,
        description="Specifies the variant caller that will used to call small variants (available are: medaka, deepvariant, and pepper_margin_deepvariant). Variant calling is only available if '--call_variants' is set and the protocol is set to `DNA`. Please note `deepvariant` and `pepper_margin_deepvariant` are only avaible if using singularity or docker.",
    ),
    'structural_variant_caller': NextflowParameter(
        type=typing.Optional[str],
        default='sniffles',
        section_title=None,
        description="Specifies the variant caller that will be used to call structural variants (available are: sniffles and cutesv). Structural variant calling is only available if '--call_variants' is set and the protocol is set to `DNA`.",
    ),
    'split_mnps': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Specifies if MNPs will be split into SNPs when using medaka.',
    ),
    'deepvariant_gpu': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Specifies whether to call variants with pepper_margin_deepvariant in GPU mode.',
    ),
    'phase_vcf': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Specifies if vcf will be phased when using medaka.',
    ),
    'skip_vc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip variant calling.',
    ),
    'skip_sv': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip structural variant calling.',
    ),
    'quantification_method': NextflowParameter(
        type=typing.Optional[str],
        default='bambu',
        section_title='Differential analysis options',
        description='Specifies the transcript quantification method to use (available are: bambu or stringtie2). Only available when protocol is cDNA or directRNA.',
    ),
    'skip_quantification': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip transcript quantification and differential analysis.',
    ),
    'skip_differential_analysis': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip differential analysis with DESeq2 and DEXSeq.',
    ),
    'jaffal_ref_dir': NextflowParameter(
        type=typing.Optional[str],
        default='for_jaffal',
        section_title='RNA fusion analysis options',
        description='Specifies the reference directory for JAFFAL.',
    ),
    'skip_fusion_analysis': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip differential analysis with DESeq2 and DEXSeq.',
    ),
    'skip_modification_analysis': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='RNA modification analysis options',
        description='Skip RNA modification analysis.',
    ),
    'skip_xpore': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip differential modification analysis with xpore.',
    ),
    'skip_m6anet': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip m6A detection with m6anet.',
    ),
    'skip_bigbed': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Process skipping options',
        description='Skip BigBed file generation.',
    ),
    'skip_bigwig': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip BigWig file generation.',
    ),
    'skip_nanoplot': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip NanoPlot.',
    ),
    'skip_fastqc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip FastQC.',
    ),
    'skip_multiqc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip MultiQC.',
    ),
    'skip_qc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip all QC steps apart from MultiQC.',
    ),
}

