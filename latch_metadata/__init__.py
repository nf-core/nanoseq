from latch.types.directory import LatchDir
from latch.types.metadata import LatchAuthor, NextflowMetadata, NextflowRuntimeResources

from .parameters import flow, generated_parameters

NextflowMetadata(
    display_name="nf-core/nanoseq",
    author=LatchAuthor(
        name="nf-core",
    ),
    repository="https://github.com/latchbio-nfcore/nanoseq",
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=100,
    ),
    flow=flow,
    log_dir=LatchDir("latch:///your_log_dir"),
)
