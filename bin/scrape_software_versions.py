#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/nanoseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'guppy': ['v_guppy.txt', r"Version (\S+)"],
    'qcat': ['v_qcat.txt',  r"qcat (\S+)"],
    'pycoQC': ['v_pycoqc.txt', r"pycoQC v(\S+)"],
    'NanoPlot': ['v_nanoplot.txt', r"NanoPlot (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'GraphMap2': ['v_graphmap2.txt', r"Version: v(\S+)"],
    'minimap2': ['v_minimap2.txt', r"(\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'StringTie': ['v_stringtie.txt', r"(\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'R': ['v_r.txt', r"R version (\S+)"],
    'DESeq2': ['v_deseq2.txt', r"(\S+)"],
    'DRIMSeq': ['v_drimseq.txt', r"(\S+)"],
    'DEXSeq': ['v_dexseq.txt', r"(\S+)"],
    'stageR': ['v_stager.txt', r"(\S+)"],
    'BSgenome': ['v_bsgenome.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"]
}

results = OrderedDict()
results['nf-core/nanoseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['guppy'] = '<span style="color:#999999;\">N/A</span>'
results['qcat'] = '<span style="color:#999999;\">N/A</span>'
results['pycoQC'] = '<span style="color:#999999;\">N/A</span>'
results['NanoPlot'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['GraphMap2'] = '<span style="color:#999999;\">N/A</span>'
results['minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['StringTie'] = '<span style="color:#999999;\">N/A</span>'
results['featureCounts'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['DESeq2'] = '<span style="color:#999999;\">N/A</span>'
results['DRIMSeq'] = '<span style="color:#999999;\">N/A</span>'
results['DEXSeq'] = '<span style="color:#999999;\">N/A</span>'
results['stageR'] = '<span style="color:#999999;\">N/A</span>'
results['BSgenome'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/nanoseq Software Versions'
section_href: 'https://github.com/nf-core/nanoseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
