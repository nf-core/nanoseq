#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/nanoseq': ['pipeline.version', r"(\S+)"],
    'Nextflow': ['nextflow.version', r"(\S+)"],
    'guppy': ['guppy.version', r"Version (\S+)"],
    'pycoQC': ['pycoqc.version', r"pycoQC v(\S+)"],
    'NanoPlot': ['nanoplot.version', r"NanoPlot (\S+)"],
    'FastQC': ['fastqc.version', r"FastQC v(\S+)"],
    'GraphMap': ['graphmap.version', r"Version: v(\S+)"],
    'minimap2': ['minimap2.version', r"(\S+)"],
    'Samtools': ['samtools.version', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'rmarkdown': ['rmarkdown.version', r"(\S+)"],
    'MultiQC': ['multiqc.version', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['nf-core/nanoseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['guppy'] = '<span style="color:#999999;\">N/A</span>'
results['pycoQC'] = '<span style="color:#999999;\">N/A</span>'
results['NanoPlot'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['GraphMap'] = '<span style="color:#999999;\">N/A</span>'
results['minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['rmarkdown'] = '<span style="color:#999999;\">N/A</span>'
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
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/nanoseq Software Versions'
section_href: 'https://github.com/nf-core/nanoseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
