#!/usr/bin/env python

import sys

outfile_name = sys.argv[1]
samples = sys.argv[2:]
outfile = open(outfile_name, "w")
outfile.write(
    "notes: Pairwise comparison without replicates with default parameter setting." + "\n" + "\n" + "data:" + "\n"
)

dict = {}
for sample in samples:
    sample = sample.strip("[").strip("]").split(",")
    path_to_sample, condition = sample[0], sample[1]
    if condition not in dict:
        dict[condition] = [path_to_sample]
    else:
        dict[condition].append(path_to_sample)

for condition in dict:
    outfile.write("    " + condition + ":" + "\n")
    r = 0
    for sample in dict[condition]:
        r += 1
        outfile.write("        rep" + str(r) + ": " + sample + "\n")

outfile.write("\n" + "out: ./diffmod_outputs" + "\n" + "\n")

outfile.write(
    "method:"
    + "\n"
    + "    prefiltering:"
    + "\n"
    + "        method: t-test"
    + "\n"
    + "        threshold: 0.1"
    + "\n"
    + "\n"
)

outfile.close()
