#!/usr/bin/env python

import pandas as pd
import argparse
import csv

argparser = argparse.ArgumentParser()
argparser.add_argument('--samplesheet', type=str)
ARGS = argparser.parse_args()
samplesheet = ARGS.samplesheet

# get idx of Data tag
data_tag_search = '[Data]'
data_index = 0
with open(samplesheet, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for idx, row in enumerate(reader):
        if data_tag_search in row:
            data_index = idx

sample_df = pd.read_csv(samplesheet, skiprows=range(0, data_index + 1))

result = sample_df[sample_df["Fastq_path"].isnull()]

basecalling_needed = "true"
if result.empty:
	basecalling_needed = "false"

result_file = open (basecalling_needed + ".txt", "w")
result_file.close()