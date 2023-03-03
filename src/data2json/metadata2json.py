#!/usr/bin/env python3

"""
Helper tool that converts meta data from the provenance table (tsv format) to json
See the current API documentation for details of the target JSON format

Note: The tool assumes that this .tsv file is in the base folder of the samples.
As an example, assuming provided filename <path>/Data_provenance_helper_for_tsv.tsv
and assuming a folder structure (typical for this project here) of
<path>/<subject_nr>/<type>/ for the actual data, and "<subject_nr>_<type>" being the sample_id,
this tool scans these subfolders for meta.txt files.

E.g. info in <path>/03/SE/meta.txt will be added for the cells of sample 03_SE.

These meta.txt files currently provide microscope information of the cells for each sample.

This tool is not generic but heavily adjusted to our current project structure.
Please adjust as needed.

v0.4 2018-04-06 pisch
"""

import json
import os
import pandas as pd
from pathlib import Path
import sys


version = 'v0.4'
output_file = 'cd_metadata.json'
meta_file = 'meta.csv'

def error_exit(text):
    print('error:', text)
    exit(1)


def parse_arguments():
    if len(sys.argv) < 2:
        return None, 'missing arguments'

    context = dict()
    context['file'] = os.path.abspath(sys.argv[1])  # abs path needed as dirname splits this name
    context['path'] = os.path.dirname(context['file'])
    return context, None


def process_provenance_row(context, row):
    data = dict()
    internal = str(row['internal_sample_name'])
    alias = str(row['official_alias'])
    id = str(row['sample_id'])
    color = str(row['color'])
    directory = str(row['directory'])
    clonal_alpha = str(row['clonal_alpha'])
    clonal_beta = str(row['clonal_beta'])
    data['sample_id'] = id
    data['alias'] = alias
    data['internal_name'] = internal
    data['color'] = color
    data['directory'] = directory
    data['clonal_alpha'] = clonal_alpha
    data['clonal_beta'] = clonal_beta
    return id, data


def process_cells_meta_data(context, sample, prefix, meta_name):
    print('processing', meta_name, 'for', sample['sample_id'])
    meta = pd.read_csv(meta_name, dtype={'well': str}, sep=r'\s+')  # accept any whitespace as sep
    cells = dict()
    sample['cells'] = cells
    for (idx, row) in meta.iterrows():
        id = str(row['well'])
        cell_id = prefix + id
        count = str(row['count'])
        cell = dict()
        cell['cell_id'] = cell_id
        cell['microscope'] = count
        cells[cell_id] = cell
    return None


def scan_samples(context, samples):
    level1 = os.listdir(context['path'])
    for l1 in level1:
        path = context['path'] + '/' + l1
        if not Path(path).is_dir():
            continue
        level2 = os.listdir(path)
        for l2 in level2:
            path2 = path + '/' + l2
            if not Path(path2).is_dir():
                continue
            meta_name = path2 + '/' + meta_file
            if not Path(meta_name).is_file():
                continue

            # handle cell meta information of a sample
            sample_id = l1 + '_' + l2
            sample = samples.get(sample_id, None)
            if sample is None:
                error_exit('no basic metadata (provenance) available for sample' + sample_id)
            err = process_cells_meta_data(context, sample, sample_id + '_', meta_name)
            if err is not None:
                error_exit(err)


def run(context):
    json_filename = context['path'] + '/' + output_file
    if not Path(context['file']).is_file():
        error_exit('file ' + context['file'] + ' does not exist.')
    input = pd.read_csv(context['file'], dtype={'cell': str}, sep='\t')
    samples = dict()
    for (idx, row) in input.iterrows():
        id, data = process_provenance_row(context, row)
        samples[id] = data
    scan_samples(context, samples)
    to_json = dict()
    to_json['samples'] = samples
    with open(json_filename, 'w') as f:
        json.dump(to_json, f, indent=4)     # including pretty print at the moment
    print('input ', context['file'])
    print('output', json_filename)


def main():
    print('metadata2json', version)
    context, err = parse_arguments()
    if err is None:
        run(context)
    else:
        print('usage:', sys.argv[0], '<tsv file with meta data>')
        error_exit(err)


if __name__ == "__main__":
    main()
