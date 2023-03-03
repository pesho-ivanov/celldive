#!/usr/bin/env python3

"""
Helper tool that recursively goes through the 2-layer folder structure of usz samples
and combines meta data from microscopic.csv and sequencing_reads.tsv into one
file meta.csv that will be read by CellDive.

Assuming the start path <path>, it assumes a folder structure
<path>/<subject_nr>/<type>/ for the actual data, and "<subject_nr>_<type>" being the sample_id.

E.g. info in <path>/03/SE/meta.txt will be added for the cells of sample 03_SE.

This tool is not generic but heavily adjusted to our current project structure.
Please adjust as needed.

Note: Since cell name needs to be extracted from file names in sequencing_reads.tsv
and the file naming scheme changed multiple times during the project, this tool
offers an interactive mode to define how the cell name can be extracted from file names.
Additionally, this setting is stored and can be used during re-runs without need
for additional interaction.

v0.6 2018-04-17 pisch
"""

import json
import os
import pandas as pd
from pathlib import Path
import sys


version = 'v0.6'
output_file = 'meta.csv'
microscopic_file = 'microscopic.csv'
sequencing_file = 'sequencing_reads.tsv'
config_file = 'sample_metadata2csv_config.json'


def error_exit(text):
    print('error:', text)
    exit(1)


def parse_arguments():
    if len(sys.argv) < 2:
        return None, 'missing arguments'

    context = dict()
    context['path'] = os.path.abspath(sys.argv[1])
    return context, None


def get_well_name_bounds(parts):
    print()
    print('select boundaries to define the cell name')
    index = 0
    association = dict()
    special_case = 'plate1'
    for p in parts:
        if p.startswith(special_case) and len(p) > len(special_case):
            p = p[len(special_case):]
        association[index] = p
        print(index, '::', p)
        index += 1
    first = int(input('first? > '))
    second = int(input('second? > '))
    print('producing cell name', get_well_name(parts, first, second))
    return first, second


def get_well_name(parts, first, second):
    special_case = 'plate1'
    for i in range(len(parts)):
        p = parts[i]
        if p.startswith(special_case) and len(p) > len(special_case):
            parts[i] = p[len(special_case):]
    return '_'.join(parts[first:second + 1])


def append_sequencing_reads_data(context, filename, sample_id, data):
    config = context['config']
    index_name = sample_id + '_' + sequencing_file
    file_config = config.get(index_name, None)
    if file_config is None:
        file_config = dict()
        config[index_name] = file_config

    file_tag = 'Read1 [File]'
    reads_tag = 'Read Count'
    column_tag = 'read_count'  # in the data output
    well_tag = 'well'
    print('processing', filename)
    sequencing_data = pd.read_csv(filename, dtype={'well': str}, sep='\t')  # accept tabs only
    if file_tag not in sequencing_data.columns:
        error_exit("Column '" + file_tag + "' not found in" + filename)
    if reads_tag not in sequencing_data.columns:
        error_exit("Column '" + reads_tag + "' not found in" + filename)

    read_counts = dict()
    first_part = file_config.get('first', None)
    last_part = file_config.get('last', None)
    for (idx, row) in sequencing_data.iterrows():
        fstr = row[file_tag]
        if not isinstance(fstr, str):
            print('missing file name in row', row)
            continue
        parts = row[file_tag].split('_')
        if first_part is None:
            first_part, last_part = get_well_name_bounds(parts)
            file_config['first'] = first_part
            file_config['last'] = last_part
        read_counts[get_well_name(parts, first_part, last_part)] = int(row[reads_tag])

    if well_tag not in data.columns:
        error_exit("Column '" + well_tag + "' not found in" + microscopic_file)
    data[column_tag] = 0
    for (idx, row) in data.iterrows():
        well = row[well_tag]
        count = read_counts.get(well, None)
        if count is None:
            continue
        data.loc[idx, column_tag] = int(count)
    return data


def read_cells_microscopic_data(context, filename):
    print('processing', filename)
    microscopic_data = pd.read_csv(filename, dtype={'well': str}, sep=r'\s+')  # accept any whitespace as sep
    return microscopic_data


def scan_samples(context):
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
            microscopic_name = path2 + '/' + microscopic_file
            if not Path(microscopic_name).is_file():
                continue

            # handle cell microscopic meta information of a sample
            sample_id = l1 + '_' + l2
            data = read_cells_microscopic_data(context, microscopic_name)

            # append sequencing reads information
            sequencing_name = path2 + '/' + sequencing_file
            if Path(sequencing_name).is_file():
                data = append_sequencing_reads_data(context, sequencing_name, sample_id, data)
            else:
                print('No sequencing meta data found in', path2)

            # store data
            output_name = path2 + '/' + output_file
            data.to_csv(output_name, sep='\t', index=False)


def run(context):
    if not Path(context['path']).is_dir():
        error_exit('path ' + context['path'] + ' is not correct.')

    config_filename = context['path'] + '/' + config_file
    context['config'] = dict()
    if Path(config_filename).is_file():
        with open(config_filename, 'r') as f:
            context['config'] = json.load(f)
        print('Config file read.')
    else:
        print('Config file not found; new file will be generated.')

    scan_samples(context)

    with open(config_filename, 'w') as f:
        json.dump(context['config'], f, indent=4)
    print('Config file written.')


def main():
    print('sample_metadata2csv', version)
    context, err = parse_arguments()
    if err is None:
        run(context)
    else:
        print('usage:', sys.argv[0], '<path>')
        error_exit(err)


if __name__ == "__main__":
    main()
