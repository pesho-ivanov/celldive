"""
Helper tool that converts the output of tracer into the defined CellDive API JSON format.
See the current API documentation for details of the target JSON format

Features
- Tested with tracer v0.5.1 and 0.6.0
- can use the extended data in recombinants.txt (see our patch)
- reads clonal_group from cell_data.csv

Note: While many J(D)V identifiers reported by tracer follow the simple pattern TRAV6_TGTGCCGTACGGGATA_TRAJ12,
some do not, such as TRAV38-2_DV8_GGAGCCTCGAAAC_TRAJ58. TRAV38-2_DV8 actually refers to TRAV38-2DV8.
Thus, the V(D)J splitting algorithm here has been extended for such cases.

v0.7 2018-04-05 pisch
"""

import json
import pandas as pd
from pathlib import Path
import sys


version = 'v0.7'
sub_dir = 'filtered_TCRAB_summary'
recombinants_file = 'recombinants.txt'
clonality_file = 'cell_data.csv'
output_file = 'cd_tracer.json'


def error_exit(text):
    print('error:', text)
    exit(1)


def parse_arguments():
    if len(sys.argv) < 3:
        return None, 'missing arguments'

    context = dict()
    context['dir'] = str(sys.argv[1])
    context['sample_id'] = str(sys.argv[2])
    context['prefix'] = context['sample_id'] + '_'
    return context, None


def process_clonality_row(context, row):
    id = context['prefix'] + str(row['cell_name'])
    clonal_group = -1
    if not pd.isna(row['clonal_group']):
        clonal_group = int(row['clonal_group'])
    return id, clonal_group


def is_nucleotide_sequence(identifier):
    """True if at least one nucleotide, and no other characters than nucleotides or ()"""
    ok = False
    for c in identifier.lower():
        if c in 'actug':
            ok = True
        elif c not in '()':
            return False
    return ok


def split_vdj_identifier(identifier):
    # only beta chain has a d
    vdj = identifier.split('_')
    if len(vdj) == 3:
        return vdj[0], vdj[1], vdj[2]

    if len(vdj) < 3:
        print('incomplete V(D)J identifier', identifier)
        return 'error', 'error', 'error'

    print('split_vdj_identifier: handle special case: ', identifier)
    index = 1
    while index < len(vdj) and not is_nucleotide_sequence(vdj[index]):
        index += 1
    if index >= (len(vdj) - 1):
        print('error')
        return 'error', 'error', 'error'

    v = ''.join(vdj[:index])
    d_or_specifier = vdj[index]
    j = ''.join(vdj[index+1:])
    print('  ->', v, d_or_specifier, j)
    return v, d_or_specifier, j


def process_recombinants_row(context, row):
    data = dict()
    id = context['prefix'] + str(row['cell_name'])
    if pd.isna(row['recombinant_id']):
        return id, '', data

    locus = row['locus']
    #data['locus'] = locus
    recombinant_id = row['recombinant_id']
    data['recombinant_id'] = recombinant_id
    v, d_or_specifier, j = split_vdj_identifier(recombinant_id)
    data['v'] = v
    data['specifier'] = d_or_specifier      # for alpha and beta
    if locus == 'B':
        data['d'] = d_or_specifier          # beta chains only
    data['j'] = j
    data['productive'] = str(row['productive']).lower() == 'true'
    data['reconstructed_length'] = int(row['reconstructed_length'])

    if not context['extended']:
        return id, locus, data

    # extended data available
    # TPM     V_percentident  V_evalue        V_bitscore      J_percentident  J_evalue        J_bitscore
    data['TPM'] = float(row['TPM'])
    data['V_percentident'] = float(row['V_percentident'])
    data['V_evalue'] = float(row['V_evalue'])
    data['V_bitscore'] = float(row['V_bitscore'])
    data['J_percentident'] = float(row['J_percentident'])
    data['J_evalue'] = float(row['J_evalue'])
    data['J_bitscore'] = float(row['J_bitscore'])
    return id, locus, data


def run(context):
    path = context['dir'] + '/' + sub_dir + '/'
    recombinants_filename = path + recombinants_file
    clonality_filename = path + clonality_file
    json_filename = path + output_file
    if not Path(clonality_filename).is_file():
        error_exit('file ' + clonality_filename + ' does not exist.')
    if not Path(recombinants_filename).is_file():
        error_exit('file ' + recombinants_filename + ' does not exist.')

    clonality_data = pd.read_csv(clonality_filename, dtype={'cell': str}, sep=',')
    clonal_groups = dict()
    clonal_groups_count = dict()
    max_clonal_groups_count = 0
    max_clonal_groups_count_id = 0
    for (idx, row) in clonality_data.iterrows():
        id, clonal_group = process_clonality_row(context, row)
        clonal_groups[id] = clonal_group
        count = clonal_groups_count.get(clonal_group, 0)
        count += 1
        if count > max_clonal_groups_count:
            max_clonal_groups_count = count
            max_clonal_groups_count_id = clonal_group
        clonal_groups_count[clonal_group] = count

    recombinants_data = pd.read_csv(recombinants_filename, dtype={'cell': str}, sep='\t')
    if 'TPM' in recombinants_data.columns:
        print('extended data in', recombinants_file,'detected')
        context['extended'] = True
    else:
        context['extended'] = False
    cells = dict()
    current_id = ''
    current = None
    max_clonal_group = -1
    for (idx, row) in recombinants_data.iterrows():
        id, locus, data = process_recombinants_row(context, row)
        if current_id != id:
            current_id = id
            current = dict()
            current['cell_id'] = id
            clonal_group = clonal_groups.get(id, -1)
            if clonal_group > max_clonal_group:
                max_clonal_group = clonal_group
            current['clonal_group'] = clonal_group
            current['alphas'] = []
            current['betas'] = []
            cells[current_id] = current
        if locus == 'A':
            current['alphas'].append(data)
        elif locus == 'B':
            current['betas'].append(data)

    found_clonal_groups = max_clonal_group + 1
    no_chains_count = 0
    non_clonal_count = 0
    for data in cells.values():
        if data['clonal_group'] == -1:
            if len(data['alphas']) == 0 and len(data['betas']) == 0:
                data['clonal_group'] = None
                no_chains_count += 1
            else:
                non_clonal_count += 1
                #max_clonal_group += 1
                #data['clonal_group'] = max_clonal_group

    sample = dict()
    sample['sample_id'] = context['sample_id']
    sample['cells'] = cells
    samples = dict()
    samples[context['sample_id']] = sample
    to_json = dict()
    to_json['samples'] = samples
    with open(json_filename, 'w') as f:
        json.dump(to_json, f, indent=4)     # including pretty print at the moment
    print()
    print('input ', recombinants_filename)
    print('input ', clonality_filename)
    print('output', json_filename)
    print()
    sum = 0
    print('tracer found', found_clonal_groups, 'clonal groups')
    for i in range(found_clonal_groups):
        print('group', i, 'size', clonal_groups_count[i])
        sum += clonal_groups_count[i]
    print()
    print(non_clonal_count, 'individual cells with alpha and/or beta chains (clonal group == -1)')
    print(no_chains_count, 'cells/wells without alpha nor beta chains (clonal group == null/None)')
    sum += non_clonal_count + no_chains_count
    print('total', sum, 'cells/wells analyzed by tracer')
    print()


def main():
    print('tracer2json', version)
    context, err = parse_arguments()
    if err is None:
        run(context)
    else:
        print('usage:', sys.argv[0], '<input/output dir> <sample_id>')
        error_exit(err)


if __name__ == "__main__":
    main()
