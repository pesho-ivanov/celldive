from collections import defaultdict
from itertools import combinations
import re
import json

from utils import *

class TCRs:
    NONPRODUCTIVE_TAG = "NONPRODUCTIVE"

    def __init__(self, _cell2tcrs, sample_id):
        def remove_empty(cell2tcrs):
            return { cell: tcrs for cell, tcrs in cell2tcrs.items() if tcrs['alphas'] or tcrs['betas'] }

        self._cell2tcrs = remove_empty(_cell2tcrs)   # cell -> dict['alphas', 'betas'] -> tcr list
        self.sample_id = sample_id
        #self.generate_VDJ_columns()
        
        #self.print_stats()

        
    def __repr__(self):
        return '{} TCR chains reconstructed'.format(len(self.cells()))
     
    #def generate_VDJ_columns(self):
    #    def extract_VDJ(tcr, regex):
    #        match = re.search(regex, tcr)
    #        if match:
    #            return match.groups()[0]
    #        else:
    #            return None # np.nan

    #    seqs = list(self.tcrs.index.get_level_values('seq')) #.unique())
    #    self.tcrs['V']  = [ extract_VDJ(seq, r'^(TR[AB]V.*?)_')   for seq in seqs ]
    #    self.tcrs['D']  = [ extract_VDJ(seq, r'_(D.*?)_')         for seq in seqs ]
    #    self.tcrs['junc'] = [ extract_VDJ(seq, r'_([ACGT\(\)]*?)_') for seq in seqs ]
    #    self.tcrs['J']  = [ extract_VDJ(seq, r'_(TR[AB]J.*?)$')   for seq in seqs ]
    #    
    #    #for ind, cols in self.tcrs.iterrows():
    #    #    whole = ind[2]
    #    #    all_parts = filter(lambda s: s, [ cols['V'], cols['D'], cols['junc'], cols['J'] ])
    #    #    merged = '_'.join(all_parts)
    #    #    if whole != merged:
    #    #        #print(red("WARNING: ") + "%s incorrectly split into %s" % (blue(whole), blue(merged)))
    #    #        None
    #           
    #def get_all_VDJ(self, which: "'V', 'D', 'junc' or 'J'"):
    #    return self.tcrs[which].dropna().unique()

    @classmethod
    def from_file(cls, true_index, input_file, *, sample_id, filt_tcr_chains): 
        if not check_file(input_file, "TCRs"):
            raise Exception('No TCR file for {}'.format(sample_id))

        def mask_nonproductive_inplace(cell2tcrs):
            #print('Masking nonproductive tcrs.')
            new_cell2tcrs = dict()
            for cell, tcrs in cell2tcrs.items():
                for ch in ['alphas', 'betas']:
                    for chain in tcrs[ch]:
                        if not 'productive' in chain:
                            raise Exception('No productivenes information for cell {}'.format(cell))
                        #print(chain)
                        if not chain['productive']:
                            chain['recombinant_id'] = '{}_{}_{}_{}'.format(
                                TCRs.NONPRODUCTIVE_TAG, 'A' if ch=='alpha' else 'B', chain['recombinant_id'], cell)
                            #tcrs[ch].pop(chain['recombinant_id'])

        data = json.load(open(input_file))
        #if 'samples' not in data:
        #    Warning('The TCR data does not have ``samples\'\' as a top level.')
        #    return
        #if sample_id not in data['samples']:
        #    Warning('The TCR data does not have ``{}\'\' inside ``samples\'\'.'.format(sample_id))
        #    return

        _cell2tcrs = data['samples'][sample_id]['cells']

        if filt_tcr_chains and 'nonproductive' in filt_tcr_chains and filt_tcr_chains['nonproductive']:
            mask_nonproductive_inplace(_cell2tcrs)

        #from pprint import pprint
        #pprint(_cell2tcrs)

        #print('tcrs.py: cell2tcrs: ', _cell2tcrs)
        return cls(_cell2tcrs, sample_id)
        
    #@classmethod
    #def from_file(cls, true_index, input_file, *, sample_id): 
    #    if not check_file(input_file, "TCRs"):
    #        return None

    #    # Clonotype results from TraCeR
    #    rawTCRs = pd.read_csv(input_file, dtype={'cell': str}, sep='\t')
    #    rawTCRs = rawTCRs.rename(columns={'cell_name': 'cell'})
    #    #rawTCRs = rawTCRs.set_index('cell')
    #    rawTCRs = rawTCRs[ pd.notnull(rawTCRs['recombinant_id']) ]
    #    # TODO: if there is a recombinant is not
    #    #cls.print_stats(rawTCRs)
    #    
    #    columns_rename = {
    #        'cell': 'cell',
    #        'locus': 'chain',
    #        'recombinant_id': 'seq',
    #        'reconstructed_length': 'len',
    #    }
    #    if (set(rawTCRs.columns) <= set(columns_rename)):
    #        raise ValueError("Wrong columns: expected %s, received %s" %
    #                         ( blue(iter2str(columns_rename.keys())), blue(iter2str(rawTCRs.columns)) ) )
    #    rawTCRs.rename(columns=columns_rename, inplace=True)
    #    
    #    tcrs = rawTCRs.set_index(['cell', 'chain', 'seq'], # not needed?
    #                             verify_integrity=True)
    #    tcrs.index = pd.MultiIndex.from_tuples([(str(x[0]), x[1], x[2]) for x in tcrs.index],
    #                                           names=['cell', 'chain', 'seq'])
    #    tcrs.sort_index(inplace=True)
    #    #tcrs.index.levels[0] = tcrs.index.levels[0].astype(str, copy=False)

    #    add_prefix_tag_multiindex(tcrs, sample_id) 
    #    
    #    #tcrs = to_nice_TCRs(rawTCRs)
    #    assert_same_indexes(true_index, tcrs.index.levels[0], "TCRs")
    #    #true_df = pd.DataFrame(index=true_index)
    #    # Joining a single Index to a Multi-index:
    #    # http://pandas.pydata.org/pandas-docs/stable/merging.html
    #    #tcrs = pd.concat([true_df, tcrs], axis=1, verify_integrity=True)

    #    # sieve
    #    tcrs = tcrs[ tcrs['productive'] == True ]
    #    tcrs = tcrs[ tcrs['len'] > 300 ]
    #    #display(tcrs)

    #    return cls(tcrs)
    
    def select_cells(self, tcr_seq):
        #pd.DataFrame(index=index.get_level_values('cell'))
        #subset = self.tcrs.loc[(slice(None), slice(None), [tcrs]), :]

        subset = [ cell for cell in self.cells() if tcr_seq in self.cell2seqs(cell) ]
                   # self._get_recombinants(self._cell2tcrs[c]['alphas']) | \
                   #               self._get_recombinants(self._cell2tcrs[c]['betas']) )) ]
        #display(subset)
        #return subset.index.get_level_values('cell').tolist()
        return subset
    
    #def select_nonclonal_cells(self):
    #    d = tcr2cell_dict(f):

    def subtcrs(self, *, cells):  #, chain: "A or B"=None):
        #print('tcrs.py::subtcrs cells', cells)
        # TODO
        #intersect = set(cells).intersection(self.tcrs.index.tolist())
        #assert_same_indexes(cells, intersect)
        if cells:
            intersect = set(cells).intersection(self.cells())
            #tcrs = self.tcrs.loc[list(intersect)]
            tcrs = { cell: self._cell2tcrs[cell] for cell in intersect }
        else:
            tcrs = self._cell2tcrs
        #if chain:
        #    tcrs = tcrs.xs(chain, level='chain')

        #samples = { 'samples': { self.sample_id: tcrs } }
        return TCRs(tcrs, sample_id=self.sample_id)
    
    def concat(self, other, sample_id):
        #return TCRs(pd.concat([self.tcrs, other.tcrs]))
        return TCRs({**self._cell2tcrs, **other._cell2tcrs}, sample_id)
    
    def cells(self):
        #cells_level = self.tcrs.index.get_level_values('cell').unique()
        return sorted(self._cell2tcrs.keys())

    def cell2struct(self, cell):
        return self._cell2tcrs[cell]
    
    def cell2seqs(self, cell, *, chains=['alpha', 'beta']):
        seqs = [ tcr['recombinant_id'] for chain in chains for tcr in self.cell2tcrs(cell, chain) ]
        seqs = [ seq for seq in seqs if not seq.startswith(TCRs.NONPRODUCTIVE_TAG)]
        return seqs

        #return self.cell2tcrs(cell, 'alphas') | self.cell2tcrs(cell, 'betas')
        #return self._get_recombinants(self._cell2tcrs[cell]['alphas']) | \
        #       self._get_recombinants(self._cell2tcrs[cell]['betas'])
        #cell_df = self.tcrs.xs(cell, level='cell')
        #return cell_df.index.get_level_values('seq').tolist()
    
    #def get_tcrs(self):
    #    seqs = self.tcrs.index.get_level_values('seq')
    #    return seqs.unique()
    
    #@staticmethod
    #def _chain2tcrs(df, chain: "A or B"):
    #    df2 = df.xs(chain, level='chain')
    #    seqs = df2.index.get_level_values('seq')
    #    return seqs.unique().tolist()
    
    def chain2tcrs(self, chain: "alpha or beta"):
        #return TCRs._chain2tcrs(self.tcrs, chain)
        tcrs_intersected = set()
        for cell in self._cell2tcrs.keys():
            tcrs_intersected |= set(self.cell2seqs(cell, chains=[chain])) #tcrs[chain + 's']
        return tcrs_intersected
    
    def cell2tcrs(self, cell, chain):  # both alpha and beta
        #print('tcrs.py: ', cell, self._cell2tcrs[cell])
        if chain:
            return self._cell2tcrs[cell][chain + 's']
        else:
            return self._cell2tcrs[cell]['alphas'] + self._cell2tcrs[cell]['betas']
        #return { chain for chain in chains } ['recombinant_id']
        #A = self._cell2tcrs[cell]['alphas']
        #B = self._cell2tcrs[cell]['betas']
        #df = self.tcrs.loc[cell]
        #A = TCRs._chain2tcrs(df, 'A')
        #B = TCRs._chain2tcrs(df, 'B')
        #return { 'alphas': A, 'betas': B }

    def print_stats(self):
        #print('tcr stats...')
        #display(self.tcrs)
        #print ("  -- %d cells with at least one reconstructed chain" % len(self.cells()))
        
        #unique_recs = self.tcrs['recombinant_id'].unique()
        #other_recs = [rec for rec in unique_recs
        #              if not rec.startswith('TRA') and not rec.startswith('TRB')]
        #print ("%d different recombinants of which:" % len(unique_recs))
        #print ("  -- %d alpha" % len([rec for rec in unique_recs if rec.startswith('TRA')]))
        #print ("  -- %d beta" % len([rec for rec in unique_recs if rec.startswith('TRB')]))
        #print ("  -- %d other: " % len(other_recs), other_recs)
        return
    
    def tcr2cells_dict(self):
        d = defaultdict(list)
        for cell in self.cells():
            for tcr in self.cell2seqs(cell):
                d[tcr].append(cell)
        return d

    def getEdges(self):
        """Returns a list of edges (from_cell, to_cell, TCR_seq)"""
        d = self.tcr2cells_dict()
        edges = []
        for tcr, cells in d.items():
            for a, b in combinations(cells, r=2):
                edges.append( (a, b, tcr) )
        return edges
        
    #def to_nice_TCRs(raw):
    #    used_columns = ['cell', 'recombinant_id', 'locus', 'productive']
    #    if not set(used_columns) <= set(raw.columns):
    #        raise ValueError("Wrong columns! should include %s but are only %s"
    #                         % (', '.join(used_columns), ', '.join(raw.columns)))
    #
    #    #display(raw)
    #    recombList = [ (raw['cell'][i],
    #                    raw['recombinant_id'][i],
    #                    raw['locus'][i],
    #                    raw['productive'][i]) for i in raw.index ]
    #    recombDict = {}
    #    meta = pd.DataFrame(columns=used_columns)
    #    #for cell, recomb, locus, prod in recombList:
    #    #    recombDict.setdefault(str(cell), []).append({'recomb_id': recomb,
    #    #                                                 'locus': locus,
    #    #                                                 'prod': prod})
    #    #for cell, recomb in recombDict.items():
    #    #    meta.set_value(cell, 'raw', recomb)
    #    
    #    return raw
