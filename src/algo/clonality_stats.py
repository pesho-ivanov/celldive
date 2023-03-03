from collections import OrderedDict
from scipy.stats import fisher_exact
import pandas as pd
import json

from utils import *

class ClonalityStats:
    def __init__(self, directory):
        self.directory = directory
        self.stats = self.extract_all_stats(directory)

    @staticmethod
    def extract_all_stats(directory):
        '''Returns a dataframe with samples as rows and cell types as columns.'''
        directory = Path(directory)
        directory.exists()
        files = list(directory.glob('*.json'))
        
        all_stats = pd.DataFrame()
        for fn in files:
            tag, sample_stats = ClonalityStats.read_sample(fn)
            #print('ordereddict: ', sample_stats)
            s = pd.Series(sample_stats, name=tag)
            all_stats = all_stats.append(s)[s.index]  ## ??
            #all_stats.loc[tag] = pd.Series(sample_stats)
        if all_stats is None:
            warning('Empty stats from {}'.format(self.directory))
        #print(all_stats)
        sums = pd.Series(all_stats.sum(), name='sum')
        sums[ [ type(a) is str for a in sums ] ] = ''
        all_stats = all_stats.append(sums)
        return all_stats.sort_index()

    @staticmethod
    def get_tcr_stats_cells(parts): 
        cells_with_alphas, cells_with_betas = set(), set()
        cells_with_2alphas, cells_with_2betas = set(), set()

        for part, cells in parts.items():
            if part != 'clonotype':
                for cell_name, cell in cells.items():
                    alphas = cell['alphas'] if 'alphas' in cell else []
                    betas = cell['betas'] if 'betas' in cell else []
                    seqs = alphas + betas
                    if len(alphas) >= 1:
                        cells_with_alphas.add(cell_name)
                    if len(betas) >= 1:
                        cells_with_betas.add(cell_name)
                    if len(alphas) >= 2:
                        cells_with_2alphas.add(cell_name)
                    if len(betas) >= 2:
                        cells_with_2betas.add(cell_name)

        res = OrderedDict()
        res['#cells with >=1 alphas'] = len(cells_with_alphas)
        res['#cells with >=1 betas'] = len(cells_with_betas)
        res['#cells with >=1 alphas and >=1 betas'] = len(cells_with_alphas & cells_with_betas)
        res['#cells with >=2 alphas'] = len(cells_with_2alphas)
        res['#cells with >=2 betas'] = len(cells_with_2betas)
        res['#cells with >=2 alphas and >=2 betas'] = len(cells_with_2alphas & cells_with_2betas)
        return res

    @staticmethod
    def get_problematic_cells(parts): 
        cell_without_tcrs, tcrs_without_cell = [], []

        for part, cells in parts.items():
            if part != 'clonotype':
                for cell_name, cell in cells.items():
                    count = cell['microscope']
                    alphas = cell['alphas'] if 'alphas' in cell else []
                    betas = cell['betas'] if 'betas' in cell else []
                    seqs = alphas + betas
                    if count > 0 and len(seqs) == 0:
                        cell_without_tcrs.append(cell_name)
                    elif count == 0 and len(seqs) > 0:
                        tcrs_without_cell.append(cell_name)

        res = OrderedDict()
        res['# cells w/o tcrs'] = len(cell_without_tcrs)
        res['# tcr w/o cells'] = len(tcrs_without_cell)
        return res

        #return { #'cell w/o tcrs': ', '.join(cell_without_tcrs),
        #         '# cells w/o tcrs': len(cell_without_tcrs),
        #         #'tcr w/o cells': ', '.join(tcrs_without_cell),
        #         '# tcr w/o cells': len(tcrs_without_cell) }

    @staticmethod
    def read_sample(file_json):
        with open(file_json) as fh:
            whole_json = json.load(fh, object_pairs_hook=OrderedDict)
        if len(whole_json) > 1:
            raise Exception('Wrong json format. The top level should only contain one entry.')
        tag, parts = list(whole_json.items())[0]

        sample_stats = OrderedDict()
        for key, value in parts['clonotype'].items():
            sample_stats[key] = value

        microscopic = {0: 0, 1: 0, 2: 0, 3: 0}
        for part, cells in parts.items():
            if part != 'clonotype':
                sample_stats[part] = len(cells)
                for cell_name, cell in cells.items():
                    count = cell['microscope']
                    if count > 3:
                        count = 3
                    microscopic[count] += 1

        sample_stats['X+ (main and main+)'] = sample_stats['main'] + sample_stats['main+'] if 'main' in sample_stats and 'main+' in sample_stats else '-1'
        sample_stats['Y+ (small and single)'] = sample_stats['small'] + sample_stats['single'] if 'small' in sample_stats and 'single' in sample_stats else '-1'

        for key, value in microscopic.items():
            sample_stats['{} cells in well'.format(key)] = value
        ge1 = sum([value for key, value in microscopic.items() if key >= 1])
        ge2 = sum([value for key, value in microscopic.items() if key >= 2])
        sample_stats['>=1 cells in well (ge1)'] = ge1
        sample_stats['>=2 cells in well (ge2)'] = ge2
        sample_stats['doublets+ percentage (ge2/ge1)'] = '{:.2f}'.format(100.0*ge2/ge1)

        for key, value in ClonalityStats.get_tcr_stats_cells(parts).items():
            sample_stats[key] = value
               
        for key, value in ClonalityStats.get_problematic_cells(parts).items():
            sample_stats[key] = value

        return tag, sample_stats

    def dump(self, out_file):
        print('Writing to {}'.format(blue(out_file.name)))
        self.stats.to_csv(out_file)  # astype('int64').


# ---------------------------- OLD ------------------
    #def __init__(self, batch):
    #    self.batch = batch
            
    def cells_with_alpha_beta(self, value):
        cells_with_alpha = 0
        cells_with_beta = 0
        for cell in self.batch.tcrs.cells():
            tcrs_A = self.batch.tcrs.cell2seqs(cell, chains=['alpha'])
            tcrs_B = self.batch.tcrs.cell2seqs(cell, chains=['beta'])
            if len(tcrs_A) == value:
                cells_with_alpha += 1
            if len(tcrs_B) == value:
                cells_with_beta += 1 
        return cells_with_alpha, cells_with_beta

    def cells_with_alpha(self, cl, noncl):
        four_vals = []
        pairs = []
        for self.batch, name in [(cl, 'clonal'), (noncl, 'nonclonal')]:
            cells_with_alpha = 0
            cells_without_alpha = 0
            for cell in self.batch.tcrs.cells():
                tcrs_A = self.batch.tcrs.cell2seqs(cell, chains=['alpha'])
                tcrs_B = self.batch.tcrs.cell2seqs(cell, chains=['beta'])
                if len(tcrs_B) == 1:
                    if len(tcrs_A) == 1:
                        cells_with_alpha += 1
                    elif len(tcrs_A) == 0:
                        cells_without_alpha += 1 
            pairs.append(('{}_cells_with_alpha|beta'.format(name), cells_with_alpha))
            pairs.append(('{}_cells_without_alpha|beta'.format(name), cells_without_alpha))
            four_vals.append([cells_with_alpha, cells_without_alpha])
#
        oddsratio, pvalue = fisher_exact(four_vals)
        pairs.append(('fisher_pvalue_equal', pvalue))
        return pairs

    def cells_with_beta(self, cl, noncl):
        four_vals = []
        pairs = []
        for self.batch, name in [(cl, 'clonal'), (noncl, 'nonclonal')]:
            cells_with_beta = 0
            cells_without_beta = 0
            for cell in self.batch.tcrs.cells():
                tcrs_A = self.batch.tcrs.cell2seqs(cell, chains=['alpha'])
                tcrs_B = self.batch.tcrs.cell2seqs(cell, chains=['beta'])
                if len(tcrs_A) == 1:
                    if len(tcrs_B) == 1:
                        cells_with_beta += 1
                    elif len(tcrs_B) == 0:
                        cells_without_beta += 1 
            pairs.append(('{}_cells_with_beta|alpha'.format(name), cells_with_beta))
            pairs.append(('{}_cells_without_beta|alpha'.format(name), cells_without_beta))
            four_vals.append([cells_with_beta, cells_without_beta])

        oddsratio, pvalue = fisher_exact(four_vals)
        pairs.append(('fisher_pvalue_equal', pvalue))
        return pairs

    #def fischer(cl, noncl, is_beta):
    #    # Contingency table for Fischer's exact test
    #    #     cells    | clonal | noncl |     |
    #    # with 1 alpha |   a    |   b   | a+b |
    #    # with 0 alphas|   c    |   d   | c+d |
    #    #              |  a+c   |  b+d  |  n  |
    #    a = self.cells_with_alpha_beta(cl, 1)[is_beta]
    #    b = self.cells_with_alpha_beta(noncl, 1)[is_beta]
    #    c = self.cells_with_alpha_beta(cl, 0)[is_beta]
    #    d = self.cells_with_alpha_beta(noncl, 0)[is_beta]
    #    oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
    #    return oddsratio, pvalue
    
    @staticmethod
    def sets2counts(s):
        counts = pd.Series()
        for name, cells in s.items():
            counts[name] = len(cells)
        return counts

    def run(self, out_counts_fn, out_cells_fn):
        if self.batch.cell_part is None:
            Warning('No clonality information for batch {}. No clonality stats produced.'.format(self.batch.tag))
            return
        
        cl, noncl, unk, nonrec = self.batch.get_clonality_subbatches()

        # alpha / beta level
        #cells = batch.tcrs.cells() 
        #cells_with_alpha, cells_with_beta = batch.tcrs.cells_with_alpha_beta()
        #print(cells_with_alpha, cells_with_beta)

        array = [
            ('all', self.batch.cells()),
            ('clonal', cl.cells()),
            ('nonclonal', noncl.cells()),
            #('unknown', len(unk.cells())),
            #('nonreconstructed', len(nonrec.cells())),
            #('total_cells_with_alpha', self.cells_with_alpha_beta(batch, 1)[0]),
            #('total_cells_without_alpha', self.cells_with_alpha_beta(batch, 0)[0]),
       #     ('clonal_cells_with_alpha', self.cells_with_alpha_beta(cl, 1)[0]),
       #     ('clonal_cells_without_alpha', self.cells_with_alpha_beta(cl, 0)[0]),
       #     ('nonclonal_cells_with_alpha', self.cells_with_alpha_beta(noncl, 1)[0]),
       #     ('nonclonal_cells_without_alpha', self.cells_with_alpha_beta(noncl, 0)[0]),
       #     ('pvalue(difference in alpha)', fischer(cl, noncl, 0)[1]),
       #     #('unknown_cells_with_alpha', self.cells_with_alpha_beta(unk, 1)[0]),
       #     #('unknown_cells_with_beta', self.cells_with_alpha_beta(unk, 1)[1]),
       #     #('clonal_cells_with_alpha / nonclonal_cells_with_alpha', self.cells_with_alpha_beta(cl, 1)[0]/ self.cells_with_alpha_beta(noncl, 1)[0]),
       #     #('clonal_cells_with_beta / nonclonal_cells_with_beta', self.cells_with_alpha_beta(cl, 1)[1]/ self.cells_with_alpha_beta(noncl, 1)[1]),
       #     
       #     #('total_cells_with_beta', self.cells_with_alpha_beta(batch, 1)[1]),
       #     #('total_cells_without_beta', self.cells_with_alpha_beta(batch, 0)[1]),
       #     ('clonal_cells_with_beta', self.cells_with_alpha_beta(cl, 1)[1]),
       #     ('clonal_cells_without_beta', self.cells_with_alpha_beta(cl, 0)[1]),
       #     ('nonclonal_cells_with_beta', self.cells_with_alpha_beta(noncl, 1)[1]),
       #     ('nonclonal_cells_without_beta', self.cells_with_alpha_beta(noncl, 0)[1]),
       #     ('pvalue(difference in beta)', fischer(cl, noncl, 1)[1]),
        ]
        
        counts = self.sets2counts(OrderedDict(array))
        
        #fisher_pairs = self.cells_with_alpha(cl, noncl)
        fisher_pairs = self.cells_with_beta(cl, noncl)
        row = counts.append(pd.Series(OrderedDict(fisher_pairs)))
        row.name = self.batch.tag
        #d = OrderedDict(array + fisher_pairs)
        #row = pd.DataFrame(d, index=[self.batch.tag])#, name=batch.tag)
        
        Path(out_counts_fn).parent.mkdir(parents=True, exist_ok=True)
        Path(out_cells_fn).parent.mkdir(parents=True, exist_ok=True)
        
        print('Writing to ', out_counts_fn)
        print('Writing to ', out_cells_fn)
        counts.to_csv(out_counts_fn)
        #pd.Series(OrderedDict(array)).to_csv(out_fn)
        pd.Series(OrderedDict(array)).to_json(out_cells_fn)
        return row
