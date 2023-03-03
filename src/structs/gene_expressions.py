import numpy as np
import pandas as pd
import re

from utils import *

class GE:
    def __init__(self, ge, *, sample_id, cell2kallisto_subdir):
        """Class for working with quality controls per cell.
        Initialize using fromTable or fromDir.
        """
        self.ge = ge
        self.sample_id = sample_id
        self.cell2kallisto_subdir = cell2kallisto_subdir   # dict: cell -> subdir
        #if transform:
        #    if transform == 'log2':
        #        self.ge = self.log2_transform()
        #    else:
        #        self.ge = self.transform(transform)
        """A pandas table with cells as rows and transcript expressions as columns."""
        assert_no_NaNs(self.ge)
        #if not self.ge.empty:
        #    print(red('WARNING: ') + 'Empty GE')
        
    def to_file(self, output_file, sep='\t'):
        def strip_batch_name(cell):
            match = re.search(r'.*-(\d+)$', cell)  # the number after the last '_'
            if match:
                return match.groups()[0]
            else:
                if not re.match(r'(\d+)$', cell):
                    print('Bad cell name: ', cell)
                    assert False
        #stripped_index = [ strip_batch_name(cell) for cell in self.ge.index ]
            
        print('Writing to ' + str(output_file))
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        renamed_T = self.ge.T.rename(columns=strip_batch_name, inplace=False).sort_index(axis=1)
        renamed_T.to_csv(output_file, sep=sep, float_format='%g')
        
    def __repr__(self):
        return '{} cells x {} genes'.format(self.ge.shape[0], self.ge.shape[1])

    @staticmethod
    def cell2file(cell, *, fns, sample_id):
        def suffix(a):
            return a.split('_')[-1].split('.')[-1]

        for f in cells:
            ending = f[-9:]
            if ending[0] == 'N' and ending[5] == 'S':
                ending = ending.replace('.', '_')
            else:
                num = ''
                for c in ending:
                    if c.isdigit():
                        num += c
                    else:
                        num = ''
                ending = num
            ending = '{}_{}'.format(sample_id, ending)
            #print('cell2file: ', cell, ending)
            #if ending == cell:
            #print(ending.split('_'), cell.split('_'))
            if suffix(ending) == suffix(cell):
                return f
        warning('No cell {} found in sample {}.'.format(cell, sample_id))
        return None


    @classmethod
    def from_file(cls, true_index, input_file, *, sample_id, kallisto_dir, sep='\t',
                  filt_cells, expression_threshold=None, min_cells_expressing_gene=None,
                  #transform='log2',
                  nrows=None):
        if not input_file.is_file():
            warning('The expressions file {} does not exist.'.format(blue(input_file)))
        #    in_files = list(input_file.parent.glob(input_file.stem + '*'))
        #    if len(in_files) > 1:
        #        assert not check_file(input_file, "Transcript expressions")
        #        return None
        #        #raise OSError("The TCR file %s doesn't exist and more multiple files with the prefix %s exist: %s",
        #        #              input_file, input_file.stem(), ', '.join(in_files))
        #    elif len(in_files) < 1:
        #        assert not check_file(input_file, "Transcript expressions")
        #        return None
        #    else:
        #        input_file = in_files[0]

        assert check_file(input_file, "Transcript expressions")

        if nrows:
            print(red('Warning: ') + 'Only the first {} genes are loaded.'.format(red(nrows)))
            print(red('Warning: ') + 'The filtering options may be stupid!')  ## TODO: remove 
        # Quantification table from Kallisto
        quant = pd.read_csv(input_file, index_col=0, sep=sep, nrows=nrows).T
        quant = quant.astype(int)
        assert_no_NaNs(quant)
        prev_shape = quant.shape

        # filter by spikein content %
        #warning('TODO: gene expressions.py: Spikein filtration')
        #print(list(quant)[-10:])
        spikein_counts = quant[ ['Spike1,', 'Spike4,', 'Spike7,'] ].sum(axis=1)
        spike_perc = spikein_counts.div(quant.sum(axis=1))
        #print('spike_perc: ', spike_perc)
        #print('before: ', quant.shape)
        #quant = quant.loc[ spike_perc < 0.9 ]  # TODO: set a param
        if 'max_spikein_perc' in filt_cells:
            filtered_cells = quant.loc[ spike_perc >= filt_cells['max_spikein_perc'] ].index  # TODO: set a param
        else:
            filtered_cells = []
        #print('after: ', quant.shape)
        
        if expression_threshold and min_cells_expressing_gene:
            quant = quant.loc[:, (quant > expression_threshold).sum(axis='index') >= min_cells_expressing_gene]  # delete the zero-columns/non-expressed  transcripts
        else:
            None
            #warning('No filtering! Please provide expression_threshold and min_cells_expressing_gene.')
        #if prev_shape != quant.shape:
        #    output_file = input_file.with_name(input_file.name + '_cleaned')
        #    print(red('WARNING: ') + 'Removed {}/{} genes with non-zero quantities in less then {} cells'.format(
        #        blue(prev_shape[1]-quant.shape[1]), blue(prev_shape[1]), green(min_cells_expressing_gene)))
        #    quant.T.to_csv(output_file, sep=sep, float_format='%g')
                
        #quant = np.log2(1 + quant)
        quant.index.name = 'cell'
        quant.sort_index(inplace=True)
        orig_index = quant.index

        def fn2cell(cell):
            bla = cell[-9:]
            if bla[0] == 'N' and bla[5] == 'S':
                return bla.replace('.', '_')
            else:
                num = ''
                for c in bla:
                    if c.isdigit():
                        num += c
                    else:
                        num = ''
                return num

        cell2kallisto_subdir = {}
        for fn in orig_index.tolist():
            #fn = GE.cell2file(cell, fns=orig_index.tolist(), sample_id=sample_id)
            cell = sample_id + '_' + fn2cell(fn)
            cell2kallisto_subdir[cell] = str(Path(kallisto_dir) / fn) if fn else None
        
        add_prefix_tag(quant, sample_id)
        assert_same_indexes(true_index, quant.index, 'expressions')
        filtered_cells = [ sample_id + '_' + fn2cell(fn) for fn in filtered_cells ]

        # erase all expression data in order to lighten the struct
        quant.drop(labels=quant.columns, axis='columns', inplace=True)

        return GE(quant, sample_id=sample_id, cell2kallisto_subdir=cell2kallisto_subdir), orig_index, filtered_cells
    
    #def transform(self, f):
    #    return f(self.ge)
    #
    #def log2_transform(self):
    #    print("Log2 transformed expressions!")
    #    return self.transform(lambda ge: np.log2(1+ge))
    
    def subge(self, cells, genes_list=None):
        intersect = set(cells).intersection(self.cells())
        assert_same_indexes(cells, intersect)
        if not genes_list is None:
            avail_genes = set(self.ge.columns)
            intersect_genes = avail_genes & set(genes_list)
            print('Selecting {} genes ({} not found) out of {}'.format(
                blue(len(intersect_genes)), blue(len(set(genes_list) - intersect_genes)), blue(len(avail_genes))))
            return GE(self.ge.loc[list(intersect), list(intersect_genes)])
        return GE(self.ge.loc[list(intersect)])
    
    def cells(self):
        return self.ge.index.tolist()

    def cell2kallisto(self, cell):
        d = self.cell2kallisto_subdir.get(cell)
        if not d:
            warning('No kallisto directory for cell {} in sample {}.'.format(cell, self.sample_id))
        return d
    
    def tids2df(self, transcripts):
        cols = list(set(self.ge.columns) & set(transcripts))
        return self.ge[cols]
    
    def transcripts(self):
        return self.ge.columns
    
    def concat(self, other, sample_id):
        conc = pd.concat([self.ge, other.ge])
        cell2kallisto_subdir = {**self.cell2kallisto_subdir, **other.cell2kallisto_subdir}  # if cell else None
        # TODO: print how many transcripts were filled with NaN
        return GE(conc.fillna(0.0), sample_id=sample_id, cell2kallisto_subdir=cell2kallisto_subdir)
    
    def print_stats(self):
        print("{} cells by {} transcripts".format(blue(self.ge.shape[0]), blue(self.ge.shape[1])))
