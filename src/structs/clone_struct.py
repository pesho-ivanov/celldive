from collections import OrderedDict
from numpy import mean, log2, exp2
from decimal import Decimal
from itertools import combinations

from utils import *

#import rpy2
#import rpy2.robjects as robj
#from rpy2.robjects.packages import importr
#r = robj.r
#base = importr('base')
#stats = importr('stats')
#utils = importr('utils')

#class Materials:
#    Shaded, Shiny, Transparent, Matte = range(4)
#class CellColor(Enum):
# Usage: Batch.CellColor.NONRECONSTRUCTED
class CellPart:
    CLONAL            = 'main'     # by clonality rule
    RELATED           = 'main+'    # as being related to the clone
    AMBIGUOUS         = 'big'
    BYSTANDERS        = 'small'   # by being part of a non-main clone
    NONCLONAL         = 'single'  # by not beging part of any clone

    FILTERED          = 'filtered'   # by expression QC, number of cells in a well, etc.
    NONRECONSTRUCTED  = 'nonreconstructed'  # by not having TCR info

    EARLY_CLONAL      = 'earlyclonal'
    LATE_CLONAL       = 'lateclonal'
    EARLY_NONCLONAL   = 'earlysingle'
    LATE_NONCLONAL    = 'latesingle'
#    @staticmethod
#    def str2CellPart(part):
#        [a for a in dir(obj) if not a.startswith('__') and not callable(getattr(obj,a))] 

CellColor = {
    CellPart.NONRECONSTRUCTED: 'light grey',
    CellPart.FILTERED:         'yellow',

    CellPart.CLONAL:           'red',
    CellPart.RELATED:          'light orange',
    CellPart.AMBIGUOUS:        'light green',
    CellPart.BYSTANDERS:       'light blue',
    CellPart.NONCLONAL:        'royal blue',

    CellPart.EARLY_CLONAL:     'light pink',
    CellPart.LATE_CLONAL:      'pink',
    CellPart.EARLY_NONCLONAL:  'light blue',
    CellPart.LATE_NONCLONAL:   'blue',
}

class SampleColor:
    SAMPLE = 'sample'  # one color per sample
    CLONAL = 'clonal'  # fixed clonal/non-clonal/unknown colors for all batches
    SAMPLE_AND_CLONAL = 'sampleAndClone'  # with shades for clonal/non-clonal/unknown
    RAINBOW = 'rainbow'  # with shades for clonal/non-clonal/unknown

class CellPartition:
    """Container for partitioning all cells of a batch into clonal-nonclonal-unknown."""

    def __init__(self, clonotype, part_dict, *, sample_id, theta, max_tc_size, cell_files):
        self.sample_id = sample_id
        self.part_dict = part_dict  # dict(CellPart -> set of cells)
        self.clonotype = OrderedDict(clonotype)  # chain -> tcr_seq
                                                 # param -> val
                                                 # WRONG: dict('alpha' -> dict(tcr_str -> set of cells))
                                                 # dict('beta' -> dict(tcr_str -> set of cells))
                                                 # dict('rule' -> str from {'alpha', 'beta', 'alpha or beta', 'alpha and beta'})
        self.theta = theta
        self.max_tc_size = max_tc_size
        self.cell_files = cell_files  # pandas index

        if not self.__is_correct():
            raise Exception("Not correct cloanlity partition")

    def cells(self):
        return set.union(*self.part_dict.values())

    #def chain2tcrs(self, chain):
    #    print(self.clonotype)
    #    return set(self.clonotype[chain].keys())

    #def cell2seqs(self, cell, *, chains=['alpha', 'beta']):
    #    seqs = [ seq for chain in chains for seq, cells in self.part_dict[chain].items() if cell in cells ]
    #    seqs = [ seq for seq in seqs if not seq.startswith(TCRs.NONPRODUCTIVE_TAG)]
    #    return seqs

    def get_group_sizes(self, *, clear_groups):
        if CellPart.CLONAL in self.part_dict and CellPart.NONCLONAL in self.part_dict:
            x = len(self.part_dict[CellPart.CLONAL])
            y = len(self.part_dict[CellPart.NONCLONAL])
            if not clear_groups:
                x += len(self.part_dict[CellPart.RELATED])
                y += len(self.part_dict[CellPart.BYSTANDERS])
        elif CellPart.EARLY_CLONAL in self.part_dict and CellPart.LATE_CLONAL in self.part_dict:
            x = len(self.part_dict[CellPart.EARLY_CLONAL])
            y = len(self.part_dict[CellPart.LATE_CLONAL])
        else:
            assert False, 'Not enough groups for DE'

        return (x, y)

    def __is_correct(self):
        #print('clone_struct: ', self.part_dict.values())
        for a, b in combinations(self.part_dict.items(), 2):
            groupname_a, cells_a = a
            groupname_b, cells_b = b
            if cells_a & cells_b:
                warning('Groups {} and {} intersect.'.format(groupname_a, groupname_b))
                return False
        return True
        #union = set()
        #sum_sets = 0
        #for v in self.part_dict.values():
        #    sum_sets += len(v)
        #    union |= v
        #if sum_sets != union:
        #    warning("sum_sets={} != union={}".format(sum_sets, union))
        #return sum_sets == union

    def __repr__(self):
        return str(self.part_dict)
        #cl    = 'clonal: {}'.format(red(', '.join(sorted(self.clonal))))
        #noncl = 'nonclonal: {}'.format(blue(', '.join(sorted(self.nonclonal))))
        #unkn  = 'unknown: {}'.format(', '.join(sorted(self.unknown)))
        #return ' | '.join([cl, noncl, unkn])

    @staticmethod
    def __tcrs2clonal(tcrs, clonotype):
        alpha, beta, rule = clonotype['alpha'], clonotype['beta'], clonotype['rule'], 

        alpha_cells = set(tcrs.select_cells(alpha))
        beta_cells = set(tcrs.select_cells(beta))

        if 'alpha' in rule and not alpha:
            raise Exception('No alpha clone specified.')
        if 'beta' in rule and not beta:
            raise Exception('No beta clone specified.')

        if rule == 'alpha and beta':
            clonal_cells = alpha_cells & beta_cells
        elif rule == 'alpha or beta':
            clonal_cells = alpha_cells | beta_cells
        elif rule == 'alpha':
            clonal_cells = alpha_cells
        elif rule == 'beta':
            clonal_cells = beta_cells
        else:
            raise Exception('Incorrect clonality rule {}.'.format(rule))
        
        return clonal_cells

    @staticmethod
    def __get_ambiguous(*, remaining_cells, tcrs, clonal_cells):
        '''Finds cells which are on dist=1 from clonal.'''
        clonal_tcrs = set()
        for cell in clonal_cells:
            A = tcrs.cell2seqs(cell, chains=['alpha'])
            B = tcrs.cell2seqs(cell, chains=['beta'])
            clonal_tcrs |= set(A) | set(B)
            
        tcr2cells = tcrs.tcr2cells_dict()     

        ambiguous_cells = set()
        for cell in remaining_cells:
            A = tcrs.cell2seqs(cell, chains=['alpha'])
            B = tcrs.cell2seqs(cell, chains=['beta'])
            AB = set(A) | set(B)
            if AB & clonal_tcrs:   # intersects with the clonal alpha or beta
                ambiguous_cells.add(cell)
               
        return ambiguous_cells

    @classmethod
    def __get_nonclonal_cells(cls, *, remaining_cells, tcrs, clonal_cells, bystander_clone_threshold_percent):
        absolute_threshold = len(remaining_cells | clonal_cells) * bystander_clone_threshold_percent / 100.0
        nonclonal_cells, bystander_cells = set(), set()
        max_tc_size = 0
        for cell in remaining_cells:
            tc = cls.__transitive_closure(remaining_cells=remaining_cells, tcrs=tcrs, seed_cells={cell})
            max_tc_size = max(max_tc_size, len(tc))
            if len(tc) == 1:
                nonclonal_cells.add(cell)
            elif len(tc) <= absolute_threshold:
                bystander_cells.add(cell)
        return nonclonal_cells, bystander_cells, max_tc_size

    # DEPRECATED
    @staticmethod
    def __get_nonclonal_cells_DEPRECATED(*, remaining_cells, tcrs, clonal_cells, bystander_clone_threshold_percent):
        def get_max_nonclone_size(out=False):
            # TO DEBUG on SchH
            #for cell in (set(self.batch.tcrs.cells()) - clonal_tcrs):
            reconstructed_cells = len(tcrs.cells())
            #perc = 0.05
            #max_nonclone_size = int(perc * reconstructed_cells) + 1
            max_nonclone_size = 1
            if out:
                if reconstructed_cells > 0:
                    print('maximum size of nonclones for {} is {}/{}={} (target: {})'.format(
                        batch.tag, blue(max_nonclone_size), blue(reconstructed_cells),
                        blue("%.2f%%" % (100.0*max_nonclone_size / reconstructed_cells)),
                        blue("%.2f%%" % (100.0*perc))))
            return max_nonclone_size

        clonal_tcrs = set()
        for cell in clonal_cells:
            A = tcrs.cell2seqs(cell, chains=['alpha'])
            B = tcrs.cell2seqs(cell, chains=['beta'])
            clonal_tcrs |= set(A) | set(B)
            
        tcr2cells = tcrs.tcr2cells_dict()     

        nonclonal_cells = set()
        for cell in remaining_cells:
            A = tcrs.cell2seqs(cell, chains=['alpha'])
            B = tcrs.cell2seqs(cell, chains=['beta'])
            AB = set(A) | set(B)
            if AB:   # at least alpha or beta  <=> reconstructed
                if not (AB & clonal_tcrs):   # no intersection with the clonal alpha or beta <=> no incidence to the main clone
                    #tcrs = [tcr for tcr in AB if len(tcr2cells[tcr]) > 1 ]
                    _tcrs = [tcr for tcr in AB if len(tcr2cells[tcr]) > get_max_nonclone_size() ] # clones with alpha or beta   <=> no part of other (big) clones
                    if not _tcrs:   # no bad incident tcrs
                        nonclonal_cells.add(cell)
               
        return nonclonal_cells

#    @classmethod
#    def from_clonality_rule(cls, *, all_cells, tcrs, clonotype, filtered, bystander_clone_threshold_percent):
#        """
#        Parameters
#        ----------
#        cells : iterable object
#            Cell names.
#        tcrs : TCRs object
#            Should include information about alpha and beta chains per cell.
#        clonotype: dict (alpha, beta, clonality_rule, facs_sorting -> string)
#            alpha, beta : Shortened notations of alpha or beta chain (e.g. "TRAV6_TGTGCCGTACGGGATA_TRAJ12").
#            clonality_rule : One of "alpha", "beta, "alpha or beta" or "alpha and beta".
#        filtered: iterable of cell names that were prefiltered
#        """
#        if not tcrs:
#            raise Exception('tcrs not available')
#        assert bystander_clone_threshold_percent
#
#        tcr_cells = set(tcrs.cells())
#        all_cells = set(all_cells)
#        filtered = set(filtered)
#
#        bad_cells = tcr_cells - all_cells
#        nonreconstructed = all_cells - tcr_cells
#        if bad_cells:
#            raise Exception('tcr cells not among meta cells, e.g. {}'.format(', '.join(list(bad_cells[:3]))))
#
#        subtcrs = tcrs.subtcrs(cells = tcr_cells - filtered)
#        clonal = cls.__tcrs2clonal(subtcrs, clonotype)
#        remaining = tcr_cells - filtered - clonal
#        related = cls.__get_ambiguous(remaining_cells=remaining, tcrs=subtcrs, clonal_cells=clonal)
#        remaining -= related
#        nonclonal, max_tc_size = cls.__get_nonclonal_cells(remaining_cells=remaining, tcrs=tcrs, clonal_cells=clonal, 
#                bystander_clone_threshold_percent=bystander_clone_threshold_percent)
#        theta = 1
#        remaining -= nonclonal
#        bystanders = remaining
#
#        part_dict = {
#            CellPart.CLONAL:           clonal,
#            CellPart.RELATED:          related,
#            CellPart.BYSTANDERS:       bystanders,
#            CellPart.NONCLONAL:        nonclonal,
#
#            CellPart.FILTERED:         filtered,
#            CellPart.NONRECONSTRUCTED: nonreconstructed,
#        }
#
#        return cls(clonotype, part_dict, theta=theta, max_tc_size=max_tc_size)

### New transitive closure algorithm

    @classmethod
    def __transitive_closure(cls, *, remaining_cells, tcrs, seed_cells):
        '''Finds cells which are on dist=1 from clonal.'''

        closure_cells = seed_cells.copy()
        remaining_cells = remaining_cells.copy()
        tcr2cells = tcrs.tcr2cells_dict()     

        while True:  # add to `closure' one neighbour layer at a time
            seed_tcrs = set()
            for cell in closure_cells:
                A = tcrs.cell2seqs(cell, chains=['alpha'])
                B = tcrs.cell2seqs(cell, chains=['beta'])
                seed_tcrs |= set(A) | set(B)
        
            neighbour_cells = set()
            for cell in remaining_cells:
                A = tcrs.cell2seqs(cell, chains=['alpha'])
                B = tcrs.cell2seqs(cell, chains=['beta'])
                AB = set(A) | set(B)
                if AB & seed_tcrs:   # intersects with the clonal alpha or beta
                    neighbour_cells.add(cell)

            if neighbour_cells:
                closure_cells |= neighbour_cells
                remaining_cells -= neighbour_cells
                #print('Neighbour cells: ', neighbour_cells)
            else:
                break
               
        return closure_cells

    @classmethod
    def __calc_maximum_incidence(cls, *, cells, tcrs):
        theta = 0
        tcr2cells = tcrs.tcr2cells_dict()
        for cell in cells:
            A = tcrs.cell2seqs(cell, chains=['alpha'])
            B = tcrs.cell2seqs(cell, chains=['beta'])
            AB = set(A) | set(B)

            #print('tcrs of cell {}: '.format(cell), tcrs.cell2tcrs(cell, chain=None))

            if AB:
                incident_cell_nums = [ len(tcr2cells[chain]) for chain in AB ]
                #print('incident_cell_nums for {}: {}'.format(cell, incident_cell_nums))
                theta = max(theta, max(incident_cell_nums))
            #else:
            #Exception('No (productive) alpha or beta chains for cell {}'.format(cell))
               
        return theta


    @classmethod
    def TC_with_threshold(cls, *,
            sample_id, all_cells, tcrs, clonotype, filtered, bystander_clone_threshold_percent, cell_files):
        '''
        Let all cells (C) of sa sample are split to cell sets M, m, A, b, B,
           where M/m for Main clone, B/b stand for Bystanders, and A for ambiguous.
        M := M(clonality_rule)   // in red
        B := B(max_theta) for max_theta = argmax_theta,  // in blue
          s.t. 
           1) theta is the biggest number of incident cells of a bystander cell.
           2) B* \cap M* = \emptyset
        m := M* \ M, b := B* \ B  // in light shades of red and blue
        A := C \ M* \ B*          // in green

        Note: Asterisk (*) stands for transitive closure (TC) in the bipartite graph
        Note2: There are also Nonreconstructed and Filtered cells that we excluded.
        '''

        if not tcrs:
            raise Exception('tcrs not available')
        assert bystander_clone_threshold_percent is not None

        tcr_cells = set(tcrs.cells())
        all_cells = set(all_cells)
        filtered = set(filtered)

        bad_cells = tcr_cells - all_cells
        nonreconstructed = all_cells - tcr_cells - filtered
        if bad_cells:
            raise Exception('tcr cells not among meta cells, e.g. {}'.format(', '.join(list(bad_cells[:3]))))

        subtcrs = tcrs.subtcrs(cells = tcr_cells - filtered)
        main = cls.__tcrs2clonal(subtcrs, clonotype)
        remaining = tcr_cells - filtered - main - nonreconstructed
        maybe_main = cls.__transitive_closure(remaining_cells=remaining, tcrs=subtcrs, seed_cells=main) - main
        remaining -= maybe_main
        nonclonal, bystanders, max_tc_size = cls.__get_nonclonal_cells(remaining_cells=remaining, tcrs=tcrs,
                clonal_cells=main|maybe_main, bystander_clone_threshold_percent=bystander_clone_threshold_percent)
        theta = cls.__calc_maximum_incidence(cells=remaining, tcrs=tcrs)
        # TODO: ambiguous by bystander_clone_threshold_percent
        #print('theta', theta)
        remaining -= nonclonal | bystanders
        #bystanders = remaining   # TODO: change meaning
        ambiguous = remaining

        part_dict = {
            CellPart.FILTERED:         filtered,
            CellPart.NONRECONSTRUCTED: nonreconstructed,

            CellPart.CLONAL:           main,            # M
            CellPart.RELATED:          maybe_main,      # m
            CellPart.AMBIGUOUS:        ambiguous,       # A
            CellPart.BYSTANDERS:       bystanders,      # b
            CellPart.NONCLONAL:        nonclonal,       # B
        }

        return cls(clonotype, part_dict, sample_id=sample_id, theta=theta, max_tc_size=max_tc_size, cell_files=cell_files)

    def is_part_of(cell, clas):
        return cell in self.part_cells[clas]

    def get_cell_part(self, cell):
        for part, cells in self.part_dict.items():
            if cell in cells:
                return part
        return Exception('No such cell: {}'.format(cell))
    
#    def __init__(self, clonotype, filtered, clonal, nonclonal, unknown, nonreconstructed):
#        """
#        Parameters
#        ----------
#        clonotype : dict (alpha, beta, rule, facs_sorting -> string)
#            
#        clonal, nonclonal, unknown, nonreconstructed: iterable object with string elements
#            Set of cells in the corresponding class (e.g. in `clonal`).
#        """
#        self.filtered = set(filtered)
#        self.clonal = set(clonal)
#        self.nonclonal = set(nonclonal)
#        self.unknown = set(unknown)
#        self.nonreconstructed = set(nonreconstructed)
#        self.clonotype = clonotype
#        
#        union = len(self.filtered | self.clonal | self.nonclonal | self.unknown | self.nonreconstructed)
#        al = len(self.filtered) + len(self.clonal) + len(self.nonclonal) + len(self.unknown) + len(self.nonreconstructed)
#        assert  union == al, \
#                'union={}, all={}, filtered: {}, clonal: {}, nonclonal: {}, unknown: {}, nonreconstructed: {}'.format(
#               union, al, self.filtered, self.clonal, self.nonclonal, self.unknown, self.nonreconstructed)


    #@staticmethod
    #def get_cell_partition(batch):
    #    assert batch.tcr, 'No tcr for {}'.format(batch.sample_id)
        
        #assert not (self.clonal & self.nonclonal)
        #assert not (self.nonclonal & self.unknown)
        #assert not (self.unknown & self.clonal)
        #assert not (self.nonreconstructed & self.clonal)

    #def __or__(self, b):
    #    if self.clonotype != b.clonotype:
    #        Warning('Different clonotypes: {} and {}'.format(self.clonotype, b.clonotype))
    #    return CellPartition(self.clonotype,
    #                         self.clonal | b.clonal,
    #                         self.nonclonal | b.nonclonal,
    #                         self.unknown | b.unknown,
    #                         self.nonreconstructed | b.nonreconstructed)
    #
    #def __sub__(self, cells):
    #    return CellPartition(self.clonotype,
    #                         self.clonal - set(cells),
    #                         self.nonclonal - set(cells),
    #                         self.unknown - set(cells),
    #                         self.nonreconstructed - set(cells))
    #
    #def __and__(self, cells):
    #    return CellPartition(self.clonotype,
    #                         self.clonal & set(cells),
    #                         self.nonclonal & set(cells),
    #                         self.unknown & set(cells),
    #                         self.nonreconstructed & set(cells))

    #@staticmethod
    #def merge_sets_in_dict(X, Y):
    #    # merging all the sets in the dict
    #    return { key: { *X.get(key, {}), *Y.get(key, {}) } for key in {*X.keys(), *Y.keys()} }
    
    # WARNING: strange operation
    def earlyVsLateClonal(self, b):
        assert self.clonotype['alpha'] == b.clonotype['alpha']
        assert self.clonotype['beta'] == b.clonotype['beta']
        assert self.clonotype['rule'] == b.clonotype['rule']

        return CellPartition(clonotype = self.clonotype,
                             #clonotype = { 
                             #    'alpha': self.merge_sets_in_dict(self.clonotype['alpha'], b.clonotype['alpha']),
                             #    'beta': self.merge_sets_in_dict(self.clonotype['beta'], b.clonotype['beta']),
                             #    'rule': rule },
                             part_dict = { CellPart.EARLY_CLONAL: self.part_dict[CellPart.CLONAL],
                                           CellPart.LATE_CLONAL: b.part_dict[CellPart.CLONAL] },
                             sample_id = 'clonal_' + self.sample_id + '_and_clonal_' + b.sample_id,
                             theta = 0,  # TODO: remove from __init__
                             max_tc_size = 0,
                             cell_files = self.cell_files.union(b.cell_files)  # index merging
                            )

    def earlyVsLateSingle(self, b):
        assert self.clonotype['alpha'] == b.clonotype['alpha']
        assert self.clonotype['beta'] == b.clonotype['beta']
        assert self.clonotype['rule'] == b.clonotype['rule']

        return CellPartition(clonotype = self.clonotype,
                             part_dict = { CellPart.EARLY_NONCLONAL: self.part_dict[CellPart.NONCLONAL],
                                           CellPart.LATE_NONCLONAL: b.part_dict[CellPart.NONCLONAL] },
                             sample_id = 'single_' + self.sample_id + '_and_single_' + b.sample_id,
                             theta = 0,  # TODO: remove from __init__
                             max_tc_size = 0,
                             cell_files = self.cell_files.union(b.cell_files)  # index merging
                            )

    def merge(self, b):
        assert self.clonotype['alpha'] == b.clonotype['alpha']
        assert self.clonotype['beta'] == b.clonotype['beta']
        assert self.clonotype['rule'] == b.clonotype['rule']

        return CellPartition(clonotype = self.clonotype,
                             part_dict = { **self.part_dict, **b.part_dict },
                             sample_id = 'all_' + self.sample_id + '_and_all_' + b.sample_id,
                             theta = 0,  # TODO: remove from __init__
                             max_tc_size = 0,
                             cell_files = self.cell_files.union(b.cell_files)  # index merging
                            )
    
#class Hypothesis:
#    # TODO: should preserve old p-values when subbatching
#    def __init__(self, batch, *, nonclonal_cells, clonal_cells, alpha=None, beta=None, logic=None,
#                 diffexpr, best_genes=None):
#        self.batch = batch
#        #self.cells = self.batch.tcrs.cells()
#        self.cells = self.batch.cells()
#        self.alpha = alpha
#        self.beta = beta
#        self.logic = logic
#        self.nonclonal_cells = nonclonal_cells    # TODO: move away to CellPartition
#        self.clonal_cells = clonal_cells          # TODO: move away to CellPartition
#        self.diffexpr = diffexpr
#        
#        if batch.quant is not None:
#            if best_genes:
#                self.best_genes = best_genes
#            else:
#                self.best_genes = self.calc_best_genes(batch.genes().tolist(), diffexpr)
#        
#        #if out_fn:
#        #    self.print_best_genes(out_fn)
#            
#        assert not (set(self.clonal_cells) & set(self.nonclonal_cells))
#        assert set(self.clonal_cells) | set(self.nonclonal_cells) <= set(self.cells)
#        
#    #def get_best_genes(self, diffexpr):
#    #    if self._best_genes is None:
#    #        assert diffexpr
#    #        self._best_genes = self.calc_best_genes(self.batch.genes().tolist(), diffexpr)
#    #    return self._best_genes
#        
#    def clone_tag(self):
#        return '{} {} {}'.format(self.alpha, self.beta, self.logic)
#        
#    def __repr__(self):
#        res = ''
#        res += '{} ({} cells vs {} cells)\n'.format(
#            blue(self.clone_tag()), green(len(self.nonclonal_cells)), green(len(self.clonal_cells)))
#        res += self.get_formatted_genes()
#        return res
#   
#    def print_best_genes(self, out_fn):
#        best_genes = self.best_genes
#        print('Writing {} best genes for {} to {}'.format(
#            blue(len(best_genes)), self.batch.tag, str(Path(out_fn).name)))
#        
#        #df = pd.DataFrame(self.best_genes, columns=['padj'])
#        #df.rename_axis("transcript_id")
#        
#        Path(out_fn).parent.mkdir(parents=True, exist_ok=True)
#        best_genes.to_csv(str(out_fn))
#        
#        df_styler = best_genes.rename_axis(None).style.format({
#            'log2fc_ratio': '{:.2}',
#            'log2fc_diff': '{:.2}',
#            'p_adj': '{:.2%}',
#        })
#        
#        html = sorttable_html() + df_styler.render()
#        #html = html.replace('<table', '<table class="sortable"')
#        html = html.replace('class="dataframe"', 'class="sortable"')
#
#        html_file = open(str(out_fn) + '.html', "w")
#        html_file.write(html)
#        html_file.close()
#        
#        #df = pd.DataFrame(best_genes, columns=['padj'])
#        #df.rename_axis("transcript_id")
#        #df['transcript_name'] = list(map(lambda t: tid2tname(t), df.index))
#        #df.to_csv(fn)
#    
#    def get_formatted_genes(self):
#        TO_SHOW = 30
#        res = ''
#        p_col = 'p_adj({})'.format(self.diffexpr)
#        best_genes = self.best_genes
#        best_genes.sort_values(by=p_col, inplace=True)
#        for t, p_adj in best_genes[:TO_SHOW][p_col].items():
#            if p_adj >= 1.0:
#                break
#
#            p_str = blue('{:>9}'.format(pval2str(p_adj)))
#
#            if p_adj >= MAX_P_VALUE:
#                if p_adj < 1.0:
#                    res += '    {} and worse...\n'.format(p_str)
#                break
#            else:
#                res += '    {} {} {}\n'.format(p_str, t, green(tid2tname(t)))
#        cnt = 0
#        for t, p_adj in best_genes[TO_SHOW+1:][p_col].items():
#            if p_adj <= MAX_P_VALUE:
#                cnt += 1
#        if cnt > 0:
#            res += '    ... [%s more cells with p_value < %.2f]\n' % (blue(cnt), MAX_P_VALUE)
#        return res
#    
#    #def number_of_good_genes(self):
#    #    return len(self.get_best_genes(self.diffexpr)[self.best_genes[p_col] < MAX_P_VALUE]
#    
#    def number_of_best_genes(self):
#        p_col = 'p_adj({})'.format(self.diffexpr)
#        return len(self.best_genes[p_col] < MAX_P_VALUE)
#        
#    def __lt__(self, other):
#        assert self.diffexpr == other.diffexpr
#        return number_best_genes(self) > number_best_genes(other)
#        
#    def calc_best_genes(self, genes, diffexpr):
#        log2fc_ratio, log2fc_diff = self.calc_log2fc(genes)
#        df = pd.DataFrame({
#            'log2fc_ratio': log2fc_ratio,
#            'log2fc_diff': log2fc_diff,
#        })
#        
#        xy2p = diffexpr2func[diffexpr]
#        max_p_value = diffexpr2max_p[diffexpr]
#        col = 'p_adj({})'.format(diffexpr)
#        df[col] = self.calc_diff_expr(genes, xy2p, max_p_value)
#        
#        df.dropna(axis='index', how='any', inplace=True)
#        df['transcript_name'] = self.get_transcript_names(df.index)
#        df.index.name = 'transcript_id'
#        return df
#    
#    def get_transcript_names(self, tids):
#        return pd.Series(list(map(lambda t: tid2tname(t), tids)), index=tids)
#    
#    # Bioconductor RankProd Package
#    def calc_diff_expr(self, genes, xy2p, max_p_value):
#        if min(len(self.nonclonal_cells), len(self.clonal_cells)) < MIN_CELLS_IN_CLONE:
#            return pd.Series()
#        
#        #batches = [ self.batch.subbatch(cells, tag=self.batch.tag)
#        quants = [ self.batch.quant.subge(cells) for cells in [self.nonclonal_cells, self.clonal_cells] ]
#        
#        p_vals = pd.Series(index=self.batch.quant.ge.columns)
#        for gene in genes:
#            x, y = [ robj.FloatVector(quant.ge[gene]) for quant in quants ]
#            p_vals[gene] = xy2p(x,y)
#
#        p_adjust = robj.r['p.adjust'](p_vals.tolist(), method='fdr')
#        res = pd.Series(p_adjust, index=p_vals.index).where(lambda x : x!=1).dropna()
#        return res[res < max_p_value].sort_values()
#    
#    def calc_log2fc(self, genes):
#        # log2 fold-change := log2(E(clonal) / E(nonclonal))
#        quants = [ self.batch.quant.subge(cells) for cells in [self.nonclonal_cells, self.clonal_cells] ]
#        log2fc_ratio = pd.Series(index=self.batch.quant.ge.columns)
#        log2fc_diff = pd.Series(index=self.batch.quant.ge.columns)
#        for gene in genes:
#            x, y = [ quant.ge[gene] for quant in quants ]
#            x = [a for a in x if a > 0]
#            y = [a for a in y if a > 0]
#            if len(x) > 0 and len(y) > 0:
#                #log2fc_ratio[gene] = log2(mean((x+1)/(y+1)))
#                log2fc_ratio[gene] = mean(log2(x)) - mean(log2(y))  # Limma style: https://www.biostars.org/p/100460/ 
#                #print(log2fc_ratio[gene])
#                #log2fc_ratio[gene] = log2(mean(exp2(x))) - log2(mean(exp2(y)))
#                log2fc_diff[gene]  = mean(x) - mean(y) # <=> log2(mean(exp2(x)) - mean(exp2(y)))
#            else:
#                log2fc_ratio[gene] = np.nan
#                log2fc_diff[gene]  = np.nan
#            
#        return log2fc_ratio, log2fc_diff
#        
