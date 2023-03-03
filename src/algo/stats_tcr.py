from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np

from utils import *

# from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def levenshtein(source, target):
    if len(source) < len(target):
        return levenshtein(target, source)

    # So now we have len(source) >= len(target).
    if len(target) == 0:
        return len(source)

    # We call tuple() to force strings to be used as sequences
    # ('c', 'a', 't', 's') - numpy uses them as values by default.
    source = np.array(tuple(source))
    target = np.array(tuple(target))

    # We use a dynamic programming algorithm, but with the
    # added optimization that we only need the last two rows
    # of the matrix.
    previous_row = np.arange(target.size + 1)
    for s in source:
        # Insertion (target grows longer than source):
        current_row = previous_row + 1

        # Substitution or matching:
        # Target and source items are aligned, and either
        # are different (cost of 1), or are the same (cost of 0).
        current_row[1:] = np.minimum(
                current_row[1:],
                np.add(previous_row[:-1], target != s))

        # Deletion (target grows shorter than source):
        current_row[1:] = np.minimum(
                current_row[1:],
                current_row[0:-1] + 1)

        previous_row = current_row

    return previous_row[-1]

class StatsTCR:
    def __init__(self, batch):
        self.batch = batch
        
    '''Returns a dict with (tcr1, tcr2) -> dist(tcr1, tcr2)'''
    @staticmethod
    def calc_all_edit_distances(tcrs):
        return { pair: levenshtein(*pair) for pair in combinations(tcrs, 2) }
    
    @staticmethod
    def plot_hist(ax, dists, title, ylabel):
        ax.set_title(title)
        ax.hist(dists.values(), bins=range(30))
        ax.set_ylabel(ylabel);
        ax.set_xlabel('edit distance');
        
    def plot_edit_distance_histograms(self, out_fn=None, interactive=False):
        tcrs = self.batch.tcrs
        fig, ((ax_A, ax_B), (ax_C, ax_D)) = plt.subplots(2, 2)
        fig.suptitle(self.batch.tag)
        for chain, ax in [('alpha', ax_A), ('beta', ax_B)]:
            tcrs_list = tcrs.chain2tcrs(chain)
            title = '{} {} tcrs'.format(len(tcrs_list), chain)
            ylabel = '#{} TCR pairs'.format(chain)
            dists = StatsTCR.calc_all_edit_distances(tcrs_list)
            StatsTCR.plot_hist(ax, dists, title, ylabel)
            
        for chain, ax in [('alpha', ax_C)]: #, ('beta', ax_D)]:
            #print('tag', self.batch.onlytag)
            def dist_between_cells(c1, c2):
                #print(tcrs.cell2seqs(c1))
                dists = [ levenshtein(tcr1, tcr2)
                    for tcr1 in tcrs.cell2seqs(c1)
                    for tcr2 in tcrs.cell2seqs(c2) ]
                if dists and min(dists) > 0:
                    return min(dists)
                else:
                    return None

            title = 'Cells histogram'
            ylabel = '# cell pairs'
            dists = { pair: dist_between_cells(*pair) for pair in combinations(tcrs.cells(), 2) }
            dists = { pair: d for pair, d in dists.items() if d }
            StatsTCR.plot_hist(ax, dists, title, ylabel)
        
        plot(plt, interactive, out_fn)
            
        #for cell in tcrs.cells():
        #    A, B = tcrs.cell2tcrs(cell)
        #    calc_all_edit_distances(A)
        #    calc_all_edit_distances(B)
            
