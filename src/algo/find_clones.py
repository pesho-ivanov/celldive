from collections import Counter
from collections import defaultdict

from structs.clone_struct import Hypothesis
from utils import *
    
mem_hypo = {}  # memoized hypotheses

class FindClones:
    def __init__(self, batch, onlytag, hypotheses=None, out_fn=None, *, diffexpr):
        assert onlytag is not ''
        print('find_clones, onlytag=', onlytag)
        
        #self.out_fn = '../out/{}_genes_{}.csv'.format(self.batch.tag, self.clone_tag())
        self.out_fn = out_fn
        self.batch = batch
        self.onlytag = onlytag
        self.diffexpr = diffexpr
        
        if hypotheses:
            self.hypotheses = hypotheses
        else:
            de = diffexpr if diffexpr else ['no_de']
            if batch.cell_part:
                genes_num = len(batch.genes())
            else:
                genes_num = 0
            args = (batch.tag, genes_num, len(batch.cells()), batch.onlytag, *de)
            #if args in mem_hypo:
            #    self.hypotheses = mem_hypo[args]
            #else:
            self.hypotheses = {}
            assert onlytag, 'The onlytag should be "all" or a specific clonotype'
            self.fill_clone_hypotheses(onlytag)
            mem_hypo[args] = self.hypotheses
        # self.union_important_genes = None # TODO: filter before instead of self._get_union_important_genes()
    
    def subclones(self, parent_batch, cells):
        hypotheses = {}
        for h in self.hypotheses.values():
            cl = list(set(h.clonal_cells) & set(cells))
            noncl = list(set(h.nonclonal_cells) & set(cells))
            hypotheses[h.clone_tag()] = Hypothesis(
                parent_batch, nonclonal_cells=noncl, clonal_cells=cl,
                alpha=h.alpha, beta=h.beta, logic=h.logic, #h.get_best_genes(self.diffexpr), 
                diffexpr=self.diffexpr)
            #hypotheses[h.clone_tag()].print_best_genes(out_fn=str(self.out_fn)+'.sub')
        return FindClones(parent_batch, self.onlytag, hypotheses, diffexpr=self.diffexpr)
        
    def best_hypo(self):
        if len(self.hypotheses) == 0:
            return None
        return next(iter(self.hypotheses.values()))
        
    def __getitem__(self, hypothesis_tag):
        if not hypothesis_tag in self.hypotheses:
            self.fill_clone_hypotheses(hypothesis_tag)
        return self.hypotheses[hypothesis_tag]
        
    def __repr__(self):
        if len(self.hypotheses) >= 1:
            return 'best hypothesis: {}\n'.format(self.best_hypo().__repr__()).strip()
        return ''
        
    def _get_union_important_genes(self):
        important_genes = defaultdict(lambda: 1.0)
        for tag, h in self.hypotheses.items():
            p_col = 'p_adj({})'.format(self.diffexpr[0])
            for gene, padj in h.best_genes[p_col].iteritems():
                important_genes[gene] = min(padj, important_genes[gene])
        s = pd.Series(important_genes)
        max_p_value = diffexpr2max_p[self.diffexpr]
        return s[s < max_p_value].sort_values()
    
    def get_max_nonclone_size(self, *, out=False):
        # TO DEBUG on SchH
        #for cell in (set(self.batch.tcrs.cells()) - clonal_tcrs):
        perc = 0.05
        reconstructed_cells = len(self.batch.tcrs.cells())
        max_nonclone_size = int(perc * reconstructed_cells) + 1
        if out:
            if reconstructed_cells > 0:
                print('maximum size of nonclones for {} is {}/{}={} (target: {})'.format(
                    self.batch.tag, blue(max_nonclone_size), blue(reconstructed_cells),
                    blue("%.2f%%" % (100.0*max_nonclone_size / reconstructed_cells)),
                    blue("%.2f%%" % (100.0*perc))))
        return max_nonclone_size
    
    def get_nonclonal_cells(self, clonal_cells):
        clonal_tcrs = set()
        for cell in clonal_cells:
            A = self.batch.tcrs.cell2seqs(cell, 'alpha')
            B = self.batch.tcrs.cell2seqs(cell, 'beta')
            clonal_tcrs |= set(A) | set(B)
            
        tcr2cells = self.batch.tcrs.tcr2cells_dict()     
        
        nonclonal_cells = []
        for cell in self.batch.tcrs.cells():
            A = self.batch.tcrs.cell2seqs(cell, 'alpha')
            B = self.batch.tcrs.cell2seqs(cell, 'beta')
            AB = set(A) | set(B)
            if AB:   # at least alpha or beta  <=> reconstructed
                if not (AB & clonal_tcrs):   # no intersection with the clonal alpha or beta <=> no incidence to the main clone
                    #tcrs = [tcr for tcr in AB if len(tcr2cells[tcr]) > 1 ]
                    tcrs = [tcr for tcr in AB if len(tcr2cells[tcr]) > self.get_max_nonclone_size() ] # clones with alpha or beta   <=> no part of other (big) clones
                    if not tcrs:   # no bad incident tcrs
                        nonclonal_cells.append(cell)
               
        return nonclonal_cells
        
    #def refine_nonclonal(self, cl, ncl):
    #    intersect = set.intersection(set(cl), set(ncl))
    #    if len(intersect) > 0:
    #        #print(red('Warning: ') + '%d cells are in both sets and are excluded from analysis' % len(intersect))
    #        None
    #    ncl = list(set(ncl) - intersect)
#
    #    for cells in [cl, ncl]:
    #        if len(cells) == 0:
    #            print(red('Warning: ') + 'len(cells)=0. Returning')
    #            return (None, None)
    #        if len(cells) < 3:
    #            #print(red('Warning: ') + 'len(cells)=%d' % len(cells))
    #            None
    #        if not set(cells) <= set(self.batch.cells()):
    #            print(red('Warning: ') + 'the cells [%s] not in batch [%s]' %
    #                  (iter2str(set(cells)-set(self.batch.cells())), iter2str(self.batch.cells())))
    #    return ncl

    def get_potential_tcrs(self):
        tcrs = self.batch.tcrs
        cells = tcrs.cells()

        tcrs2cellcnt = Counter()
        tcr2cellcnt = Counter()
        for cell in cells:
            alphas, betas = tcrs.cell2tcrs(cell)
            for a in alphas:
                for b in betas:
                    tcrs2cellcnt[(a,b)] += 1.0 / len(alphas) / len(betas)

            for tcr in [*alphas, *betas]:
                tcr2cellcnt[tcr] += 1

        #nonclonal = []
        #for cell in cells:
        #    alphas, betas = tcrs.cell2tcrs(cell)
        #    al, bl = len(alphas), len(betas)
        #    if (al>0 or bl>0):
        #        nonclonal.append((cell, (alphas, betas)))            
        #        #nonclonal.append()
#
        #    #if len(alphas) == 1 and len(betas) == 1:
        #    #if al == 0 or (al == 1 and tcr2cellcnt[alphas[0]] == 1):
        #    #    if bl == 0 or (bl == 1 and tcr2cellcnt[betas[0]] == 1):
        #    #        if al == 0 or bl == 0 or tcrs2cellcnt[(alphas[0], betas[0])] <= 1.0:
        #    #            nonclonal.append((cell, (alphas, betas)))            

        # TODO
        #print(tcrs2cellcnt)
        clonal = []
        for tcr, cnt in tcrs2cellcnt.most_common():
            #if cnt > 1.0:
                clonal.append((cnt, tcr))

        return clonal

    def tcrs2cells(self, tcrs):
        cells_sets = [ set(self.batch.tcrs.select_cells(tcr)) for tcr in tcrs ]
        return list( set.intersection(*cells_sets) )

    def tcrs2clone(self, alpha, beta, logic=''):
        if logic == 'and':
            assert alpha != '' and beta != ''
            clonal_cells = self.tcrs2cells([alpha, beta])
            return clonal_cells
        
        if logic == 'ond':
            and_clonal_cells = self.tcrs2clone(alpha, beta, 'or')
            clonal_cells = []
            for cell in and_clonal_cells:
                A, B = self.batch.tcrs.cell2tcrs(cell)
                if len((set(A) | set(B)) - set([alpha, beta])) == 0:
                    clonal_cells.append(cell)
            return clonal_cells

        alpha_cells = self.tcrs2cells([alpha])
        beta_cells = self.tcrs2cells([beta])
        if alpha != '' and beta != '': 
            assert logic == 'or', 'logic should be "or", instead it is "{}"'.format(logic)
            if min(len(beta_cells), len(alpha_cells)) <= 0:
                return []
            clonal_cells = list(set(alpha_cells).union(set(beta_cells)))
        elif alpha != '':
            return alpha_cells
        else:
            assert beta != ''
            return beta_cells

        return clonal_cells

    def process(self, alpha, beta, logic=''):
        #print('Pesho: find_clones, process', alpha, beta, logic)
        
        clonal_cells = self.tcrs2clone(alpha, beta, logic=logic)
        #nonclonal_cells = self.refine_nonclonal(clonal_cells, nonclonal_cells)
        nonclonal_cells = self.get_nonclonal_cells(clonal_cells)
        #if not clonal_cells or not nonclonal_cells:
        #    return
        #
        #if min(len(nonclonal_cells), len(clonal_cells)) < 5:
        #    return

        h = Hypothesis(self.batch, nonclonal_cells=nonclonal_cells, clonal_cells=clonal_cells,
                       alpha=alpha, beta=beta, logic=logic, diffexpr=self.diffexpr) 
        #h.print_best_genes(out_fn=self.out_fn)
        #print(h.get_formatted_genes())
        self.hypotheses[ h.clone_tag() ] = h
        
        if self.batch.cell_part:
            noncl = '{:>3}'.format(len(nonclonal_cells))
            cl = '{:>3}'.format(len(clonal_cells))
            line = ' {} noncl, {} cl: {:>40} {} {:<30} --> {} genes'.format(
                blue(noncl), red(cl), alpha, green('{:3}'.format(logic)), beta, h.number_of_best_genes())
            print(line)

        
    def split_rule(self, onlytag):
        arr = onlytag.split()
        
        if len(arr) == 3:
            alpha, logic, beta = arr
            return (alpha, beta, logic)
        elif len(arr) == 1:
            alpha = arr[0]
            beta = ''
            if 'TRB' in alpha:
                alpha, beta = beta, alpha
                assert not 'TRA' in beta
            return (alpha, beta, '')
        assert 'Wrong # of parts: {}'.format(arr)

    def fill_clone_hypotheses(self, onlytag=''):
        color = self.batch.color if type(self.batch.color) is str else 'black'
        if self.batch.cell_part is not None:
            gene_num = len(self.batch.genes())
        else:
            gene_num = '?'
        print('{} ({} cells x {} genes)'.format(
              xkcd2console(self.batch.tag, color),
              blue(len(self.batch.cells())), blue(gene_num)))
        self.get_max_nonclone_size(out=True)
        
        if onlytag == 'auto':
            alphas, betas = set(), set()
            for cnt, tcr in self.get_potential_tcrs():
                alpha, beta = tcr
                for logic in ['and', 'or', 'ond']:
                    self.process(alpha, beta, logic=logic)

                alphas.add(alpha)
                betas.add(beta)

            # check agains clones defined by (a single alpha), and
            # check agains clones defined by (a single beta)
            for alpha in alphas:
                self.process(alpha, '')
            for beta in betas:
                self.process('', beta)

            best_hypo = self.best_hypo()
            if best_hypo:
                print()
                print('best hypo: {} out of {}'.format(best_hypo, blue(len(self.hypotheses))))
            #self.hypotheses.sort()
            
        if onlytag:
            self.process(*self.split_rule(onlytag))
            return
