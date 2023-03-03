import pandas as pd
import collections

from algo.find_clones import FindClones
from vis.quant_heatmap import QuantHeatmap
from utils import *

class GridCrossValidation:
    def __init__(self, *, clf, features, diffexprs, Cs, kernels, gammas):
        self.clf = clf
        self.features = features
        self.diffexprs = diffexprs
        self.Cs = Cs
        self.kernels = kernels
        self.gammas = gammas
        
        self.params2res = {}   # (sample_name, kernel, C)
       
    def _run_classifier(self, *, batch, kernel, c, gamma, diffexpr, feature):
        assert c > 0 and gamma >= 0
        
        clonotype = batch.onlytag# if batch.onlytag else 'NoClone'
        assert clonotype
        #hypo = batch.clones.best_hypo()
        hypo = FindClones(batch, onlytag=clonotype, diffexpr=diffexpr).best_hypo()
        res = self.clf(hypo, kernel, c)
        mean = res.scores.mean() if not res.scores is None else 0.0
        std2 = res.scores.std() * 2 if not res.scores is None else 0.0
        #if kernel == 'rbf':
        #    c = 'c={}, g={}'.format(c, gamma)
        key = (batch.tag, clonotype, diffexpr, feature, kernel, c, gamma, mean, std2)   ####### !!!! ######
        
        if key in self.params2res:
            self.params2res.pop(key, None)
            print("Warning: overwriting the classifier for ", key)
            #if self.lamp:
            #    print("Warning: overwriting the classifier for ", key)
            #    self.lamp = False
            assert not key in self.params2res
        self.params2res[key] = res
        
    def calc_batch_grid_cv(self, batch, feature):
        if not isinstance(self.Cs, collections.Iterable):
            self.Cs = [self.Cs]
        if isinstance(self.kernels, str):
            self.kernels = [self.kernels]
            
        #QuantHeatmap.draw(batch, feature)

        for diffexpr in self.diffexprs:
        #    print('diffexpr:', green(diffexpr))
            for kernel in self.kernels:
                for c in self.Cs:
                    _gammas = self.gammas if kernel == 'rbf' else [1]
                    for gamma in _gammas:
                        self._run_classifier(batch=batch, diffexpr=diffexpr, kernel=kernel,
                                             c=c, gamma=gamma, feature=feature)
    
    def results2html(self, out_fn=None):
        html = '<table width="100%">'
        if not self.params2res is None:
            data = list(self.params2res.keys())
            df = pd.DataFrame(data, columns=['batch', 'clonotype', 'diffexpr', 'features', 'kernel',
                                             'C', 'gamma', 'mean', 'std2'])
            #unique_features = df['clonotype'].unique()
            #for i, row in df.iterrows():
            #    print(row['a'])
            unique_features = { (row['batch'], row['clonotype']) for i, row in df.iterrows() }
            dfs = { f: pd.DataFrame for f in unique_features }
            for b, clonotype in dfs.keys():
                not_pivoted = df[:][ (df['batch'] == b) & (df['clonotype'] == clonotype) ]
                pivoted = not_pivoted.pivot_table(index=['C'],
                                                  columns=['batch', 'clonotype', 'diffexpr', 'kernel', 'features', 'gamma'],
                                                  values='mean')
                table = 'Table for <b>{}: {}</b> <br>'.format(b, clonotype)
                table += pivoted.to_html() + '<br>'
                row = '<td>{}</td><br>'.format(table)
                
                #for bi_fn in self.bipartite_files:
                #    assert (self.out_prefix / bi_fn).is_file()
                #    if b in str(bi_fn) and clonotype in str(bi_fn):
                #        row += '<td><img src="{}" height="500"></td>\n'.format(str(bi_fn))
                #        break
                        
                html += '<tr>{}</tr>'.format(row)
        html += '</table>'
        
        # write to file
        if out_fn:
            print('Writing cross-validation results to {}'.format(out_fn))
            html_file = open(out_fn, "w")
            html_file.write(html)
            html_file.close()
        
        return html