import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib
from pandas.core.indexing import _maybe_numeric_slice, _non_reducing_slice
from IPython.display import HTML
from itertools import chain
from collections import defaultdict
import openpyxl

from utils import *

#idx = pd.IndexSlice

#class MidpointNormalize(colors.Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        colors.Normalize.__init__(self, vmin, vmax, clip)
#
#    def __call__(self, value, clip=None):
#        # I'm ignoring masked values and all kinds of edge cases to make a
#        # simple example...
#        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y))
    
class GeneIntersect:
    def __init__(self, directory, *, out_fn=None):
        self.dir = Path(directory)
        Path(directory).exists()
        self.files = list(Path(directory).glob('*.csv'))  #'genes_*.csv'
        
        if self.files:
            print('Merging all gene files from {}:'.format(blue(self.dir)))
            print(', '.join([blue(fn.name) for fn in self.files]))
            #self.genes_vs_samples_table(Path(directory) / "genes.html")
        else:
            warning('No gene files in {}'.format(blue(self.dir)))

        self.fn2sample = {}
        for fn in self.files:
            sample = self.read_sample(fn)
            if len(sample.index) > 0:
                self.fn2sample[fn.stem] = sample
            else:
                self.fn2sample[fn.stem] = pd.DataFrame(columns=['symbol', 'pvalue', 'logFC'])
                print(red('Warning: Zero genes in {}'.format(str(Path(fn).name))))

    def read_sample(self, fn):
        #print('Loading genes from {}...'.format(blue(str(Path(fn).name))))
        # names=['p_adj', 'transcript_name', 'target_id'], 
        df = pd.read_csv(fn)  #, index_col='transcript_id') #, index_col='target_id')
        if 'Unnamed: 0' in df:
            df.drop('Unnamed: 0', axis=1, inplace=True)
        return df

    def get_fn2sample(self, fn_should_include=''):
        if fn_should_include:
            fn2sample = { fn: sample for fn, sample in self.fn2sample.items() if (fn_should_include in fn) }
        else:
            fn2sample = self.fn2sample
        return fn2sample
        
    def dump_gene_tables(self, out_fn_tmpl=None, *, fn2sample):
        Path(out_fn_tmpl).parent.mkdir(parents=True, exist_ok=True)
        if not self.fn2sample:
            warning('Zero gene files.')
            return

        for column in [ "logFC" ]: #, "FDR", "pvalue" ]:
            self.dump_table_from_column(
                    column, out_fn=str(out_fn_tmpl).format(column), out_bool_fn=str(out_fn_tmpl).format("bool"),
                    fn2sample=fn2sample)

    def dump_table_from_column(self, column, *, out_fn, out_bool_fn, fn2sample):
        # https://stackoverflow.com/questions/38931566/pandas-style-background-gradient-both-rows-and-columns
        def background_gradient(s, m, M, cmap, low=0, high=0):
            rng = M - m
            norm = colors.Normalize(m - (rng * low),
                                    M + (rng * high))
            normed = norm(s.values)
            c = [colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
            return ['background-color: %s' % color for color in c]

        def sample2serie(sample, name):
            for col in [ column, 'symbol' ]:
                if col not in sample:
                    return None
            return pd.Series(data=list(sample[column]), index=sample['symbol'], name=name)
            #s.index.name = name

        #print(self.fn2sample.values())
        all_samples = [ sample2serie(sample, fn) for fn, sample in fn2sample.items() if sample is not None ]


        # get gene descriptions
        descr = dict()
        for _, sample in fn2sample.items():
            if 'description' in sample:
                for symbol, d in sample[ ['symbol', 'description'] ].values:
                    descr[symbol] = d

        if not all_samples or all_samples[0] is None:
            warning('Missing {} column'.format(column))
            return

        #for serie in all_samples:
        #    #print(serie)
        #    dupl = serie.index.duplicated()
        #    if sum(dupl):
        #        warning('Duplicated symbols.')
        #        print(serie.index.name)
        #        print(serie[dupl])
        #        assert(False)

        if all_samples:
            df = pd.concat(all_samples, axis=1, sort=True)
        else:
            Warning('No samples to concatenate for ', fn2sample)
            df = pd.DataFrame()
        df.index.name = 'symbol'
        df.sort_index(axis='columns', inplace=True)
        sample_columns = df.columns

        df = df.assign(
            #rank = lambda d: len(d.columns) - d.isnull().sum(axis='columns'),
            count     = lambda d: d.notna().sum(axis='columns'),
            count_pos = lambda d: d[d>0].notna().sum(axis='columns'),
            count_neg = lambda d: d[d<0].notna().sum(axis='columns'),
            count_diff= lambda d: d[d>0].notna().sum(axis='columns') - d[d<0].notna().sum(axis='columns'),
            sum_pos   = lambda d: d[d>0].sum(axis='columns'),
            sum_neg   = lambda d: d[d<0].sum(axis='columns'),
            sum       = lambda d: d.sum(axis='columns'),
        )
        #df = df[ df['count'] > 1 ]
        #df['only_pos?'] = (df['count_pos'] > 0) & (df['count_neg'] == 0)
        #df['only_neg?'] = (df['count_pos'] == 0) & (df['count_neg'] > 0)
        df['same_sign?'] = (df['count_pos'] == 0) | (df['count_neg'] == 0)
        df['description'] = pd.Series(df.index, index=df.index).map(descr)

        df.sort_values(by=['count_diff', 'sum'], ascending=[False, False], inplace=True)
        #df = df.round(2) 
        #print(df.iloc[:5, :3])
        df.to_csv(out_fn, float_format='%.3g')
        out_xlsx = Path(out_fn).with_suffix('.xlsx')
        print('Writing to {}'.format(blue(out_xlsx.name)))
        #print(df.head())

        df_styler = df.style.apply(background_gradient,
                           #cmap='RdBu_r',  # TODO: colorscheme `bwr'
                           cmap='bwr',  # TODO: colorscheme `bwr'
                           m=-7, #df.min().min(),
                           M=+7, #df.max().max(),
                           low=0,
                           high=0,
                           subset=sample_columns).highlight_null(null_color='white')

        #print(df)

        df_styler.to_excel(out_xlsx, engine='openpyxl')

        # boolean matrix
        #if not Path(out_bool_fn).exists():
        #    df_isin = df.copy()
        #    df_isin[ ~df_isin.isnull() ] = 1
        #    df_isin[ df_isin.isnull() ] = 0
        #    df_isin.drop(columns=['count'], inplace=True)
        #    #df_isin.fillna(0)
        #    df_isin.astype(int).to_csv(out_bool_fn)
    
    def dump_stats(self, out_fn, *, fn2sample):
        Path(out_fn).parent.mkdir(parents=True, exist_ok=True)
        fn2len = { fn: len(sample) for fn, sample in fn2sample.items() if sample is not None }
        stats = pd.DataFrame(list(fn2len.items()), columns=['DE experiment', 'genes'])
        #stats['DE genes'] = 
        stats.set_index('DE experiment')
        stats.sort_index(inplace=True)
        
        print('Writing to {}'.format(blue(Path(out_fn).name)))
        stats.to_csv(out_fn, index=False)  # astype('int64')
            
    def dump_experiment_comparison_table(self, out_fn_tmpl, *, fn2sample_X, fn2sample_Y):
        Path(out_fn_tmpl).parent.mkdir(parents=True, exist_ok=True)

        def print_csv(df, tag):
            out_fn = str(out_fn_tmpl).format(tag)
            print('Writing to {}'.format(blue(Path(out_fn).name)))
            if df.empty:
                warning('Empty dataframe printed to {}.'.format(out_fn))
                df = pd.DataFrame(['-1'])
            df.to_csv(out_fn)

        def print_xlsx(df, tag):
            out_xlsx = Path(str(out_fn_tmpl).format(tag)).with_suffix('.xlsx')
            print('Writing to {}'.format(blue(out_xlsx.name)))
            if df.data.empty:
                warning('Empty dataframe printed to {}.'.format(out_xlsx))
                df.data = pd.DataFrame(['-1'])
            df.to_excel(out_xlsx, engine='openpyxl')

        # https://stackoverflow.com/questions/38931566/pandas-style-background-gradient-both-rows-and-columns
        def background_gradient2(s, m, M, cmap, low=0, high=0):
            rng = M - m
            norm = colors.Normalize(m - (rng * low),
                                    M + (rng * high))
            vals = [ int(num_str.split(':')[0]) for num_str in s.values ]
            normed = norm(vals)
            c = [colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
            return ['background-color: %s' % color for color in c]

        def heat(df, color):
            cm = sns.light_palette(color, as_cmap=True) #, reverse=True)
            return df.style.apply(background_gradient2,
                           cmap=cm, #'BuGn', 
                           m=0, #df.min().min(),
                           M=20, #df.max().max(),
                           low=0,
                           high=0) 

        def genes2str(genes):
            return '{}: {}'.format(len(genes), ', '.join(genes))

        def fn2sample_2_fn2genes(fn2sample):
            return { fn: (sample['symbol'], sample['logFC']) for fn, sample in fn2sample.items() if sample is not None }

        tag2genes_X = fn2sample_2_fn2genes(fn2sample_X)
        tag2genes_Y = fn2sample_2_fn2genes(fn2sample_Y)

        sect_genes_df_both = None
        for upOrDown, sign, color in [ ('up', +1, 'red'), ('down', -1, 'blue') ]:
            sect_genes_df = pd.DataFrame().astype('object')
            #sect_num_df, sect_genes_df = pd.DataFrame(), pd.DataFrame()
            for tag1, genes1 in tag2genes_X.items():
                for tag2, genes2 in tag2genes_Y.items():
                    d = defaultdict(list)
                    for symbol, logFC in list(zip(*genes1)) + list(zip(*genes2)):
                        d[symbol].append(logFC)
                    #common_genes = [ symbol for symbol, logFCs in d.items() if len(logFCs) == 2 and logFCs[0]*logFCs[1] > 0 ]
                    common_genes = sorted([ symbol for symbol, logFCs in d.items() if len(logFCs) == 2 and sign*logFCs[0] >= 0 and sign*logFCs[1] >= 0 ])
                    #common_genes = sorted(_)

                    #common_genes = sorted(list(symbols1 & symbols2))
                    #sect_num_df.loc[tag1, tag2] = len(common_genes)
                    #if tag1 != tag2:
                    #print(tag1, tag2, common_genes)
                    sect_genes_df.at[tag1, tag2] = 'bla'
                    sect_genes_df = sect_genes_df.astype('object')
                    sect_genes_df.at[tag1, tag2] = common_genes

            #stats = pd.DataFrame(list(fn2len.items()), columns=['DE experiment', 'genes'])
            #stats.set_index('DE experiment')

            for axis in [0, 1]:
                #sect_num_df.sort_index(axis=axis, inplace=True)
                sect_genes_df.sort_index(axis=axis, inplace=True)

            # remove half
            #for i in range(sect_genes_df.shape[0]):
            #    for j in range(sect_genes_df.shape[1]):
            #        if i < j:
            #            sect_genes_df.iat[i, j] = []

            if sect_genes_df_both is None:
                sect_genes_df_both = sect_genes_df
            else:
                sect_genes_df_both += sect_genes_df
                sect_genes_df_both = sect_genes_df_both.applymap(lambda x: sorted(list(set(x))))

            print_csv(sect_genes_df.applymap(len).astype('int64'), '{}_num'.format(upOrDown))
            print_csv(sect_genes_df.applymap(', '.join), '{}_genes'.format(upOrDown))
            #print_xlsx(heat(sect_genes_df.applymap(len).astype('int64'), color), '{}_num'.format(upOrDown))
            print_xlsx(heat(sect_genes_df.applymap(genes2str), color), '{}_genes'.format(upOrDown))

        #print_csv(sect_genes_df_both.applymap(len).astype('int64'), 'upAndDown_num')
        #print_csv(sect_genes_df_both.applymap(', '.join), 'upAndDown_genes')

#        sect_genes_df_both = sect_genes_df_both.style.background_gradient(cmap=cm)

        #df = sect_genes_df_both.applymap(genes2str)
        print_xlsx(heat(sect_genes_df_both.applymap(genes2str), 'purple'), 'all_genes')

    #def two_lists_to_table(self, hristo_fn, jechko_fn):
    #    hristo = self.read_sample(hristo_fn)
    #    jechko = self.read_sample(jechko_fn)

    #    if hristo is None or jechko is None:
    #        return None

    #    res = pd.merge(hristo, jechko, suffixes=('_'+hristo_fn.name.split('_')[:10], '_'+jechko_fn.name),
    #                   left_index=True, right_index=True, how='inner')
    #    #print(res)
    #    #if not res.empty:
    #    #    print('res:', res)
    #    #    intersect1(hristo, jechko)
    #    return res

    #    try:
    #        hristo_genes = hristo.index.tolist()
    #        jechko_genes = jechko.index.tolist()
    #    except:
    #        return []

    #    return list(set(hristo_genes) & set(jechko_genes))
        
    #@staticmethod
    #def color_pvals(val):
    #    if val > 0.05:
    #        return 'background-color: white'
    #    else:
    #        return ''
    #    #color = xkcd2rgb('light green') if val != "" and float(val) < 0.05 else 'white'
    #    #return 'background-color: %s' % color
    #
    #@staticmethod
    #def color_log2fc_ratio(val):
    #    color = xkcd2rgb('light blue') if val != "" and abs(float(val)) >= 1.0 else 'white'
    #    return 'background-color: %s' % color
    #
    #@staticmethod
    #def color_nan(val):
    #    if pd.isnull(val):
    #        return 'color: white'
    #    else:
    #        return 'color: black'
    #
    #@staticmethod
    #def my_background_gradient(style, cmap='PuBu', low=0, high=0, axis=0, subset=None):
    #    subset = _maybe_numeric_slice(style.data, subset)
    #    subset = _non_reducing_slice(subset)
    #    style.apply(GeneIntersect._my_background_gradient, cmap=cmap, subset=subset,
    #                axis=axis, low=low, high=high)
    #    return style
    #    
    #@staticmethod
    #def _my_background_gradient(s, cmap='PuBu', low=0, mid=0, high=0):
    #    #rng = s.max() - s.min()
    #    # extend lower / upper bounds, compresses color range
    #    #norm = colors.Normalize(s.min() - (rng * low),
    #    #                        s.max() + (rng * high))
    #    ## matplotlib modifies inplace?
    #    ## https://github.com/matplotlib/matplotlib/issues/5427
    #    norm = MidpointNormalize(midpoint=mid, vmin=low, vmax=high)
    #    normed = norm(s.values)
    #    c = [colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
    #    return ['background-color: %s' % color for color in c]

    #@staticmethod
    #def set_style(style, p_cols):
    #    cm = sns.blend_palette(["green", "light green"], input='xkcd', as_cmap=True)
    #    div_cm = sns.diverging_palette(10, 240, n=9, as_cmap=True)
    #    #cm = sns.diverging_palette(10, 10, n=9, as_cmap=True)
    #    
    #    style.background_gradient(cmap=cm, low=.0, high=0.1, subset=(idx[:], idx[:, p_cols]))
    #    #GeneIntersect.my_background_gradient(style, cmap=cm, low=.0, high=0.1,
    #    #                                     subset=(idx[:], idx[:, p_cols]))
    #    style.applymap(GeneIntersect.color_pvals, subset=(idx[:], idx[:, p_cols]))
    #    GeneIntersect.my_background_gradient(style, cmap=div_cm, low=-5.0, high=5.0,
    #                              subset=(idx[:], idx[:, ['log2fc_ratio', 'log2fc_diff', 'avg_logFC']]))
    #    #TODO: uncommnet
    #    #style.format('{:.3}', subset=(idx[:], idx[:, [*p_cols, 'log2fc_ratio', 'log2fc_diff']]))
    #    
    #    style.highlight_null(null_color='white')
    #    style.applymap(GeneIntersect.color_nan)
    #    
    #    #def magnify():
    #    #    return [dict(selector="th",
    #    #                 props=[("font-size", "4pt")]),
    #    #            dict(selector="td",
    #    #                 props=[('padding', "0em 0em")]),
    #    #            dict(selector="th:hover",
    #    #                 props=[("font-size", "12pt")]),
    #    #            dict(selector="tr:hover td:hover",
    #    #                 props=[('max-width', '200px'),
    #    #                        ('font-size', '12pt')])
    #    #]
    #    #style.set_properties(**{'table-layout': 'fixed', 'width': '10px', 'font-size': '1pt'})\
    #    #    .set_precision(2)\
    #    #    .set_table_styles(magnify())
    #
    #    return style

    #def genes_vs_samples_table_old(self, out_fn=None):
    #    #df = pd.concat([sample['p_adj'].rename(fn.stem.split('_')[0])
    #    #                for fn, sample in self.fn2sample.items() if sample is not None], axis=1)
    #    if not self.fn2sample:
    #        print(red('Warning: Zero gene files.'))
    #        return
    #    
    #    all_samples = []
    #    for fn, sample in self.fn2sample.items():
    #        if sample is not None:
    #            s = sample.copy()
    #            if 'transcript_name' in s:
    #                del s['traescript_name']
    #            good_p_adj = False
    #            for p_adj_cand in ['p_adj(wilcox)', 'p_adj(ks)', 'p_val_adj']:
    #                if p_adj_cand in s:
    #                    good_p_adj = s[p_adj_cand] < 0.9
    #            good_fc = False
    #            if 'log2fc_diff' in s and 'log2fc_ratio' in s:
    #                good_fc = abs(s.log2fc_diff) > 1.0 & abs(s.log2fc_ratio) > 1.0
    #            elif 'avg_logFC' in s:
    #                good_fc = abs(s.avg_logFC) > 1.0
    #            s = s[ good_fc & good_p_adj ]
    #            col_name = fn.stem #.split('_')[1]
    #            s.columns = [ [col_name]*len(s.columns), s.columns ]
    #            all_samples.append(s)
    #            
    #    assert all_samples
    #        
    #    df = pd.concat(all_samples, axis=1)
    #    
    #    df.sort_index(axis='columns', inplace=True)
    #    p_cols = [ col for sample, col in df.columns if col.startswith('p_adj') or col.startswith('p_val_adj') ]
    #    #list(set(df.columns) - set(['transcript_name', 'rank', 'p_val5']))
    #    df = df.assign(
    #        #rank = lambda d: len(d.columns) - d.isnull().sum(axis='columns'),
    #        rank   = lambda d: len(d.columns) - d.isnull().sum(axis='columns'),
    #        p_max  = lambda d: d.loc[:, idx[:, p_cols]].max(axis='columns'),
    #        p_val5 = lambda d: (d.loc[:, idx[:, p_cols]] < 0.05).sum(axis='columns').astype(int),
    #    )
    #    df.insert(0, 'transcript_name', list(map(lambda tid: tid2tname(tid), df.index)))
    #    df.sort_values(by=['p_val5', 'p_max'], ascending=[False, True], inplace=True)
    #    
    #    if len(df) == 0:
    #        Warning("Empty table!")
    #        return
    #    
    #    #cm = sns.light_palette("green", as_cmap=True, reverse=True)
    #    #df = df.head().style.background_gradient(cmap=cm)
    #    
    #    #df_toshow = df.head(n=100)
    #    #df = df[:2000]
    #    #df = df.fillna('')
    #    #df = df.iloc[:100]
    #    df.sort_index(ascending=True, axis='columns', inplace=True)
    #    df_styler = GeneIntersect.set_style(df.style, p_cols)
    #    
    #    df_html_small = df[:20].style
    #    df_styler_small = GeneIntersect.set_style(df_html_small, p_cols)
    #    #display(HTML(df_html_small))
    #    
    #    df_html = df_styler.render()  #escape=False)
    #    
    #    #if not out_fn:
    #    #    out_fn = str(self.dir.parent / 'genes_{}.html'.format(self.dir.parent))
    #    
    #    if not out_fn:
    #        out_fn = self.dir/'genes.html'
    #        
    #    html = sorttable_html() + df_html
    #    html = html.replace('<table', '<table class="sortable"')
    #    #html = html.replace('class="dataframe"', 'class="sortable"')

    #    print('Writing genes_vs_samples table to {}'.format(blue(Path(out_fn).name)))
    #    html_file = open(out_fn, "w")
    #    html_file.write(html)
    #    html_file.close()

    #    if len(df) > 0:
    #        out_xlsx = str(Path(out_fn).with_suffix('.xlsx'))
    #        Path(out_xlsx).parent.mkdir(parents=True, exist_ok=True)
    #        print('Writing genes_vs_samples table to {}'.format(blue(Path(out_xlsx).name)))
    #        df_styler.to_excel(out_xlsx, engine='openpyxl')
    
    #def samples_vs_samples_table(self):
    #    for hristo_fn in self.files:
    #        for jechko_fn in self.files:
    #            if hristo_fn < jechko_fn:
    #                isect = self.two_lists_to_table(hristo_fn, jechko_fn)
    #                if isect is not None and not isect.empty:
    #                    #names = tid2tname(isect)

    #                    print('{} & {}: {} genes:'.format(blue(Path(hristo_fn).name), green(Path(jechko_fn).name), len(isect)))
    #                    display(isect)
    #                    print()
    #                    #print('\n   '.join(green(names)))
