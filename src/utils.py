import pandas as pd
import random
#import warnings
from IPython.display import display
import collections
from libs import colortrans
import seaborn as sns
from pathlib import Path
#from enum import Enum
from matplotlib.colors import ColorConverter

from file_utils import *
from structs.genomics_convert import GenomicsConvert
from structs.vb_trbv_helper import VbTrbvHelper

### CONTSTANTS ###
#GENES = 0  # 0 <=> all
#MAX_P_VALUE = 0.9
MAX_P_VALUE = 0.9
MIN_CELLS_IN_CLONE = 5

def colorify(msg, ansi_color):
    # usage colorify('colortext', '#121212')
    if ansi_color.startswith('#'):
        short, rgb = colortrans.rgb2short(ansi_color)
        ansi_color = '\033[38;5;%sm' % short
    
    if not isinstance(msg, str) and isinstance(msg, collections.Iterable):
        return [ colorify(elem, ansi_color) for elem in msg ]
    else:
        return (ansi_color + str(msg) + "\x1b[0m")

def xkcd2console(msg, xkcd):
    return colorify(msg, xkcd2rgb(xkcd))

# returns a string '#rrggbb'
def xkcd2rgb(xkcd):
    c=sns.xkcd_rgb[xkcd]  
    return c

# returns RGBA
#def set_alpha(c, coef):
#    c = ColorConverter.to_rgba(c)
#    c = c[0], c[1], c[2], coef
#    return c

def add_more_alpha(c, coef):
    c = ColorConverter.to_rgba(c)
    c = (c[0], c[1], c[2], c[3]*coef)
    return c

def xkcd2rgba(xkcd):
    return tuple(ColorConverter.to_rgba(xkcd2rgb(xkcd)))

def to_3rgb_floats(color):
    if not color.startswith('#'):
        #if not color.startswith('xkcd:'):
        #    color = 'xkcd:' + color
       #mcd.XKCD_COLORS[color]
        #return color
        #return ColorConverter.to_rgb(color)
        return xkcd2rgb(color)
    return color

def blue(msg):
    return colorify(msg, '\x1b[34m')

def red(msg):
    return colorify(msg, '\x1b[31m')

def green(msg):
    return colorify(msg, '\x1b[32m')

def check_file(f, msg="", silent=False):
    if not Path(f).is_file():
        if not silent:
            print(red("NOT FOUND: ") + "%s %s" % (msg, blue(str(f))))
        return False
    return True

def to_str_list(d):
    return [ str(x) for x in d ]

def iter2str(c):
    return ', '.join([str(x) for x in c])

def print_df(df, full=False):
    if full:
        display(df)
    else:
        display(df.head(n=10))

def random_color(seed=None):
    if seed: random.seed(seed)
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())

def add_prefix_tag(df, tag):
    rename_cells = lambda cell: tag + '_' + cell
    df.rename(index=rename_cells, inplace=True)
    
#def add_prefix_tag_multiindex(df, tag):
#    rename_cells = lambda cell: tag + '-' + cell
#    cells = df.index._get_levels()[0].map(rename_cells)
#    df.index.set_levels(levels=cells, level=0, inplace=True)
    

def assert_no_NaNs(df, *, fail=True):
    if df.isnull().any(1).any():
        print('Cells with NaNs: %s' % blue(str(df.isnull().any(1).loc[lambda b: b == True].index.tolist())))
        display(df.head)
        if fail:
            raise Exception('NaNs in the GE matrix')
        
def assert_same_indexes(should_be, tested, tag=None, subset_ok=False):
    #err_msg = 'The index {%s} is not of type "object"/"string" but of %s.'
    #assert should_be.dtype_str == 'object', \
    #       err_msg % (iter2str(should_be), type(should_be[0]))
    #assert tested.dtype_str == 'object', \
    #       err_msg % (iter2str(tested), type(tested))
    
    should_be_set = set(to_str_list(should_be))
    tested_set = set(to_str_list(tested))
    
    if not tag: tag = ''
        
    if should_be_set != tested_set:
        if set(tested_set) <= set(should_be_set):
            if subset_ok:
                return True
            msg = ("%s covers only %s cell names out of the %s known."
                   % (blue(tag), blue(len(tested_set)),
                      blue(len(should_be_set))))
            return False
        else:
            unknown = tested_set - should_be_set
            msg = ("%s has %s unknown cell names: %s"
                   % (blue(tag), blue(len(unknown)), blue(iter2str(unknown))))
            #msg += (" Should be (len=%3d):  %s\n"
            #        "Instead of (len=%3d):  %s\n"
            #        % (len(should_be_set), iter2str(should_be_set),
            #           len(tested_set), iter2str(tested_set)))
            #print(red('WARNING: ') + msg)
            return False
    return True

def warning(msg):
    print('{}: {}'.format(red('Warning: '), msg))
        
def refresh_index(df):
    #df = pd.DataFrame({'a': [1,2,3], 'b': [5,8,10]}, index=['08','2','14'])
    df.index = pd.Index([int(x) for x in df.index])
    df.sort_index(inplace=True)
    df.index = pd.Index([str(x) for x in df.index])
    
def pval2str(pval):
    if pval < 0.001:
        return "%.2e" % pval
    else:
        return "%.2f%%" % (100.0*pval)
    
def pval2colstr(pval):
    s = pval2str(pval)
    if pval < 0.05:
        return green(s)
    else:
        return red(s)

dirs = Dirs('..')
gene_db = GenomicsConvert.from_file(dirs.inp)
vb_trbv_helper = VbTrbvHelper.init_from_file(dirs.inp / 'vb-table.tsv')

def tid2tname(tid):
    if not isinstance(tid, str) and isinstance(tid, collections.Iterable):
        return [ tid2tname(t) for t in tid ]
    else:
        return gene_db.tid2tname(tid)
    
def sorttable_html():
    with open(dirs['src/libs/sorttable.js'], 'r', encoding='utf-8', errors='ignore') as sorttable_file:
        sorttable  = '<html><body><script type="text/javascript">'
        sorttable += sorttable_file.read()
        #sorttable += '/* alabala */'
        sorttable += '</script></body></html>'
        return sorttable
    raise Exception("hm..")

def plot(plt, interactive, out_fn):
    # doens't have to be commented out
    #if interactive:
    #    plt.show()
    if out_fn:
        Path(out_fn).parent.mkdir(parents=True, exist_ok=True)
        #print("Figure saved in {}".format(blue(Path(out_fn).relative_to(dirs['out']))))
        print("Figure saved as {}".format(blue(Path(out_fn).name)))
        plt.savefig(str(out_fn), dpi=300)
    #plt.gcf().clear()
    plt.close()
    
def add_new_ax(fig):
    # http://matplotlib.1069221.n5.nabble.com/dynamically-add-subplots-to-figure-td23571.html
    # now later you get a new subplot; change the geometry of the existing
    for i in range(len(fig.axes)):
        fig.axes[i].change_geometry(n+1, 1, i+1)

    # add the new
    ax = fig.add_subplot(n+1, 1, n+1)
    return ax
    
from running import Run
run = Run(dirs)

