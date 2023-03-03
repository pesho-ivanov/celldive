import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

from utils import vb_trbv_helper
from utils import *

class VbTrbv:
    def __init__(self, *, tag, vb2perc):
        # TODO
        vb2perc.loc['other'] = 100.0 - vb2perc['perc'].sum()
        
        self.tag = tag
        self.vb2perc = vb2perc
        '''Returns series vb -> perc'''
        
    @classmethod
    def from_file(cls, fn):
        if not fn.is_file() :
            Warning("Vb file %s does not exist!" % fn)
            return
        column_types = {'vb': str, 'perc': float}
        vb2perc = pd.read_csv(fn, delim_whitespace=True, comment='#',
                              index_col=0,
                              dtype=column_types)
        return cls(tag=fn, vb2perc=vb2perc)
    
    @classmethod
    def from_tcrs(cls, tcrs, tag):
        if not tcrs:
            Warning('No TCRs')
            return
        if len(tcrs.cells()) == 0:
            Warning('No cells in {}'.format(tag))
            return
        vb2perc = pd.DataFrame(index=vb_trbv_helper.vbs(), columns={'perc': float})

        def get_v(betas):
            return { beta['v'] for beta in betas }
        
        for vb in vb_trbv_helper.vbs():
            trbvs = vb_trbv_helper.vb2trbvs(vb)
            #subcells = tcrs.tcrs[tcrs.tcrs['V'].isin(trbvs)].index.get_level_values('cell').unique()
            subcells = [ cell for cell in tcrs.cells() if get_v(tcrs.cell2tcrs(cell, 'beta')) & set(trbvs) ]
            val = 100.0 * len(subcells) / len(tcrs.cells())
            #vb2perc.set_value(vb, 'perc', val)
            vb2perc.at[vb, 'perc'] = val
        return cls(tag=tag, vb2perc=vb2perc)
    
    @staticmethod
    def draw_all(batch, *, out_fn, interactive=False, ax=None):
        vb_trbvs = batch.vbs
        if vb_trbvs is None:
            Warning("No Vb not loaded.")
            return
        all_df = None
        for tag, vb_trbv in vb_trbvs.items():
            df = vb_trbv.vb2perc.rename(columns={'perc': tag})
            if all_df is None:
                all_df = df
            else:
                all_df = all_df.merge(df, left_index=True, right_index=True, how='outer')
        
        if not ax:
            plt.ioff()
            fig, ax = plt.subplots(1, 1)
            fig.suptitle(batch.tag, fontsize=12, weight='bold')
            #fig = plt.figure() #figsize=(20,10))
            #ax = fig.add_subplot(111) 
        #print('Vb statistics plotted to {}'.format(blue(out_fn.name)))
        #plt.close()
        all_df.plot(kind='bar', ax=ax)
        plot(plt, interactive, out_fn)
