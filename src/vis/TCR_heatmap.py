import numpy as np
import matplotlib.pyplot as plt
from utils import *
import seaborn as sns

class TCRheatmap:
    @classmethod
    def create_heatmap(cls, tcrs):
        Vs = tcrs.get_all_VDJ('V')
        Js = tcrs.get_all_VDJ('J')
        
        #display(tcrs.tcrs)
        #display(Vs)
        
        v2ind = { v: i for (i, v) in enumerate(Vs) }
        j2ind = { j: i for (i, j) in enumerate(Js) }
        
        img = np.zeros((len(Vs), len(Js)), dtype=np.int)
        for ind, cols in tcrs.tcrs.iterrows():
            if cols['V'] and cols['J']:
                i = v2ind[ cols['V'] ]
                j = j2ind[ cols['J'] ]
                img[i, j] += 1
        return img, Vs, Js
                
    @classmethod
    def plot(cls, ax, fig, img, Vs, Js, interactive, out_fn, *, chain):
        # global
        if chain == 'A':
            title = 'alpha'
        elif chain == 'B':
            title = 'beta'
        else:
            assert False
        ax.set_title('{} chain, {} seqs'.format(title, img.sum().sum()), fontsize=10)
        ax.grid(False)
        
        # image
        cmap = plt.cm.summer_r
        cmap.set_under(color=xkcd2rgb('white'))  # 'cream'
        im = ax.imshow(img, cmap=cmap, vmin=1, interpolation='none')
        
        # ticks
        ax.set_yticks(range(len(Vs)))
        ax.set_yticklabels(Vs) #, rotation=-30) #, ha='center')
        ax.set_xticks(range(len(Js)))
        ax.set_xticklabels(Js, rotation=-90) #, ha='center')
        
        # numbers in cells
        for (j,i), label in np.ndenumerate(img):
            if label > 0:
                ax.text(i, j, label, ha='center', va='center')
            
        #fig.colorbar(im)
        
    @classmethod
    def draw(cls, batch, *, interactive, out_fn=None):
        plt.ioff()
        #plt.show()
        
        #sns.set()
        #g = sns.clustermap(img)
        #return
        
        fig, (axa, axb) = plt.subplots(1, 2)
        fig.suptitle('V-J heatmap for ' + batch.tag, fontsize=12, weight='bold')
        fig.tight_layout()
        
        for chain, ax in [('A', axa), ('B', axb)]:
            tcrs = batch.tcrs.subtcrs(cells=None, chain=chain)
            display(tcrs)
            img, Vs, Js = cls.create_heatmap(tcrs) 
            cls.plot(ax, fig, img, Vs, Js, interactive, out_fn, chain=chain)
        plot(plt, interactive, out_fn)
        