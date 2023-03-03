#from ipywidgets import interact, interactive, fixed, interact_manual

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib import collections  as mc
import matplotlib.lines as lines

import seaborn as sns; sns.set()
from .mpl_annote_finder import AnnoteFinder

import algos
from utils import *

class ScatterTCR:
    def __init__(self, batch, *, dim_reduct, perplexity=30, interactive=False, out_fn=None,
                 get_color, points2d=None):
        self.points2d = None
        if not points2d is None:
            self.points2d = points2d
        
        q = batch.quant
        if not q or len(q.cells()) < 2:
            warning('Not enough cells for ScatterTCR for batch ' + batch.tag)
            return

        #if genes is None:
        #    ge = q.ge
        #else:
        #    ge = q.ge[genes]
        
        ge = q.ge

        if ge.shape[1] < 2:
            print(red('Warning: Cannot draw ScatterTCR. Only {} (<2) genes for batch {}'.format(ge.shape[1], batch.tag)))
            return
            
        plt.ioff()
        # 2D
        if self.points2d is None:
            print("Processing %s cells with %s dimensions..." %
                  (blue(len(q.cells())), blue(len(q.transcripts()))))

            self.points2d, var = algos.getDimReduct(ge, dim_reduct, perplexity, 2)
        
        if not self.points2d is None:
            edges = batch.tcrs.getEdges()
            fig = plt.figure() #figsize=(20,10))
            ax2d = fig.add_subplot(111)

            self.scatter2D(ax2d, batch, fig, self.points2d, edges, get_color)
            self.add_annot(fig, ax2d, self.points2d, var, dim_reduct, ge, initOn=False)

            # 3D
            #points3d = pesho_algos.getDimReduct(q, dim_reduct, perplexity, 3)
            #ax3d = fig.add_subplot(122, projection='3d')
            #pesho_visual.scatter3D(ax3d, points3d, dim_reduct + ' over ' + filter_expr, edges)

            plot(plt, interactive, out_fn)
        else:
            warning('points2d is None for batch ' + batch.tag)


    def scatter2D(self, ax, batch, fig, points, edges, get_color):
        colormap = { 'TCRA': sns.xkcd_rgb["medium blue"],   #2c6fbb
                     'TCRB': sns.xkcd_rgb["greyish green"] }#82a67d

        for cell in points.index:
            c = to_3rgb_floats(get_color(cell))
            annot = get_color(cell, annot=True)   ## clonal/nonclonal/unknown
            b = batch.meta.cell2batch(cell)
            s = batch.meta.cell2sample(cell)
            label = b
            if not b.endswith(s):
                label += '-' + s
            if annot:
                label += '-' + annot
            
            points.loc[cell, 'color'] = c
            points.loc[cell, 'label'] = label

        # nodes
        for label in sorted(points['label'].unique()):
            subset = points[points.label == label]
            label += ' ({})'.format(len(subset.index))
            scatter = ax.scatter(subset['x'], subset['y'], c=subset['color'], label=label,
                                 s=60, edgecolors='face', alpha=0.7, zorder=2);
        ax.legend(shadow=True)

        # edges
        if len(edges) > 0:
            lines = []
            colors = []
            linewidths = []
            #title += '+TraCeR'
            random.shuffle(edges)
            for a, b, tcr in edges:
                assert(tcr.startswith('TRA') or tcr.startswith('TRB')), 'Bad tcr: %s' % tcr
                c = colormap['TCRA'] if tcr.startswith('TRA') else colormap['TCRB']
                linewidth = 0.6 if c == colormap['TCRA'] else 0.3
                if a in points['x'] and b in points['y']:
                    x0 = points.loc[a, 'x']
                    x1 = points.loc[b, 'x']
                    y0 = points.loc[a, 'y']
                    y1 = points.loc[b, 'y']
                    
                    lines.append([(x0, y0), (x1, y1)])
                    colors.append(c)
                    linewidths.append(linewidth)
                    
            lc = mc.LineCollection(lines, colors=colors, linewidths=linewidths, alpha=0.3, zorder=1)
            ax.add_collection(lc)

    def add_annot(self, fig, ax, points, var, dim_reduct, ge, initOn):
        title = '{} over {} cells with {} genes'.format(dim_reduct, len(points), len(ge.columns))
        if not var is None:
            perc = map(lambda x: str(int(x*100))+'%', var) 
            title += ' expl_var=[{}]'.format(','.join(perc))
        ax.set_title(title)
        ax.grid(False)
        ax.axis('off')
        
        # interactive annotations on click
        # labels
        # TODO: add TCR seq
        points['label'] = [cell for cell in points.index]
        #tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
        #mpld3.plugins.connect(fig, tooltip)
        
        # legend
        #red_dot, = plt.plot(z, "ro", markersize=15)
        #red_patch = mpatches.Patch(color='red', label='The red data')
        #ax.legend([red_patch, ], ["Red dot", ])

        af = AnnoteFinder(points['x'], points['y'], points['label'], ax=ax, initOn=initOn)
        fig.canvas.mpl_connect('button_press_event', af)

#def interactive_scatter_tcr_graph(batch):
#    filter_expr = ''
#    reduct_methods = ['PCA', 'DUMB', 'SNE', 'MDS']
#    interact(processAndDraw, batch=fixed(batch),
#             dimreduct=reduct_methods, filter_expr=filter_expr,
#             perplexity=30);
