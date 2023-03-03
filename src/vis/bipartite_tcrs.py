from collections import Counter, defaultdict, OrderedDict
import numpy as np
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
#from statistics import mode
from scipy.stats import mode

#from .mpl_annote_finder import AnnoteFinder
from structs.clone_struct import SampleColor, CellPart, CellColor
import math
from utils import *

# TODO: remove time
import time
tt = {}

class BipartiteTCR:
    _source = 'no β'
    _sink = 'no α'
    _rainbow_colors = ['red', 'yellow', 'green', 'cyan', 'blue']
    
    black = add_more_alpha('black', 0.7)

    @classmethod
    def _construct_graph(cls, batch, use_cellcolors):
        '''Returns a multigraph'''
        G = nx.MultiGraph()

        tcrs = batch.tcrs
        cells = tcrs.cells()
        tcrs_A = tcrs.chain2tcrs('alpha')
        tcrs_B = tcrs.chain2tcrs('beta')

        # Nodes
        G.add_nodes_from(tcrs_A, side='left')
        G.add_nodes_from(tcrs_B, side='right')
        G.add_nodes_from([cls._source, cls._sink], side='no')
            
        # background edges from No beta to No alpha 
        for cell in set(batch.cells()) - set(cells):
            #c_xkcd = batch.get_color(cell, coloring=SampleColor.CLONAL)
            if batch.cell_part:
                part = batch.cell_part.get_cell_part(cell)
            else:
                raise Exception("No cell partitioning.")
                part = 'UNKNOWN'
            G.add_edge(cls._source, cls._sink, 
                       color=xkcd2rgb(CellColor[CellPart.NONRECONSTRUCTED]),
                       cell=cell, clonal=part)
            
        def bipartite_cell_color(cell):
            if batch.cell_part is not None:
                clonal_t = batch.cell_part.get_cell_part(cell)
            else:
                clonal_t = 'UNKNOWN'
            is_certain = True
            if use_cellcolors == SampleColor.CLONAL:
                color = xkcd2rgba(batch.get_color(cell, coloring=SampleColor.CLONAL))
                #if lena > 1 or lenb > 1:
                #    if True or clonal_t == 'CLONAL':
                #        is_certain = False
                #        color = add_more_alpha(color, 0.5)
            elif use_cellcolors == SampleColor.SAMPLE:
                color = batch.meta.get_color_rgb(cell)
            else:
                if use_cellcolors != SampleColor.RAINBOW:
                    raise Exception('{} instead'.format(use_cellcolors))
                if lena <= 1 and lenb <= 1:
                    color = cls.black
                else:
                    color = random_color()
            color = add_more_alpha(color, 0.7)
            return is_certain, clonal_t, color

        # main edges
        tcr2cellcnt = Counter()
        clones = defaultdict(list)
        for cell in cells:
            alphas = tcrs.cell2seqs(cell, chains=['alpha'])
            betas = tcrs.cell2seqs(cell, chains=['beta'])
            lena, lenb = len(alphas), len(betas)
            
            for tcr in [*alphas, *betas]:
                tcr2cellcnt[tcr] += 1
            
            if lena > 0 and lenb > 0:
                edges = [ (a, b) for b in betas for a in alphas]
            elif lena > 0:
                edges = [ (cls._source, a) for a in alphas ]
                tcr2cellcnt[cls._source] += 1
            else:
                edges = [ (b, cls._sink) for b in betas ]
                tcr2cellcnt[cls._sink] += 1
            
            is_certain, clonal_t, color = bipartite_cell_color(cell)
            if CellPart.FILTERED != batch.cell_part.get_cell_part(cell):
                G.add_edges_from(edges, color=color, cell=cell)
                clones[clonal_t].append((color, is_certain))
            
        return G, tcr2cellcnt, clones
                
    @classmethod
    def _get_node_pos(cls, G):
        # positions
        left =  [ node for node in G.nodes() if G.node[node]['side'] == 'left' ]
        right = [ node for node in G.nodes() if G.node[node]['side'] == 'right' ]

        left.sort()
        right.sort()
        
        w = 10
        h = 20
        d = 0.2
        
        max_len = max(len(left), len(right))
        while (max_len * 2 * d  > h):
            h *= 2
            w *= 2

        pos = {}  # node -> [x,y]
        pos[cls._source] = [0, h/2.0]
        pos[cls._sink] = [w, h/2.0]
        for side, x in [(left, 1*w/3.0), (right, 2*w/3.0)]:
            for i, node in enumerate(side):
                y = h * (i+1.0) / (len(side)+1)
                pos[node] = [x, y]
                
        if len(G.nodes()) != len(pos):
            raise Exception("G.nodes() != len(pos)")
        return pos, w, h, d
    
    @classmethod
    def rainbow_edge(cls, ax, p, q, colors=_rainbow_colors, spacing=0.1):
        n = 14
        base = np.linspace(0.1, np.pi-0.1, n)
        shape = np.sin(base)
        x = np.linspace(p[0], q[0], n)
        y = np.linspace(p[1], q[1], n)
        
        delta = len(colors) * spacing / 2.0
        lines = []
        # TODO: group by color to speed up 3-5 times
        for i in range(len(colors)+1):
            newline_data = y - 0.00002*spacing*i ## ?!?!
            coef = len(colors)/2.0 - i
            newline_data += coef * spacing * shape
            lines.append(newline_data)
        for i, c in enumerate(colors):
            ax.fill_between(x, lines[i], lines[i+1], facecolor=c)

    @classmethod
    def _plot_bipartite(cls, ax, G, pos, w, h, d, alpha, beta, tcr2cellcnt):
        # TODO: colors
        # http://matplotlib.org/users/transforms_tutorial.html#using-offset-transforms-to-create-a-shadow-effect
        # http://stackoverflow.com/questions/42165631/in-matplotlib-how-can-i-plot-a-multi-colored-line-like-a-rainbow/42190453#42190453
        ax.grid(False)
        ax.axis('equal')
        ax.axis('off')
        #for text in texts:
        #    text.set_color('grey')
    
        # draw nodes
        for node in G.nodes():
            x, y = pos[node]
            #print('node: ', node)
            if node == alpha or node == beta:
                c = 'red'  
            elif node.startswith('NONPRODUCTIVE_A'):
                #ax.text(x, y, 'NP', c='black', ha='right', va='center', size=9) # TODO: brakes
                c = 'pink'
            elif node.startswith('NONPRODUCTIVE_B'):
                #ax.text(x, y, 'NP')  # TODO: brakes
                #ax.text(x, y, 'NP', c='black', ha='left', va='center', size=9)  # TODO: brakes
                c = 'pink'
            else:
                c = cls.black
            circle = Circle((x, y), d, color=c)
            ax.add_patch(circle)
            #if node == beta:
            #    ax.text(x, y, 'clone', c='r', ha='left', va='center', size=9)
            #if node == alpha:
            #    ax.text(x, y, 'clone', c='r', ha='right', va='center', size=9)

            #cnt = tcr2cellcnt[node]
            #if cnt > 1:
            #    t = ax.text(x, y, str(cnt), color='w', ha='center', va='center', weight='bold')
            #    t.set_bbox(dict(facecolor='black', alpha=0.2))
                
        # draw edges
        def cmp(edge):
            a, b = edge
            if a == cls._source and b == cls._sink:
                return -10000
            return min(tcr2cellcnt[a], tcr2cellcnt[b])  # bigger clones on top
        edgelist = sorted(G.edges(), key=cmp)
        #edge_cnt = Counter(edgelist)
        
        start = time.time()
        for a, b in edgelist:
            edges=G[a][b]
            edgecolors = [ e['color'] for e in edges.values() ]
            edgecolors = sorted(edgecolors)
            cls.rainbow_edge(ax, pos[a], pos[b], colors=edgecolors)
            # writing number of edges in circles
            #if len(edges) > 1:
            #    x, y = [ sum(p)/2 for p in zip(pos[a], pos[b]) ]
            #    if a == cls._source and b == cls._sink:
            #        x = pos[a][0] + (pos[b][0]-pos[a][0]) / 6.0 # move the non-reconstructed number of cells to the left
            #    #color = mode(edgecolors).mode[0]
            #    color = max(set(edgecolors), key=edgecolors.count)
            #    t = ax.text(x, y, str(len(edges)),
            #                color='w', ha='center', va='center', size=9, 
            #                bbox=dict(boxstyle="circle", fc=color, ls=None)) #, weight='bold')

        #ax.scatter = ax.scatter(points['x'], points['y'], c=points['color'], s=60, edgecolors='face', alpha=0.7, zorder=2);
        
        tt['plot_bipartite_edges'] = time.time() - start
        start = time.time()
        
        # column annotations
        x, y = pos[cls._source]
        size = 9
        c = 'grey'
        ax.text(x, y*1.05, r'no $\beta$', ha='right', va='center', size=size, color=c)
        x, y = pos[cls._sink]
        ax.text(x, y*1.05, r'no $\alpha$', ha='left', va='center', size=size, color=c)
        
        ax.text(1*w/3.0, 0, r'$\alpha$' '\nchains', ha='center', va='center', size=size, color=c)
        ax.text(2*w/3.0, 0, r'$\beta$' '\nchains', ha='center', va='center', size=size, color=c)
            
    @classmethod
    def _add_annot(cls, fig, ax, G, tcr2cellcnt, pos):
        # interactive annotations on click
        X = []
        Y = []
        annotes = []
        someOn = []
        for tcr, coord in pos.items():
            x, y = coord
            if 'TRA' in tcr: x *= 0.95
            if 'TRB' in tcr: x *= 1.03
            if tcr == cls._source: x *= 0.85
            if tcr == cls._sink: x *= 1.05
            X.append(x)
            Y.append(y)
            annotes.append(tcr)
            someOn.append(tcr2cellcnt[tcr] > 1 or tcr == cls._source or tcr == cls._sink)
            
        af = AnnoteFinder(X, Y, annotes, ax=ax, someOn=someOn)
        fig.canvas.mpl_connect('button_press_event', af)
        
    @classmethod
    def _plot_clonality_piechart(cls, clones, ax):
        def make_pie(ax, d):
            def make_autopct(values):
                def my_autopct(pct):
                    total = sum(values)
                    val = int(round(pct*total/100.0))
                    return '{v:d}'.format(v=val)
                    #return '{p:.0f}%'.format(p=pct)
                    #return '{p:.0f}% ({v:d})'.format(p=pct,v=val)
                return my_autopct
            
            colors, cellcnt, explode, labels = [], [], [], []
            members = [ getattr(CellPart, attr) for attr in dir(CellPart) if not callable(getattr(CellPart, attr)) and not attr.startswith("__") ]
            for clonality in members:
            #for clonality in [CellPart.CLONAL, CellPart.RELATED, CellPart.AMBIGUOUS, CellPart.BYSTANDERS, CellPart.NONCLONAL]:
                if clonality in d:
                    # start with certain main clone and end with certain non-clonal
                    #certain = [ t for t in d[clonality].items() if t[0][1] ]
                    #noncertain = [ t for t in d[clonality].items() if not t[0][1] ]
                    #if clonality == CellPart.NONCLONAL:
                    #    certain, noncertain = noncertain, certain
                        
                    #for key, cnt in certain + noncertain:
                    #print(d[clonality])
                    for key, cnt in d[clonality].items():
                        color, is_certain = key
                        label = clonality
                        #if clonality == CellPart.CLONAL:
                        #    label = 'main clone'
                        #elif clonality == CellPart.RELATED:
                        #    label = 'TC(main)'
                        #elif clonality == CellPart.AMBIGUOUS:
                        #    label = 'big bystanding clones'
                        #elif clonality == CellPart.BYSTANDERS:
                        #    #print('bipartite_tcrs.py: bystander: ')
                        #    label = 'small bystanding clones'
                        #elif clonality == CellPart.NONCLONAL:
                        #    label = 'bystander non-clonal'
                        #else:
                        #    label = ''
                        #if not is_certain:
                        #    label += '?'
                        labels.append(label)
                        cellcnt.append(cnt)
                        colors.append(color)
                        #explode.append(0.1 if clonality == 'UNKNOWN' else 0.05)
                        explode.append(0.04)
                        
            series = pd.Series(cellcnt, index=labels)
            unk = [ i for i, l in enumerate(labels) if l == '' ]
            if len(unk) > 0:
                labels[sum(unk) // len(unk)] = 'related to clone' #'other clones'
            if not series.name:
                series.name = ''
            patches, texts, autotexts = ax.pie(series,
                explode=explode,
                #subplots=True,
                colors=colors,
                #colormap=cmap, # 'Set1',
                labels=labels,
                #size=9,
                #use_index=False,
                #figsize=(7, 6),
                pctdistance=0.7, 
                startangle=90,
                counterclock=False,
                autopct=make_autopct(series.values))
            #centre_circle = plt.Circle((0,0),0.70,fc='white')
            #fig = plt.gcf()
            #fig.gca().add_artist(centre_circle)
            circle = Circle((0, 0), 0.4, fc='white') #, color='black')
            ax.add_patch(circle)
            ax.annotate('{}\ncells'.format(sum(series)), #weight='bold', 
                        xy=(0.0, 0.0), va='center', ha='center', color=cls.black) #, xycoords='axes fraction')
            
            
            #texts[0].set_fontsize(9)

            # text color
            #for text in texts:
            #    text.set_color('grey')
            #for autotext in autotexts:
            #    autotext.set_color('grey')
            
        series = pd.Series([len(cl) for cl in clones.values()], index=clones.keys())
        d = defaultdict(lambda: defaultdict(int))
        for clonality, cell in clones.items():
            for edge_color in cell:
                d[clonality][edge_color] += 1
        make_pie(ax, d) # TODO: use greek: αβ

    @classmethod
    def _plot_tcr_histogram(cls, batch, ax, *, chain, clonal_tcr):
        def get_tcr_distr(batch, chain, merge_nonclonal):
            def filter_nonclonal(distr):
                new_distr = defaultdict(int)
                nonclonal = 0
                for tcr, cnt in distr.items():
                    if merge_nonclonal and cnt <= 1:
                        nonclonal += cnt
                    else:
                        if tcr in new_distr.keys():
                            raise Exception('tcr in new_distr')
                        new_distr[tcr] = cnt
                return new_distr, nonclonal

            tcrs = batch.tcrs
            cells = tcrs.cells()
            orig_distr = defaultdict(int)
            for cell in cells:
                tcr_list = tcrs.cell2seqs(cell, chains=[chain])
                for tcr in tcr_list:
                    orig_distr[tcr] += 1

            distr, noncl = filter_nonclonal(orig_distr)
            series = pd.Series(list(distr.values()), index=list(distr.keys()), name=chain)
            series = series.sort_values(ascending=False)
            if merge_nonclonal:
                series = series.append(pd.Series(noncl, index=['non-clonal']))
            return series
    
        series = { chain: get_tcr_distr(batch, chain=chain, merge_nonclonal=False) for chain in ['alpha', 'beta'] }
    
        if len(series[chain]) == 0:
            Warning('No bipartite graph produced for {}: No values to plot.'.format(batch.tag))
            return
        tcr_cnt = [ (tcr, cnt) for tcr, cnt in series[chain].items() if cnt > 1 ]
        if not series[chain].empty:
            colors = [ cls.black if clonal_tcr != tcr else 'red' for tcr, cnt in series[chain].items() ]
            series[chain].plot.bar(ax=ax, color=colors)
            for i, (tcr, cnt) in enumerate(tcr_cnt):
                short_tcr = re.sub('_.*_', '..', tcr)
                c = 'grey' if clonal_tcr != tcr else 'red'
                if clonal_tcr == tcr or i == 0:
                    ax.text(i+0.25, cnt, short_tcr, color=c, ha='left', va='bottom', size=7)
        
        ax.set_xlabel(r'$\{}$-chains'.format(chain), color=cls.black)
        #ax.annotate(r'$\{}$-chains'.format(chain), xy=(0.5, 0.7), xycoords='axes fraction')
        ax.set_xlim(right=max(10, series['alpha'].size, series['beta'].size))
        ax.set_ylim(top=max(10, series['alpha'][0], series['beta'][0]))
        maxy = series[chain][0]
        ax.yaxis.set_ticks([maxy])
        ax.xaxis.set_ticks([]) #ax.xaxis.set_visible(False)  # removes the label too
        ax.set_facecolor('w')
        ax.grid(False)

    @classmethod
    def draw(cls, batch, *, out_fn, interactive, use_cellcolors):
        if not batch.tcrs:
            Warning('No bipartite graph for {}: no tcrs found.'.format(batch.tag))
            return

        start = time.time()
        
        plt.ioff()
        #fig = plt.figure() #ax = fig.add_subplot(111)
        import matplotlib.gridspec as gridspec
        gs = gridspec.GridSpec(3, 3)

        fig = plt.figure()
        #fig.tight_layout()
        ax = fig.add_subplot(gs[:, 0:2])
        ax_pie = fig.add_subplot(gs[0, 2])
        ax_hist_alpha = fig.add_subplot(gs[1, 2])
        ax_hist_beta = fig.add_subplot(gs[2, 2])
        
        #fig, (ax, ax_pie_alpha, ax_pie_beta) = plt.subplots(1, 2)
        #ax.set_xlim([0, 3])
        #ax.set_ylim([0, h])
        if batch.cell_part:
            cl = batch.cell_part.clonotype
            #print('bipartite_tcrs: clonotype: ', batch.cell_part)
            alpha, beta, rule, facs = cl['alpha'], cl['beta'], cl['rule'], cl['facs']
        else:
            alpha, beta, rule, facs = 'NoAlpha', 'NoBeta', 'NoRule', 'NoFACS'
        
        G, tcr2cellcnt, clones = cls._construct_graph(batch, use_cellcolors)
        
        tt['construct_graph'] = time.time() - start
        start = time.time()
        pos, w, h, d = cls._get_node_pos(G)
        cls._plot_bipartite(ax, G, pos, w, h, d, alpha, beta, tcr2cellcnt)
        #cls._add_annot(fig, ax, G, tcr2cellcnt, pos)
        
        tt['plot_bipartite'] = time.time() - start
        start = time.time()
        cls._plot_clonality_piechart(clones, ax_pie)
        tt['plot_piechart'] = time.time() - start
        start = time.time()
        cls._plot_tcr_histogram(batch, ax_hist_alpha, chain='alpha', clonal_tcr=alpha)
        cls._plot_tcr_histogram(batch, ax_hist_beta, chain='beta', clonal_tcr=beta)
        tt['plot_tcr_histograms'] = time.time() - start
        start = time.time()
       
        #ax.axis('equal')
        #ax.relim()
        #ax.autoscale_view()
        
        fig.suptitle('{} ({})'.format(batch.tag, batch.internal_name), fontsize=12, weight='bold')
        if batch.cell_part:
            alpha_letter = r'$\alpha$'
            beta_letter = r'$\beta$'
            tag = 'clonality rule: {}\n{}: {}, {}: {}\nFACS: {}'.format(rule, alpha_letter, alpha, beta_letter, beta, facs)
        else:
            tag = 'No clone specified'
        fig.text(0.5, 0.89, tag, ha="center", va="center", color='r', size=7)
            #ax.set_title('({}) main clone:\n{}'.format(chain, tag), fontsize=10, color='r')
        
        #Path(out_fn).parent.mkdir(parents=True, exist_ok=True)
        #print("Figure saved as {}".format(blue(Path(out_fn).name)))
        #fig.savefig(str(out_fn))
        #plt.close()
        plot(plt, interactive, out_fn)
