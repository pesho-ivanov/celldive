import numpy as np
import matplotlib.pyplot as plt

from utils import *

class DiversityTCR:
    @classmethod
    def _construct_graph(cls, batch, use_cellcolors):
        '''Returns a multigraph'''
        G = nx.MultiGraph()

        tcrs = batch.tcrs
        cells = tcrs.cells()
        tcrs_A = tcrs.chain2tcrs('A')
        tcrs_B = tcrs.chain2tcrs('B')

        G.add_nodes_from(tcrs_A, side='left')
        G.add_nodes_from(tcrs_B, side='right')
        G.add_nodes_from([cls._source, cls._sink], side='no')

        cnt = 0
        tcr2cellcnt = Counter()
        for cell in cells:
            alphas, betas = tcrs.cell2tcrs(cell)
            lena, lenb = len(alphas), len(betas)
            
            for tcr in [*alphas, *betas]:
                tcr2cellcnt[tcr] += 1
            
            if lena > 0 and lenb > 0:
                edges = [ (a, b) for b in betas for a in alphas]
            elif lena > 0:
                edges = [ (cls._source, a) for a in alphas ]
            else:
                edges = [ (b, cls._sink) for b in betas ]
            cnt += len(edges)
            
            if use_cellcolors:
                color = batch.get_color_rgb(cell)
            else:
                if lena == 1 and lenb == 1:
                    color = 'black'
                else:
                    color = random_color()
            G.add_edges_from(edges, color=color, cell=cell)
            
        edges = G.edges()
        return G, tcr2cellcnt
    
    @classmethod
    def _get_node_pos(cls, G, w, h, d):
        # positions
        left =  [ node for node in G.nodes() if G.node[node]['side'] == 'left' ]
        right = [ node for node in G.nodes() if G.node[node]['side'] == 'right' ]
        
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
                
        assert len(G.nodes()) == len(pos)
        return pos
    
    @classmethod
    def rainbow_edge(cls, ax, p, q, colors=_rainbow_colors, spacing=0.1):
        n = 10
        base = np.linspace(0.1, np.pi-0.1, n)
        shape = np.sin(base)
        x = np.linspace(p[0], q[0], n)
        y = np.linspace(p[1], q[1], n)
        
        delta = len(colors) * spacing / 2.0
        lines = []
        for i in range(len(colors)+1):
            newline_data = y - 0.00002*spacing*i ## ?!?!
            coef = len(colors)/2.0 - i
            newline_data += coef * spacing * shape
            lines.append(newline_data)
        for i, c in enumerate(colors):
            ax.fill_between(x, lines[i], lines[i+1], facecolor=c)

    @classmethod
    def _mpl_plot(cls, ax, G, pos, w, h, d, tcr2cellcnt):
        # TODO: colors
        # http://matplotlib.org/users/transforms_tutorial.html#using-offset-transforms-to-create-a-shadow-effect
        # http://stackoverflow.com/questions/42165631/in-matplotlib-how-can-i-plot-a-multi-colored-line-like-a-rainbow/42190453#42190453
        ax.grid(False)
        ax.axis('off')
    
        for node in G.nodes():
            x, y = pos[node]
            circle = Circle((x, y), d)
            ax.add_patch(circle)
            cnt = tcr2cellcnt[node]
            if cnt > 1:
                t = ax.text(x, y, str(cnt), color='w', ha='center', va='center', weight='bold')
                t.set_bbox(dict(facecolor='black', alpha=0.2))
                
        # column annotations
        x, y = pos[cls._source]
        ax.text(x, y*1.05, 'no β', ha='center', va='top')
        x, y = pos[cls._sink]
        ax.text(x, y*1.05, 'no α', ha='center', va='top')
        ax.text(1*w/3.0, h, 'TCR-α', ha='center', va='top')
        ax.text(2*w/3.0, h, 'TCR-β', ha='center', va='top')
            
        edgelist = G.edges()
        edge_cnt = Counter(edgelist)
        for a, b in edgelist:
            edges=G[a][b]
            edgecolors = [ e['color'] for e in G[a][b].values() ]
            edgecolors = np.sort(edgecolors)
            cls.rainbow_edge(ax, pos[a], pos[b], colors=edgecolors)

        #ax.scatter = ax.scatter(points['x'], points['y'], c=points['color'], s=60, edgecolors='face', alpha=0.7, zorder=2);
        
    @classmethod
    def _add_annot(cls, fig, ax, G, pos):
        # interactive annotations on click
        X = []
        Y = []
        annotes = []
        for tcr, coord in pos.items():
            x, y = coord
            X.append(x)
            Y.append(y)
            annotes.append(tcr)
            
        af = AnnoteFinder(X, Y, annotes, ax=ax)
        fig.canvas.mpl_connect('button_press_event', af)

    @classmethod
    def draw(cls, batch, use_cellcolors=False):
        fig, ax = plt.subplots(1, 1)
        #ax.set_xlim([0, 3])
        #ax.set_ylim([0, h])
        
        w = 10
        h = 10
        d = 0.2
        
        G, tcr2cellcnt = cls._construct_graph(batch, use_cellcolors)
        pos = cls._get_node_pos(G, w, h, d)
        cls._mpl_plot(ax, G, pos, w, h, d, tcr2cellcnt)
        cls._add_annot(fig, ax, G, pos)
        
        ax.axis('equal')
        #ax.relim()
        #ax.autoscale_view()

        plt.show()
