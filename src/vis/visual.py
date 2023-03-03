#from bokeh.charts import Scatter, output_file, show
#from bokeh.sampledata.autompg import autompg as df
#from bokeh.plotting import figure, output_notebook, show
#from bokeh.models import ColumnDataSource, Range1d, LabelSet, Label
#from bokeh.models import HoverTool
#output_notebook()  # bokeh plotting in the notebook

# matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# 3D; register Axes3D class with matplotlib by importing Axes3D
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

# Hover tool
#from mpld3 import plugins
#import mpld3
#mpld3.display()

import seaborn as sns
import random

from utils import *

#fig.savefig(out_dir + '/%s.svg' % (name + '-' + str(len(points))))
#plt.show()

def get_3d_cell(ax3d, x, y, z, r, c):
    c = Circle((x, y), radius=r, color=c, alpha=0.5)
    patch = ax3d.add_patch(c)
    art3d.patch_2d_to_3d(patch, z=z, zdir="z")
    
def plot_3d_clonogram(ax3d, coords_2d):
    alphas, betas = zip(*coords_2d)
    
    ax3d.set_xlim3d(0, 100)
    ax3d.set_ylim3d(0, 100)
    ax3d.set_zlim3d(0, 40)
    #p(max(alphas + betas))
    #p(alphas)
    
    used = []
    patches = [] # [Circle((3,4))]
    for alpha, beta in coords_2d: 
        cnt = used.count((alpha, beta))
        used.append((alpha, beta))
        #patches.append(Circle((alpha, beta), radius=r, color='b', alpha=0.5))
        get_3d_cell(ax3d, x=beta, y=alpha, z=cnt, r=2, c='b')
        
    #patches = ax3d.scatter(x=alphas, y=betas, s=200, c=betas, linewidths=0.1, cmap=plt.get_cmap('nipy_spectral'), alpha=0.4)
    #art3d.patch_collection_2d_to_3d(patches, z=30, zdir='z')
    #collection = art3d.Patch3DCollection(patches, zs=0, zdir="z")
    #p(patches)
    #collection = PatchCollection(patches)
    #ax3d.add_collection3d(collection, zs=0, zdir='z')
    ax3d.set(title=str(len(alphas)) + ' cells\' TCRs', xlabel='beta chain', ylabel='alpha chain', xticks=[], yticks=[])
    ax3d.set_facecolor('white')
    handles, labels = ax3d.get_legend_handles_labels()
    #ax3d.view_init(elev=20., azim=-35)

def get_coords(recombs):
    tcrs = []
    for cell_chains in recombs:
        for chain in cell_chains:
            tcrs.append(chain)
    tcra = [chain['recomb_id'] for chain in tcrs if chain['locus']=='A' and chain['prod']]
    tcra = sorted(set(tcra))
    tcrb = [chain['recomb_id'] for chain in tcrs if chain['locus']=='B' and chain['prod']]
    tcrb = sorted(set(tcrb))

    coords = []  # a pair of positions in tcra and tcrb arrays
    for cell_chains in recombs:
        if len(cell_chains) == 2:
            a, b = -1, -1
            for chain in cell_chains:
                if chain['recomb_id'] in tcra:
                    a = tcra.index(chain['recomb_id'])
                if chain['recomb_id'] in tcrb:
                    b = tcrb.index(chain['recomb_id'])
            if min(a, b) > -1:
                coords.append((a, b))
    return coords

def calc_and_plot_clonogram(recombs):
    coords = get_coords(recombs)
    fig = plt.figure(figsize=(20,10))
    ax2d = fig.add_subplot(121)
    plot_2d_clonogram(ax2d, coords)
    ax3d = fig.add_subplot(122, projection='3d')
    #ax3d = Axes3D(fig)
    #plot_25d_clonogram(ax3d, coords)
    plot_3d_clonogram(ax3d, coords)

def getEdges(points, recombinants):
    edges = []
    recomb2cell = dict()
    for cell in points.index:
        for recomb in recombinants[cell]:
            recomb2cell.setdefault(recomb['recomb_id'], []).append(cell)
    for recomb, cells in recomb2cell.items():
        for i, cellA in enumerate(cells):
            for cellB in cells[:i]:
                assert cellA in points.index
                edges.append(tuple(points.loc[cellA]) + tuple(points.loc[cellB]) + (recomb,))
    return edges

def getTCRs(cell_meta):
    if cell_meta.recombinants:
        return ' | '.join([recomb['recomb_id'] for recomb in cell_meta.recombinants])
    else:
        return ''
    
def scatter3D(ax3d, points, meta, name, edges):
    colormap = { 'BE': sns.xkcd_rgb["light brown"],     #ad8150
                 'BL': sns.xkcd_rgb["dark brown"],      #341c02
                 'SE': sns.xkcd_rgb["pastel yellow"],   #fffe71
                 'SL': sns.xkcd_rgb["goldenrod"],       #fac205
                 'TCRA': sns.xkcd_rgb["medium blue"],   #2c6fbb
                 'TCRB': sns.xkcd_rgb["greyish green"] }#82a67d
    #colormap = { 'BE': '#ff9999', 'BL': '#ff3236', 'SE': '#ffde99', 'SL': '#ffc132' }
    colors = [colormap[x] for x in meta.ix[points.index]['sample']]
    
    ax3d.set_title(name + ' over ' + str(len(points)) + ' cells')
    #ax3d.grid(False)
    #ax3d.axis('off')
    
    ax3d.scatter(points['x'], points['y'], zs=points['z'], zdir='z', c=colors, s=60, edgecolors='face', alpha=0.7) #, label='points in (x,z)')
    ax3d.view_init(elev=20., azim=-35)
    
    # labels
    labels = ['#'+cell+': '+getTCRs(meta.loc[cell]) for cell in points.index]

def dendroplot(D):
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('dendrogram.png')
    
def plot_2d_clonogram(ax, coords):
    alphas, betas = zip(*coords)
    ahand = ax.scatter(x=alphas, y=betas, s=200, c=betas, linewidths=0.1, cmap=plt.get_cmap('nipy_spectral'), alpha=0.4)
    bhand = ax.scatter(x=alphas, y=betas, s=50, c=alphas, linewidths=0, cmap=plt.get_cmap('gist_ncar'), alpha=0.4)
    ax.set(title=str(len(alphas)) + ' cells\' tcrs', xlabel='beta chain', ylabel='alpha chain', xticks=[], yticks=[])
    ax.set_facecolor('white')
    #ax.set_xticks(range(len(tcrb))[::4], tcra[::4], rotation='vertical')
    #ax.set_yticks(range(len(tcra))[::4], tcra[::4], rotation='horizontal')
    #plt.legend()#[p2, p1], ["line 2", "line 1"])
    handles, labels = ax.get_legend_handles_labels()
    #p(ahand.get_array())
    #p(mpl.artist.getp(ahand))
    #p(handles, labels)
    #legend = ax.legend((ahand,bhand), ('ala','bala'))
    # cell labels
    for x, y in set(coords):
        cnt = coords.count((x,y))
        if cnt > 1:
            ax.annotate(str(cnt), (x+0.5,y+0.5), size='large')
    return ax

#def plot_25d_clonogram(ax, coords):
#    alphas, betas = zip(*coords)
#    ahand = ax.scatter(xs=alphas, ys=betas, zs=0, zdir='z', s=200, c=betas, linewidths=0.1, cmap=plt.get_cmap('nipy_spectral'), alpha=0.4)
#    #bhand = ax.scatter(x=alphas, y=betas, zs=0, zdir='z', s=50, c=alphas, linewidths=0, cmap=plt.get_cmap('gist_ncar'), alpha=0.4)
#    ax.set(title=str(len(alphas)) + ' cells\' tcrs', xlabel='beta chain', ylabel='alpha chain', xticks=[], yticks=[])
#    ax.set_facecolor('white')
#    return ax

def get_3d_cell(ax3d, x, y, z, r, c):
    c = Circle((x, y), radius=r, color=c, alpha=0.5)
    patch = ax3d.add_patch(c)
    art3d.patch_2d_to_3d(patch, z=z, zdir="z")
    
def plot_3d_clonogram(ax3d, coords_2d):
    alphas, betas = zip(*coords_2d)
    
    ax3d.set_xlim3d(0, 100)
    ax3d.set_ylim3d(0, 100)
    ax3d.set_zlim3d(0, 40)
    #p(max(alphas + betas))
    #p(alphas)
    
    used = []
    patches = [] # [Circle((3,4))]
    for alpha, beta in coords_2d: 
        cnt = used.count((alpha, beta))
        used.append((alpha, beta))
        #patches.append(Circle((alpha, beta), radius=r, color='b', alpha=0.5))
        get_3d_cell(ax3d, x=beta, y=alpha, z=cnt, r=2, c='b')
        
    #patches = ax3d.scatter(x=alphas, y=betas, s=200, c=betas, linewidths=0.1, cmap=plt.get_cmap('nipy_spectral'), alpha=0.4)
    #art3d.patch_collection_2d_to_3d(patches, z=30, zdir='z')
    #collection = art3d.Patch3DCollection(patches, zs=0, zdir="z")
    #p(patches)
    #collection = PatchCollection(patches)
    #ax3d.add_collection3d(collection, zs=0, zdir='z')
    ax3d.set(title=str(len(alphas)) + ' cells\' TCRs', xlabel='beta chain', ylabel='alpha chain', xticks=[], yticks=[])
    ax3d.set_facecolor('white')
    handles, labels = ax3d.get_legend_handles_labels()
    #ax3d.view_init(elev=20., azim=-35)

def get_coords(recombs):
    tcrs = []
    for cell_chains in recombs:
        for chain in cell_chains:
            tcrs.append(chain)
    tcra = [chain['recomb_id'] for chain in tcrs if chain['locus']=='A' and chain['prod']]
    tcra = sorted(set(tcra))
    tcrb = [chain['recomb_id'] for chain in tcrs if chain['locus']=='B' and chain['prod']]
    tcrb = sorted(set(tcrb))

    coords = []  # a pair of positions in tcra and tcrb arrays
    for cell_chains in recombs:
        if len(cell_chains) == 2:
            a, b = -1, -1
            for chain in cell_chains:
                if chain['recomb_id'] in tcra:
                    a = tcra.index(chain['recomb_id'])
                if chain['recomb_id'] in tcrb:
                    b = tcrb.index(chain['recomb_id'])
            if min(a, b) > -1:
                coords.append((a, b))
    return coords

def calc_and_plot_clonogram(recombs):
    coords = get_coords(recombs)
    fig = plt.figure(figsize=(20,10))
    ax2d = fig.add_subplot(121)
    plot_2d_clonogram(ax2d, coords)
    ax3d = fig.add_subplot(122, projection='3d')
    #ax3d = Axes3D(fig)
    #plot_25d_clonogram(ax3d, coords)
    plot_3d_clonogram(ax3d, coords)
   


# flights = sns.load_dataset("flights")
# flights = flights.pivot("month", "year", "passengers")
# g = sns.clustermap(flights)
# #pesho_visual.dendroplot(batch3.quant)
#     
# a = batch3.quant.ix[:100, :100]
# #p(a)
# #p(a.values.tolist())
# ##a = [item for sublist in a for item in sublist]
# #K = 1000
# g = sns.clustermap(a)
# #g.savefig('dendrosns.png')
# #sns.distplot(a, kde=False, rug=True);