import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# rpy2
import rpy2.robjects as robj
from rpy2.robjects import numpy2ri
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

from utils import *

stats = importr("stats")
gplots = importr("gplots")

numpy2ri.activate()
pandas2ri.activate()
R = robj.r
rprint = robj.globalenv.get("print")

class QuantHeatmap:
#    @staticmethod
#    def draw_dendrogram(batch):
#        #ge = pd.DataFrame(np.random.randn(6,4))
#        ge = batch.quant.ge.iloc[:5,:5]  # (cells, genes)
#        #print(ge)
#        corrm = R.cor(R.t(ge), method="spearman")
#
#        # no idea how to do the 1-matrix automagically
#        robj.globalenv['corrm'] = corrm
#        R('corrm <- corrm[!is.na(corrm)]')
#        rprint(corrm)
#        distm = stats.as_dist(R("1-corrm"))
#
#        #rprint(distm)
#        
#        # do the clustering
#        hr = stats.hclust(distm, method = "complete", members = R("NULL"))
#        
#        R.png(file='/tmp/pesho/out_dendro.png')
#        
#        #R.heatmap(ge.values)
#        #R.par(mfrow = R.c(1,2))
#        R.plot(hr, hang = 0.1)
#        #R.plot(hr, hang = -0.1)
#
#        R("dev.off()")
    
    @staticmethod
    def draw(batch, out_fn=None, *, get_color):
        ge = batch.quant.ge
            
        #ge = batch.quant.ge.iloc[:,:100]
        #corrm = R.cor(R.t(ge), method="spearman")
        #print(corrm)
        #distm = stats.as_dist(R("1-corrm"))
        #print(distm)
        
        #if len(ge.columns) > 100:
        #    print('Warning: Too many good genes (%d). Cutting!' % len(ge.columns))
        #    ge = ge[:100]
        
        #genes = batch.clones.best_hypo().best_genes  # series: gene->padj
        
        if len(ge.index) <= 2 or len(ge.columns) <= 2:
            print('Warning: Quant heatmap cannot be produced! Too small matrix {}x{}.'.format(
                len(ge.index), len(ge.columns)))
            return
        
        # pandas -> R
        #important_genes_with_pvals = batch.get_important_genes2pvals()
        r_matrix = R.matrix(ge.values, nrow=len(ge.index))
        r_matrix.rownames = robj.StrVector(ge.index.tolist())
        #def to_gene_names(t): tid2tname(t) + ': ' + pval2str(important_genes_with_pvals[t])
        #gene_names = list(map(to_gene_names, ge.columns.tolist()))
        gene_names = ge.columns.tolist()
        
        r_matrix.colnames = robj.StrVector(gene_names)
        
        #fn = '../out/{}_heatmap_{}_{}.png'.format(batch.tag, len(batch.genes()), batch.get_onlytag())
        if out_fn:
            R.png(file=str(out_fn), width=1000, height=1000, res=400, pointsize=5)
            print("QuantHeatmap is saved in {}".format(green(str(out_fn))))
            #cell_colors = list(map(lambda cell: batch.meta.get_color_rgb(cell), ge.index.tolist()))
            cell_colors = list(map(lambda cell: to_3rgb_floats(get_color(cell)), ge.index.tolist()))
            categories = list(map(lambda cell: get_color(cell, annot=True), ge.index.tolist()))
            #print(cell_colors)
            
            tag = batch.get_onlytag() if batch.get_onlytag else 'No clone'
            title = '{}\n{}'.format(batch.tag, tag)
            
            gplots.heatmap_2(
                r_matrix,
                main = title,
                xlab = 'genes', ylab='cells',
                RowSideColors = robj.StrVector(cell_colors),
                legend = robj.StrVector(categories),
                #Colv = robj.NA_Logical,
                scale = 'column', # 'none',  # 'row' 
                dendrogram='both',
                col = R("rev(heat.colors(10))"),  # R('greenred(10)') # heat.colors(n),

                # specific to heatmap_2
                trace='none',
                #distfun=R('function(x) as.dist((1 - cor(x))/2)'),
            )  # R('c("red","blue","blue")')
            R("dev.off()")
