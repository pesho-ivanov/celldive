import pandas as pd
from utils import *

# Machine learning
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.manifold import TSNE
from sklearn.manifold import MDS

#from ZIFA import ZIFA
#from ZIFA import block_ZIFA

# Assumes that the probability density functions for every class are both normally distributed
def getLda2D(X, y):
    lda = LinearDiscriminantAnalysis(n_components=2)
    X_r2 = lda.fit(X, y).transform(X)
    return pd.DataFrame(lda.transform(X), index=X.index, columns=['x', 'y'])

#labels = meta['sample']
#lda_all = getLda2D(quant, labels)

#lda = LinearDiscriminantAnalysis(n_components=2)
#X_r2 = lda.fit(quant, labels).transform(quant)
#processAndDraw2D(lda_all, meta, 'LDA')

# Reduct the dimentionality from D-dimentional array of n points X to d-dimentional
def getDimReduct(X, method, perplexity, dim):
    assert len(X) > 0
    print("DimReduct ({}) of {} cells with {} dimensions...".format(
        green(method), blue(X.shape[0]), blue(X.shape[1])))
    #print("shape: ", X.shape)
    
    var = None
    if method == 'PCA':  # Principle component analysis: O((n+d)D^2)
        latent = PCA(n_components=dim).fit(X)
        var = latent.explained_variance_ratio_
        perc = map(lambda x: blue(int(x*100))+'%', var)
        print ("PCA explained variance: ", " ".join(perc))
        res = latent.transform(X)
        
    elif method == 'MDS':  # Multi-dimentional scaling: O((D+d)n^2)
        res = MDS(n_components=dim).fit_transform(X)
        #latent = MDS(n_components=dim).fit(X)
        #print ("components: ", latent.components_)
        #res = latent.transform(X)
        
    elif method == 't-SNE':  # t-distributed stochastic neighbor embedding
        # tutorial: http://distill.pub/2016/misread-tsne
        res = TSNE(n_components=dim, random_state=0, perplexity=perplexity).fit_transform(X)
        
    elif method == 'DUMB':  # t-distributed stochastic neighbor embedding
        res = X.ix[:, :dim].values
        
    #elif method == 'ZIFA':  # domain specific: Zero Inflated Factor Analysi
        # TODO: run
        #logX = np.log2(1 + X.as_matrix())    # as recommended on https://github.com/epierson9/ZIFA
        #print(logX.shape)
        #print(logX)
        
        #print(X)
        #res, model_params = ZIFA.fitModel(X, 2)
        #print(model_params)
        
        #Z, model_params = block_ZIFA.fitModel(logX, 2, n_blocks=10)
    else:
        print("No such dimentionality reduction method!")
        assert(False)
   
    # TODO
    df = pd.DataFrame(res, index=X.index)
    if dim <= 3:
       df.columns = ['x', 'y', 'z'][:dim]
    return df, var

#processAndDraw2D(quant.ix[:,:100], 'ZIFA', ' bla')

# TODO: just get the top_n genes with maximal variance. No PCA!
output_base = '~/work/BioViz-gl/output/'
def getTopGenesPCA(quant, top_n, output_file=None):
    pca = PCA(n_components=2).fit(quant)
    abs_var_sum = pd.DataFrame(pca.components_.T, index=quant.columns).abs().sum(axis=1)
    abs_var_sum.sort_values(0, ascending=False, inplace=True)
    top_genes = abs_var_sum[:top_n].index
    if output_file:
        top_quant = quant[top_genes]#.T
        top_quant.to_csv(os.path.join(output_base, output_file))
    return top_genes

def getGenesWithMostVar(quant, top_n):
    v = np.var(np.array(quant), 0)
    assert(len(v) == quant.shape[1])
    return quant.columns[np.argsort(v)[-top_n:]]

def printSampleCells():
    sample = meta['sample']
    SL = sample[sample.isin(['SL'])]
    print('\n'.join(SL.index))
    
def printSynthesisParametersFile(output_base):
    #genes = getTopGenesPCA(quant, 60, None)
    genes = getGenesWithMostVar(quant, 10)
    print("genes: ", genes)
    params = pd.DataFrame([[2,2,100]]*len(genes),
                          index=genes,
                          columns=['MaxActivators', 'MaxRepressors', 'Threshold'])
    print(params.head())
    top_quant = quant[genes]#.T
    params.to_csv(os.path.join(output_base, 'params.csv'))
    top_quant.to_csv(os.path.join(output_base, 'top-quant.csv'))
    
    #fout = open(outfile, 'w')
    #print 'Gene,MaxActivators,MaxRepressors,Threshold'
    #rows = [ cell + ',1,1,100' for cell in quant.index ]
    #return '\n'.join(rows)

#printSynthesisParametersFile(output_base)

#top_genes = getTopGenesPCA(quant, 60, 'top60-t.csv')
#print top_genes
#print quant.head()
#print quant[top_genes].head().T

#def isCellBE(metaElem):
#    return metaElem['sample'] == 'BE'