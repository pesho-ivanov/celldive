import matplotlib.pyplot as plt

from .violin import Violin
from utils import *

class BatchStats:
    
    def __init__(self, batch):
        self.batch = batch
        self.ge = batch.quant.ge
    
    def draw_genes_vs_expression(self, ax, *, expression_threshold):
        ax.set_yscale('log')
        ax.set_ylabel('genes', fontsize=10)
        ax.set_xlabel('expression', fontsize=10)
        s = pd.Series(self.ge.values.ravel())
        ax.axvline(expression_threshold, color='r', linestyle='--')
        s.hist(ax=ax, bins=20, color='b')
        ax.text(0.95*expression_threshold, ax.get_ylim()[1],
                'min expression: {}'.format(expression_threshold),
                verticalalignment='top', horizontalalignment='right', rotation=90, color='r')
        percentage = 100.0 * s[s >= expression_threshold].sum() / s.sum()
        ax.text(1.05*expression_threshold, 0.95*ax.get_ylim()[1],
                '{:.2f}%'.format(percentage),
                verticalalignment='top', horizontalalignment='left')
    
    def draw_cells_vs_expressed_genes(self, ax, *, expression_threshold, expressed_genes):
        ax.set_ylabel('cells', fontsize=10)
        ax.set_xlabel('expressed genes (>={})'.format(expression_threshold), fontsize=10)
        ax.axvline(expressed_genes, color='r', linestyle='--')
        s = (self.ge >= expression_threshold).sum(axis=1)
        s.hist(ax=ax, bins=20, color='g')
        ax.text(0.95*expressed_genes, ax.get_ylim()[1],
                'min expressed genes: {}'.format(expressed_genes),
                verticalalignment='top', horizontalalignment='right', rotation=90, color='r')
        percentage = 100.0 * s[s >= expressed_genes].size / len(self.ge.index)
        ax.text(1.05*expressed_genes, 0.95*ax.get_ylim()[1],
                #'{} cells'.format(s.size), horizontalalignment='left')
                '{:.2f}% ({} cells)'.format(percentage, s.size),
                verticalalignment='top', horizontalalignment='left')
        
    def draw_genes_vs_expressing_cells(self, ax, *, expression_threshold, min_cells_expressing_gene):
        ax.set_yscale('log')
        ax.set_ylabel('expressed genes (>={})'.format(expression_threshold), fontsize=10)
        ax.set_xlabel('cells', fontsize=10)
        ax.axvline(min_cells_expressing_gene, color='r', linestyle='--')
        s = (self.ge >= expression_threshold).sum(axis=0)
        s.hist(ax=ax, bins=20, color='xkcd:lavender')
        ax.text(0.95*min_cells_expressing_gene, ax.get_ylim()[1],
                'min expressed genes: {}'.format(min_cells_expressing_gene),
                verticalalignment='top', horizontalalignment='right', rotation=90, color='r')
        ax.text(1.05*min_cells_expressing_gene, 0.95*ax.get_ylim()[1],
                #'{} cells'.format(s.size), horizontalalignment='left')
                '{} genes'.format(s[s>=min_cells_expressing_gene].size), verticalalignment='top', horizontalalignment='left')
        
    ### Coeficient of variation:  var(X) / mean(X)
    
    @staticmethod
    def coef_of_var(ge, axis):
        if axis == 'genes':
            axis = 'rows'
        elif axis == 'cells':
            axis = 'columns'
            
        return ge.var(axis=axis) / ge.mean(axis=axis)
    
    @staticmethod
    def remove_ticks(ax):
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') 
        
    def draw_cells_var(self, ax, *, ge):
        ax.set_ylabel('sorted cell vars (nonzero expression)')
        ax.set_xlabel('cells')
        ge_nonzero = ge.replace(0, np.NaN)
        sorted_genes = self.coef_of_var(ge_nonzero, axis='cells').sort_values()
        del sorted_genes.index.name
        sorted_genes.plot(ax=ax) 
        self.remove_ticks(ax)
        
    def draw_genes_var(self, ax, *, ge):
        ax.set_ylabel('sorted gene vars (nonzero expression)')
        ax.set_xlabel('genes')
        ge_nonzero = ge.replace(0, np.NaN)
        sorted_genes = self.coef_of_var(ge_nonzero, axis='genes').sort_values()
        del sorted_genes.index.name
        sorted_genes.plot(ax=ax) 
        self.remove_ticks(ax)
        
    def draw_cells_per_number_of_nonzero_genes(self, ax, *, ge):
        ax.set_ylabel('#cells with nonzero expression of a gene')
        ax.set_xlabel('genes sorted by number of containing cells')
        cells = (ge > 0).sum(axis='rows').sort_values()
        del cells.index.name
        cells.plot(ax=ax)
        self.remove_ticks(ax)
        
    def draw_var_of_genes_expressed_in_all_cells(self, ax, *, ge):
        number_of_cells = (ge > 0).sum(axis='rows')
        maxval = number_of_cells.sort_values()[-1]
        genes = number_of_cells[number_of_cells == maxval].index
        sorted_genes = self.coef_of_var(ge[genes], axis='genes').sort_values()
        del sorted_genes.index.name
        sorted_genes.plot(ax=ax) 
        ax.set_ylabel('var')
        xlabel = '{} genes expressed in maximum #cells ({})'.format(green(len(genes)), green(len(ge.index)))
        ax.set_xlabel(xlabel)
        
    def draw_all(self, *, out_fn=None, interactive=False, filt_cells, filt_genes):
        if self.ge is None or self.ge.empty:
            warning('Empty gene list. No stats outputted.')
            return
        
        expression_threshold = filt_genes['expression_threshold']
        min_cells_expressing_gene = filt_genes['min_cells_expressing_gene']
        expressed_genes = filt_cells['expressed_genes']
        
        #plt.style.use('fivethirtyeight')
        plt.style.use('ggplot')
        #plt.style.use('bmh')
        #plt.axis('off')
 
        plt.ioff()
        
        fig, ((ax11, ax12, ax13, ax14), (ax21, ax22, ax23, ax24)) = plt.subplots(2, 4)
        title = '{} ({} cells x {} genes)'.format(self.batch.tag, len(self.batch.cells()), len(self.batch.genes()))
        fig.suptitle(title, fontsize=12, weight='bold')
        
        # row I
        self.draw_genes_vs_expression(ax11, expression_threshold=expression_threshold)
        self.draw_cells_vs_expressed_genes(ax12, expression_threshold=expression_threshold, expressed_genes=expressed_genes)
        self.draw_genes_vs_expressing_cells(ax13, expression_threshold=expression_threshold, min_cells_expressing_gene=min_cells_expressing_gene)
        
        # row II
        self.draw_cells_var(ax21, ge=self.ge)
        self.draw_genes_var(ax22, ge=self.ge)
        self.draw_cells_per_number_of_nonzero_genes(ax23, ge=self.ge)
        self.draw_var_of_genes_expressed_in_all_cells(ax24, ge=self.ge)
        
        spikein_genes = [ gene for gene in self.batch.genes() if gene.startswith('Spike')]
        if spikein_genes:
            spikein_genes_df = self.ge[spikein_genes]
            Violin.draw(spikein_genes, out_fn='violin.png')  # TODO
        else:
            Warning('No spikein genes found!')
        
        plot(plt, interactive, out_fn) 
        