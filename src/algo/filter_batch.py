class FilterBatch:
    @staticmethod
    def filter_cells(batch, filt_cells, filt_genes):
        expressed_genes = filt_cells.get('expressed_genes')
        if expressed_genes:
            s = (self.ge >= expression_threshold).sum(axis=0)
            print(s.head())
            #cells = batch.quant.ge[ filt_cells ]
        
        return batch.subbatch(cells=cells, features=genes)
        
        
    #def filter_genes(filt_genes):
        