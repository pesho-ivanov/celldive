import numpy as np
from utils import *

class QualityControl:
    cell_suffixes = ['', '_1', '_2']
    """To aggregate pair read information."""
    
    def __init__(self, table):
        """Class for working with quality controls per cell.
        Initialize using from table or fromDir.
        """
        self.table = table
        """A pandas table with pairs of cells as rows and quality control information as columns."""
        
        assert len(self.table) > 0
        
    @classmethod
    def from_table(cls, table):
        """Sample  adapter_content Sequences flagged as poor quality       sequence_duplication_levels \
        avg_sequence_length     Encoding        kmer_content    per_base_sequence_quality \
        sequence_length_distribution    Sequence length File type       basic_statistics \
        per_sequence_gc_content Total Sequences per_base_n_content      per_base_sequence_content \
        overrepresented_sequences       %GC     total_deduplicated_percentage   Filename \
        per_tile_sequence_quality       per_sequence_quality_scores
        
        Example: 100_1   fail    0.0     fail    151.0   Sanger / Illumina 1.9   fail    pass    pass    151.0   Conventional base calls pass    fail    1789948.0       pass    warn    warn    50.0    22.9682959994   100_1.      fastq.gz        pass    pass
        """
        return cls(table)
    
    @staticmethod
    def merge_pair_reads(pairs_table): 
        def cell_name_without_orientation_number(cell):
            if cell.endswith('_1') or cell.endswith('_2'):
                assert len(cell) > 2
                return cell[:-2]
            else:
                return cell
            
        groups = pairs_table.groupby(cell_name_without_orientation_number)
        return groups.first()
    
    @classmethod
    def from_file(cls, true_index, input_file, sep='\t', *, tag):
        if not check_file(input_file, "Quality control", silent=True):
            return None
        
        #column_types = { 'Sample': str, }
        pairs_table = pd.read_csv(input_file, sep=sep, index_col=0)
                                  #dtype=column_types)
        pairs_table.index.name = 'cell'
        pairs_table.index = pairs_table.index.astype(str)
        #columns_rename = { }
        #pairs_table.rename(columns=columns_rename, inplace=True)
        pairs_table.sort_index(inplace=True)
        
        table = QualityControl.merge_pair_reads(pairs_table)
        
        add_prefix_tag(table, tag)
        assert_same_indexes(true_index, table.index, 'Quality control')
        #true_df = pd.DataFrame(index=true_index)
        #pairs_table = pd.concat([true_df, pairs_table], axis=1)
        return cls(table=table)
    
    def subqc(self, cells):
        # TODO!
        #    return QualityControl(self.table.loc[cells])
        return QualityControl(self.table)
    
    def concat(self, other, self_rename, other_rename):
        renamed_qc = []
        for table, rename_f in [(self.ge, self_rename),
                                (other.ge, other_rename)]:
            a = table.rename(index=rename_f)
            renamed_qcs.append(a)
        
        return QualityControl(pd.concat(renamed_qcs))
    
    def _allcells(self, cell):
        """Aggregating pair-end reads."""
        cells = [ cell + suffix for suffix in self.cell_suffixes if cell + suffix in self.table.index ]
        if len(cells) == 0: raise ValueError('No entry for cell ' + str(cell))
        return cells
        
    def reads(self, cell):
        """Returns the number of single-end reads (pair-end reads count twice)."""
        return sum([ self.table.loc[c, 'Total Sequences'] for c in self._allcells(cell) ])
    
    def unique_reads_ratio(self, cell):
        """Returns the number of single-end reads (pair-end reads count twice)."""
        return np.mean([ self.table.loc[c, 'total_deduplicated_percentage']/100.0 for c in self._allcells(cell) ])
    
    def get_full(self, cell):
        return [ self.table.loc[c] for c in self._allcells(cell) ][0]
    
    def print_stats(self):
        #TODO: add mappability stats
        #print(self.table.index)
        #display(table.head())
        return True