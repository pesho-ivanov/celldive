from file_utils import *

import pandas as pd

class GenomicsConvert:
    def __init__(self, table):
        """Class for converting gene ids to gene names to transcript ids to transcript names.
        Downloaded from Ensembl. TODO: automate download
        Initialize using fromTable or fromDir.
        """
        self.table = table
        
        def lists2dict(keys, values):
            d = {}
            for k, v in zip(keys, values):
                d.setdefault(k, []).append(v)
            return d
        
        self.gname2tids_dict = lists2dict(table['geneName'], table['transcriptId'])
        #self.gname2tids_dict = dict(zip(table['geneName'], table['transcriptId']))
        self.tid2tname_dict = dict(zip(table['transcriptId'], table['transcriptName']))
        """A pandas table with a row for a transcript entry."""
        assert len(self.table) > 0
        
    @classmethod
    def from_table(cls, table):
        """          Gene ID,        Transcript ID,   Associated Gene Name, Associated Transcript Name
            Example: ENSG00000252303,ENST00000516494, RNU6-280P,            RNU6-280P-202
        """
        return cls(table)
    
    @classmethod
    def from_file(cls, dirorfile, infile='id2name.csv', sep=','):
        dirorfile = Path(dirorfile)
        assert_exists(dirorfile)
        if dirorfile.is_dir():
            fn = dirorfile / infile
        else:
            fn = dirorfile
        assert_file_exists(fn)
        columns = ['geneId', 'transcriptId', 'geneName', 'transcriptName']
        table = pd.read_csv(fn, sep=sep, header=0, names=columns)
        return cls(table)
        
    def gid2gname(self, gid):
        """geneId -> gene name"""
        row = self.table.loc[self.table['geneId'] == gid]
        # TODO raise an exception
        if not len(row) == 1:
            Warning(str(len(row)) + ' gene names found for the geneId: ' + str(gid))
        return row['geneName'].tolist()
    
    def gname2gid(self, gname):
        """gene name -> geneId"""
        rows = self.table.loc[self.table['geneName'] == gname]
        return rows['geneId'].tolist()
    
    def tid2tname(self, tid):
        """transcriptId -> transcript name"""
        if tid in self.tid2tname_dict:
            return self.tid2tname_dict[tid]
        else:
            return 'NO_NAME'
        
        #row = self.table.loc[self.table['transcriptId'] == tid]
        ## TODO raise an exception
        #if not len(row) == 1:
        #    #print(red('Warning: ') + str(len(row)) + ' transcript names found for the transcriptId ' + str(tid))
        #    return 'NO_NAME'
        #return row['transcriptName'].tolist()[0]
    
    def gname2tids(self, tname):
        """transcriptId -> gene name"""
        if tname in self.gname2tids_dict:
            return self.gname2tids_dict[tname]
        else:
            return 'NO_GENE'
        #row = self.table.loc[self.table['transcriptName'] == tname]
        #if not len(row) == 1: raise Warning(str(len(row)) + ' transcript IDs found for the transcript name ' + str(tname))
        #return row['transcriptId'].tolist()[0]
    
    def gene_names(self):
        return sorted(self.table['geneName'].unique())
    
    def transcript_names(self):
        return sorted(self.table['transcriptName'].unique())