import pandas as pd
from collections import defaultdict

class VbTrbvHelper:
    def __init__(self, vb2trbv_df):
        self._vb2trbvs = defaultdict(list)
        self._trbv2vbs = defaultdict(list)
        for row in vb2trbv_df.iterrows():
            vb, trbvs = row[1]['vb'], row[1]['trbv'].split(',')
            self._vb2trbvs[vb] = trbvs
            for trbv in trbvs:
                self._trbv2vbs[trbv].append(vb)
        
    @classmethod
    def init_from_file(cls, vb2trbv_fn):
        '''Read the correspondence table.'''
        if not vb2trbv_fn.is_file() :
            raise FileNotFoundError("Vb file %s should exist!" % vb2trbv_fn)
        column_types = {'vb': str, 'trbv': str}
        vb2trbv_df = pd.read_csv(vb2trbv_fn, delim_whitespace=True, comment='#',
                                 dtype=column_types)
        return cls(vb2trbv_df)
    
    def vbs(self):
        return self._vb2trbvs.keys()
    
    def trbvs(self):
        return self._trbv2vbs.keys()
    
    def vb2trbvs(self, vb):
        return self._vb2trbvs[vb]
    
    def trbv2bvs(self, trbv):
        return self._trbv2bvs[trbv]
    