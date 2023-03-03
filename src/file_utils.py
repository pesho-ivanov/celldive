from pathlib import Path
import numpy as np

class Dirs:
    def __init__(self, base):
        self.base = Path(base).resolve()
        assert_dir_exists(self.base)
        self.inp = self.base / 'in/'
        assert_dir_exists(self.base)
        #print("Batches: %s" % ', '.join(self.list_batches()))

    def __getitem__(self, subdir):
        res = self.base / subdir
        #assert_exists(res)
        return res
    
    #def list_batches(self):
    #    return np.sort([str(x.name) for x in self.inp.iterdir() \
    #                    if x.is_dir() and not str(x.name).startswith('.')])
        
def assert_exists(dorf):
    if not Path(dorf).exists(): raise FileNotFoundError(str(dorf) + ' not found')
        
def assert_dir_exists(d):
    if not Path(d).is_dir(): raise FileNotFoundError('The dir ' + str(d) + ' not found')
        
def assert_file_exists(f):
    if not Path(f).is_file(): raise FileNotFoundError('The file ' + str(f) + ' not found')

def to_full_filenames(directory, fns):
    d = Path(directory)
    #assert_dir_exists(d)

    def fn2full(fn):
        if type(fn) is list:
            return [d/f for f in fn]
        else:
            return d/fn
    file_paths = {source: fn2full(fn) for source, fn in fns.items()}

    #for source, fn in file_paths.items():
    #    assert_file_exists(fn)

    #if not fn.is_file():
    #    if fn.suffix == '.gz' or fn.suffix == '.gzip':
    #        fn = fn.with_suffix('')
    #    else:
    #        fn = fn.with_suffix('')

    return file_paths
