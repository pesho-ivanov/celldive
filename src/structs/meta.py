import pandas as pd
import json
from utils import *

class Meta:
    _default_color = '#000000'
    
    def __init__(self, meta):
        self.meta = meta
        '''Columns:
        * 'count' and 'count_int' (may be missing)
        * 'sample'
        * 'batch'
        * 'color'
        '''

    def dump(self, fn, batch_tag):
        to_json = {
            #batch_tag: self.meta.to_json(orient='index'),
            batch_tag: self.meta.to_dict('index'),
        }

        if Path(fn).suffix != 'json':
            Warning('{} is in json format without a json extension.'.format(fn))

        Path(fn).parent.mkdir(parents=True, exist_ok=True)
        with open(fn, 'w') as f:
            json.dump(to_json, f, indent=4)

        print('Metadata JSON dumped: {}'.format(Path(fn).name))
        
    @staticmethod
    def cells_int(cnt_str):
        if len(cnt_str) == 0 or not str.isdigit(cnt_str[0]):
            Warning('Cell count not a number: {}'.format(cnt_str))
            return 0
        else:
            return int(cnt_str[0])
    
    @classmethod
    def from_file(cls, input_file, *, tag):
        if not check_file(input_file, "  Metadata"):
            #return None
            raise FileNotFoundError("Metadata file %s should exist!" % blue(input_file))
        
        # Simple list of the microscope observations over the single-cell places
        column_types = {'well': str, 'count': str, 'sample': str, 'read_count': int}
        meta = pd.read_csv(input_file, delim_whitespace=True,
                           comment='#', dtype=column_types)
        meta.set_index('well', inplace=True, verify_integrity=True)
        meta.sort_index(inplace=True)
        if 'count' in meta.columns:
            meta['count_int'] = meta['count'].apply(Meta.cells_int)
        #sc = meta.loc[ (meta['count'] == '1')
        #                  | (meta['count'] == '1d')
        #                  | (meta['count'] == '1c')]
        ##sc = microscope_table.loc[(microscope_table['count'] != '0')]
        #print ("  -- %d wells with single cells (may contain debri): " % len(sc))

        #stats = get_stats(meta)
        #meta = pd.concat([stats, meta], axis=1)
        #return recombinants, meta
        #sample = get_sample(true_idx, batch3_which_sample)
        #meta['sample'] = sample['sample'] 
        
        meta['batch'] = tag  # set the tag for every cell
        add_prefix_tag(meta, tag)
        meta['color'] = cls._default_color
        
        ret = cls(meta=meta)
        return ret

    def well2count(self, well):
        if not well in self.meta.index:
            raise Exception('Well name {} is not found among the cells of the sample.'.format(well))
        return self.meta.loc[well, 'count_int']
    
    def set_color(self, color, cells=None):
        #color = set_alpha(color, 0.7)
        #color = add_more_alpha(color, 0.5)
        
        if cells != None:
            color = pd.Series(data=[color]*len(cells), index=cells)
            
        if type(color) is str:
            self.meta['color'] = color
        else:
            #assert type(color) is pd.Series
            self.meta.loc[color.index, 'color'] = color
            
    def get_colors(self):
        return self.meta['color']
    
    def get_color(self, cell):
        return self.meta.loc[cell, 'color']
    
    def get_color_rgb(self, cell):
        color = self.get_color(cell)
        return to_3rgb_floats(color)
        #if not color.startswith('#'):
        #    #if not color.startswith('xkcd:'):
        #    #    color = 'xkcd:' + color
        #    #mcd.XKCD_COLORS[color]
        #    #return color
        #    #return ColorConverter.to_rgb(color)
        #    return xkcd2rgb(color)
        #return color

    
    def submeta(self, cells):
        intersect = set(cells).intersection(self.cells())
        assert_same_indexes(cells, intersect)
        return Meta(self.meta.loc[list(intersect)])
    
    def concat(self, other):
        return Meta(pd.concat([self.meta, other.meta]).drop_duplicates())
        #renamed_metas = []
        #for table, rename_f in [(self.meta, self_rename),
        #                        (other.meta, other_rename)]:
        #    a = table.rename(index=rename_f)
        #    renamed_metas.append(a)
        #
        #return Meta(pd.concat(renamed_metas))
    
    def index(self):
        return self.meta.index
    
    def cells(self):
        return self.index().tolist()
    
    def print_stats(self):
        #display(self.meta)
        #assert 'color' in self.meta.columns
        
        cells = self.cells()
        #ret = self.meta.loc[cells[0], 'color']
        #print(cells)
        #print('ku: ', self.meta.loc['96', 'color'])
        colored_cells = cells
        #colored_cells = [ colorify(cell, self.meta.loc[cell, 'color']) for cell in cells ]
        #print(colored_cells)
        #print("    cells: %s" % iter2str(blue(self.cells())))
        maxn = 5
        print('    wells ({})'.format(blue(len(cells))))
        #print(': %s%s' % (', '.join(colored_cells[:maxn]), ", ..." if len(cells) > maxn else ''))
        samples = [ '\'%s\'(%s)' % (sample, blue(cnt)) for sample, cnt
                   in self.meta['sample'].value_counts().iteritems() ]
        print('    samples: ' + ', '.join(samples))
        try:
            counts = [ '\'%s\'(%s)' % (cell_counts, blue(cnt)) for cell_counts, cnt
                      in self.meta['count'].value_counts().sort_index().iteritems() ]
            print('    well counts: ' + ', '.join(counts)) 
        except:
            print('    well counts: %s' % red('???'))
    
    def cell2sample(self, cell):
        if not cell in self.meta.index:
            return None
        else:
            return self.meta.loc[cell, 'sample']
        
    def cell2batch(self, cell):
        if not cell in self.meta.index:
            return None
        else:
            return self.meta.loc[cell, 'batch']
    
    def get_index(self):
        return self.meta.index.copy()
    
