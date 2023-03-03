import random
from matplotlib.colors import to_rgb
from matplotlib.colors import to_hex
#from enum import Enum

from utils import *
from .genomics_convert import GenomicsConvert
from .quality_control import QualityControl
from .gene_expressions import GE
from .meta import Meta
from .tcrs import TCRs
from .vb_trbv import VbTrbv

#from algo.find_clones import FindClones
from .clone_struct import CellPartition, CellPart, CellColor, SampleColor

from algo.normalize import NormalizeByTotalSum, NormalizeBySpikeIn

class Batch:
 
    default_infiles = {
        #'expression_matrix': "expressions.tsv",
        'expression_matrix': "expressions.tsv",
        'quality_control': "multiqc_fastqc.txt",
        'TCRs': "filtered_TCRAB_summary/cd_tracer.json",
        #'TCRs': "filtered_TCRAB_summary/recombinants.txt",
        #'TCRs': "TCRs-tracer.txt",
        'meta': "meta.csv",
        #'vbs': ["Vb-distribution-CD4.tsv", "Vb-distribution-CD8.tsv", "Vb-distribution-CD4neg-CD8neg.tsv"],
        'vbs': ["Vb-distribution-CD3.tsv"],
    }

    def __init__(self, *, meta, quant, tcrs, qc, vbs, sample_id, internal_name, color,
                 diffexpr, cell_part, print_stats=False, out_dir=None):
        
        """Class for working with all data available for a batch of single cells.
        Usage: initialize using from_file
        """
        self.meta = meta
        """Meta information: cells in a well, cell samples, etc."""
        self.quant = quant
        """Gene/transcript expressions."""
        self.tcrs = tcrs
        """Table of TraCeR reconstructed TCRs."""
        self.qc = qc
        """Quality conrol class."""
        self.vbs = vbs
        self.tag = sample_id.strip()
        self.internal_name = internal_name.strip()
        self.color = color
        self.out_dir = out_dir
        
        self.diffexpr = diffexpr  # TODO: remove
        self.onlytag = None       # TODO: remove
        self.cell_part = cell_part
        
        self.do_color()
        if not self.is_valid():
            raise Exception('Batch not valid')
        
        if print_stats:
            self.print_stats()
    
    def do_color(self): 
        if not self.color:
            return
        
        #if self.color == 'random':
        #    self.color = random_color(seed=tag)

        if type(self.color) is str:
            color_dict = { '*': self.color }
        else:
            color_dict = self.color
            
        if not type(color_dict) is dict:
            raise Exception('color_dict is not a dict')
        self.select_and_set_color(color_dict)
            
        #if not self.clones is None:
        #    self.shade_cells(self.clonal_cells, 0.2)
        #    self.shade_cells(self.unknown_cells, 2.0)
        
    #def shade_cells(self, clonal_cells, coef):
    #    darker_color = [ min(0.9, coef * a) for a in to_rgb(xkcd2rgb(self.color))]
    #    self.meta.set_color(to_hex(darker_color), clonal_cells)
    
    #def copy(self, *, tag):
    #    if tag == self.tag:
    #        print(red('Warning') + ' copying a batch with the same tag %s' % blue(tag))
    #    return Batch(meta=self.meta, quant=self.quant, tcrs=self.tcrs, qc=self.qc, tag=tag)
    
    @classmethod
    def get_full_filenames(cls, batch_dir, input_files):
        if not set(input_files) <= set(cls.default_infiles):
            diff = set(input_files.keys()) - set(cls.default_infiles.keys())
            raise RuntimeError("The dict keys {%s} are not expected among {%s}"
                               % (', '.join(diff), ', '.join(cls.default_infiles.keys()))) 
        _ = input_files.copy()
        input_files = cls.default_infiles.copy()
        input_files.update(_)
        
        file_paths = to_full_filenames(batch_dir, input_files)
        return file_paths

    @staticmethod
    def __get_microscopic_filtered(*, tcrs, meta, filt_cells, tag):
        def should_filter(cell):
            if 'min_cells_in_well' in filt_cells:
                if count < filt_cells['min_cells_in_well']:
                    return True
            if 'max_cells_in_well' in filt_cells:
                if count > filt_cells['max_cells_in_well']:
                    return True
            return False

        filtered = []
        for cell in tcrs.cells():
            count = meta.well2count(cell)
            if should_filter(cell):
                filtered.append(cell)

        return filtered
    
    @classmethod
    def from_file(cls, batch_dir, *, sample_id, internal_name, color, diffexpr, input_files={}, clonality, kallisto_dir, #conf_de,
                  clonotype, expression_threshold, min_cells_expressing_gene, filt_cells, filt_tcr_chains, norm, sep='\t', nrows=None, out_dir=None):
        """input_files should be a dict with filenames for (expression_matrix, quality_control, TCRs, meta) """
        c = color if type(color) is str else 'black'
        id_with_color = xkcd2console(sample_id if sample_id else "", c)
        #print(xkcd2console(" -----------------------------", c))
        print("{}, ".format(id_with_color), end='')
        #print("Loading batch {} ({}) from directory {}...".format(
        #      id_with_color, green(onlytag), blue(batch_dir)))
        file_paths = cls.get_full_filenames(batch_dir, input_files)
        
        meta = Meta.from_file(file_paths['meta'], tag=sample_id)
        true_idx = meta.get_index()
        tcrs = TCRs.from_file(true_idx, file_paths['TCRs'], sample_id=sample_id, filt_tcr_chains=filt_tcr_chains)
        microscopic_filtered = cls.__get_microscopic_filtered(tcrs=tcrs, meta=meta, filt_cells=filt_cells, tag=sample_id)
        qc = QualityControl.from_file(true_idx, file_paths['quality_control'], tag=sample_id)

### BEGINNING OF WORKAROUND
        #quant = None
        quant, orig_index, filtered = GE.from_file(
            true_idx,
            file_paths['expression_matrix'],
            kallisto_dir=kallisto_dir,
            #transform='log2', 
            sample_id=sample_id,
            #expression_threshold=expression_threshold,
            #min_cells_expressing_gene=min_cells_expressing_gene, 
            filt_cells=filt_cells,
            nrows=nrows)

        print('Filtering out {} cells in sample {}.'.format(len(filtered), sample_id))
        if quant is None:
            warning('No gene expression provided for sample {}'.format(sample_id))
        #quant = None
#### END OF WORKAROUND

        #print("batch.py: ")
        #print(quant.ge.sum(axis=1))
        filtered += microscopic_filtered

        if tcrs and not (set(meta.cells()) & set(tcrs.cells())):
            Warning('No common cells between the meta and tcrs. Most probably the format is different.')
            print('example meta cells: ', meta.cells()[:3])
            print('example tcrs cells: ', tcrs.cells()[:3])

        # TODO: move away
        vbs = {}
        vbs['TraCeR'] = VbTrbv.from_tcrs(tcrs, 'tcrs')
        for fn in file_paths['vbs']:
            vb_trbv = VbTrbv.from_file(fn)
            if vb_trbv:
                vbs[ str(fn.stem) ] = vb_trbv
        
        if clonality:
            if 'bystander_clone_threshold_percent' not in clonality:
                raise Exception('clonality->bystander_clone_threshold_percent number should also be provided.')
            bystander_clone_threshold_percent = clonality['bystander_clone_threshold_percent']
            cell_part = CellPartition.TC_with_threshold(
                    sample_id=sample_id,
                    cell_files=orig_index,
                    all_cells=meta.cells(), tcrs=tcrs,
                    clonotype=clonotype, filtered=filtered, 
                    bystander_clone_threshold_percent=bystander_clone_threshold_percent)
            #if 'method' not in clonality:
            #    raise Exception('clonality->method number should also be provided.')
            #method = clonality['method']
            #if method == 'old':
            #    cell_part = CellPartition.from_clonality_rule(all_cells=meta.cells(), tcrs=tcrs,
            #            clonotype=clonotype, filtered=filtered,
            #            bystander_clone_threshold_percent=bystander_clone_threshold_percent)
            #if method == 'tc':
            #else:
            #    raise Exception('Unknown clonality_method: {}'.format(method))
        
        batch = cls(meta=meta, quant=quant, tcrs=tcrs, qc=qc, vbs=vbs, sample_id=sample_id, cell_part=cell_part,
                    internal_name=internal_name, color=color, diffexpr=diffexpr, out_dir=out_dir)
        #print("batch.py::from_file batch.tcrs=", batch.tcrs)
        #batch.onlytag = onlytag

        return batch
#        #ret = batch.sieved(filt_cells)
#
#        #ret.get_clones_lazy(onlytag=onlytag, diffexpr=diffexpr)
#        #print("batch.py::from_file sieved tcrs=", ret.tcrs)
#        
#        if norm == 'byTotalSum':
#            ret.quant.ge = NormalizeByTotalSum(ret.quant.ge)
#        elif norm == 'bySpikeIn':
#            ret.quant.ge = NormalizeBySpikeIn(ret.quant.ge)
#        else:
#            assert not norm
#        #quant = Normalize(quant) 
#        
#        return ret
    
    @staticmethod
    def transform_query(expr, obj, silent=False):
        if not silent:
            print("Selecting the cells by: %s" % blue(expr))
        
        if expr == '' or expr == '*':
            expr = "self.meta.meta.index"
            
        #if len(obj) > 0 and obj[-1] != '.':
        #    obj += '.'
            
        subst = { 'cells_in_well': obj + "['count_int']",
                         'sample': obj + "['sample']",
                          'reads': obj + "['Total Sequences']",
                         'mapped': obj + "['mapped']",
                              'M': "000000",
                              'K': "000",
        }
        
        for old, new in subst.items():
            expr = expr.replace(old, new)
        #if not silent:
        #    print('adapted: ' + blue(expr))
            
        return expr
    
    def select_cells(self, row_expr, silent=False):
        adapted_row_expr = self.transform_query(row_expr, 'merged', silent)
        
        #merged = pd.concat([self.meta.meta, self.qc.table], axis='columns')
        if self.qc:
            merged = pd.merge(self.meta.meta, self.qc.table,
                              left_index=True, right_index=True, how='left')
        else:
            merged = self.meta.meta
            
        if self.quant is not None:
            merged['mapped'] = self.quant.ge.sum(axis='columns')
        
        try:
            submeta = merged.loc[ eval(adapted_row_expr) ]
        except:
            print(red("Cannot parse."))
            submeta = self.meta.meta.drop(self.meta.meta.index)
        return submeta.index.tolist()
    
    def is_valid(self):
        meta_cells = set(self.meta.cells())
        if self.tcrs:  assert_same_indexes(meta_cells, self.tcrs.cells(), tag='tcrs')
        if self.quant: assert_same_indexes(meta_cells, self.quant.cells(), tag='quant')
        #if self.qc:
        
        return True
    
    def cells(self):
        if not self.is_valid():
            raise Exception("Different cells")
        #return self.tcrs.cells()

        return self.meta.cells()
    
    def genes(self):
        if self.quant is None:
            raise Warning('No quantification for batch {}'.format(self.tag))
        return self.quant.transcripts()
    
    def index(self):
        return self.meta.index()
    
    def set_color(self, cell2color_series):
        self.meta.set_color(cell2color_series)
        
    def _get_clonal_color(self, cell, annot=False):
        if not self.cell_part:
            part = 'unknown'
        else:
            part = self.cell_part.get_cell_part(cell)

        return part if annot else CellColor[part]

            #except AttributeError:
            #else: raise Exception("Not existing cell {} in the batch {}".format(cell, self.tag))

        #if not self.cell_part or self.cell_part.is_unknown(cell):  #cell in self.unknown_cells:
        #    if annot: return 'unknown'
        #    return self.ClonalityColor.UNKNOWN
        #elif self.cell_part.is_clonal(cell):
        #    if annot: return 'clonal'
        #    return self.ClonalityColor.CLONAL
        #elif self.cell_part.is_nonclonal(cell):  #cell in self.nonclonal_cells:
        #    if annot: return 'nonclonal'
        #    return self.ClonalityColor.NONCLONAL
        #elif self.cell_part.is_nonreconstructed(cell):  #cell in self.unknown_cells:
        #    if annot: return 'nonreconstructed'
        #    return self.ClonalityColor.NONRECONSTRUCTED
        #elif self.cell_part.is_filtered(cell):  #cell in self.unknown_cells:
        #    if annot: return 'filtered'
        #    return self.ClonalityColor.FILTERED
        #else:
        #    raise Exception("Not existing cell {} in the batch {}".format(cell, self.tag))
        
    def get_color(self, cell, *, annot=False, coloring):
        cell_color = self.meta.get_color(cell)
        
        if coloring == SampleColor.SAMPLE:
            if annot: return ''
            return cell_color
        elif coloring == SampleColor.CLONAL:
            return self._get_clonal_color(cell, annot)
        elif coloring == SampleColor.SAMPLE_AND_CLONAL:
            if self.cell_part.is_clonal(cell):  #cell in self.clonal_cells:
                if annot: return 'clonal'
                coef = 0.2
            elif self.cell_part.is_unknown(cell):
                if annot: return 'unknown'
                coef = 2.0
            elif self.cell_part.is_nonreconstructed(cell):
                if annot: return 'nonreconstructed'
                coef = 2.0
            else:
                if annot: return 'nonclonal'
                return cell_color
            
            coef = 1.0
            rgb3 = to_rgb(to_3rgb_floats(xkcd2rgb(cell_color)))
            shifted_color = [min(0.9, coef * channel) for channel in rgb3]
            return to_hex(shifted_color)
        else:
            raise Exception('No such coloring: {}'.format(coloring))
        
    def select_and_set_color(self, expr2color: "expr -> color"):
        all_cells = set()
        for expr, color in expr2color.items():
            cells = self.select_cells(expr, silent=True)
            intersection = all_cells & set(cells)
            if intersection:
                print(red('WARNING: ') + 'Some cells are not uniquely colored: ' + iter2str(blue(intersection)))
            all_cells.update(cells)
            
            if not cells:
                print(red('WARNING: ') + 'No cells are selected to be colored as ' + color)
            self.meta.set_color(color, cells)
        if len(all_cells) != len(self.cells()):
            print(red('Warning: ') + 'Colored {} cells out of {} in total.'.format(len(all_cells), len(self.cells())))
        
    def set_tag(self, tag):
        self.tag = tag
    
    def head(self, n=10, *, tag):
        return self.subbatch(self.cells()[:n], tag=tag)
   
    def sample(self, n=10, *, tag):
        indices = random.sample(range(len(self.cells())), n)
        cells = [self.cells()[i] for i in sorted(indices)]
        return self.subbatch(cells, tag=tag)

    def get_genes_from_file(self, genes_file):
        print('Loading genes from {}'.format(blue(genes_file)))
        genes_df = pd.read_csv(genes_file, index_col=0)
        return genes_df.index.tolist()
        
    def subbatch(self, cells='all', *, tag=None, features=None, color=False, print_stats=False):
        raise Exception("No subbatching!")
        if type(cells) is str:
            cells = self.select_cells(cells)
            
        if not cells:
            print(red('Warning: subbatch with no cells'))
            
        if not color:
            color = self.color
        if not tag:
            tag = self.tag
        if not features:
            features = 'all'
            
        assert_same_indexes(self.meta.cells(), cells, tag='subbatch', subset_ok=True) # TODO: tag?
        #if not set(cells) <= set(self.meta.cells()):
        #    raise Warning('')
        if self.quant is None:
            print('Warning: quant not existing')
            subquant = None
        elif features == 'all':
            subquant = self.quant.subge(cells)
        elif features == 'important':
            #print(self.quant.ge.head())
            genes = self.get_important_genes()
            subquant = self.quant.subge(cells, genes)
        elif type(features) is int:
            genes = self.get_important_genes()
            if len(genes) < features:
                print('Warning: Not enough genes: {} / {}'.format(len(genes), features))
                subquant = self.quant.subge(cells, genes)
            else:
                # TODO!!!!!!!!!!!!!!!!!! 
                #subquant = self.quant.subge(cells, genes[:features])
                #indices = random.sample(range(len(genes)), features)
                #random_genes = [ genes[i] for i in sorted(indices) ]
                subquant = self.quant.subge(cells, random_genes)
        elif Path(features).is_file():
            genes = self.get_genes_from_file(features)
            subquant = self.quant.subge(cells, genes)
        elif Path(features).is_dir():  ## unite genes from all files in the dir
            files = list(Path(features).glob('**/genes_*.csv'))
            if not files:
                warning('No genes files!')
            all_genes = set()
            for f in files:
                genes = self.get_genes_from_file(f)
                all_genes |= set(genes)
            subquant = self.quant.subge(cells, all_genes)
        else:
            assert False, 'which genes to use?! given features={}'.format(features)
            #subquant = self.quant.subge(cells, features)
            
        submeta   = self.meta.submeta(cells)        if self.meta else None
        subtcrs   = self.tcrs.subtcrs(cells=cells)  if self.tcrs else None
        subqc     = self.qc.subqc(cells)            if self.qc else None
        subvbs    = self.vbs
        subcell_part = self.cell_part & cells       if self.cell_part else None
        if not subvbs is None:
            subvbs['TraCeR'] = VbTrbv.from_tcrs(subtcrs, 'tcrs')
        
        subbatch = Batch(meta=submeta, quant=subquant, tcrs=subtcrs, qc=subqc, vbs=subvbs, sample_id=self.tag,
                         internal_name=self.internal_name, color=color, diffexpr=self.diffexpr,
                         cell_part=subcell_part, print_stats=print_stats, out_dir=self.out_dir)
        #if self.clones is not None:
        #    subclones = self.clones.subclones(subbatch, cells)
        #    #subbatch.get_clones_lazy(clones=subclones, onlytag=subclones.onlytag, diffexpr=self.diffexpr)
        #    subbatch.get_clones_lazy(clones=subclones, onlytag=self.onlytag, diffexpr=self.diffexpr)
        #else:
        #    subclones = None
        
        return subbatch
    
    def sieved(self, filt):
        if not 'cells' in filt:
            filt['cells'] = '*'
        if not 'genes' in filt:
            filt['genes'] = 'all'
        
        return self.subbatch(cells=filt['cells'], features=filt['genes'])

    
    #def __add__(self, other):
    #    sect = set(self.cells()).intersection(other.cells())
    #    assert len(sect) == 0, red("WARNING:") + " Both batches contain the cells %s" % iter2str(blue(sect))
    #    
    #    summeta = self.meta.concat(other.meta)
    #    sumquant = self.quant.concat(other.quant)
    #    sumtcrs = self.tcrs.concat(other.tcrs)
    #    sumqc = self.qc
    #            
    #    return Batch(meta=summeta, quant=sumquant, tcrs=sumtcrs, qc=sumqc, tag='', color=False, print_stats=False)

    @classmethod
    def merge(cls, a, b, *other, sample_id, internal_name, color=False, cell_part):
        #raise Exception('No merging!')

        #if a.tag == b.tag:
        #    print(red('Warning') + ' concatinating batches with same %s tags' % blue(a.tag))
        #    assert False
        
        sect = set(a.cells()) & set(b.cells())
        #assert len(sect) == 0, red("WARNING:") + " Both batches contain the cells " + iter2str(blue(sect))
        warning(" Both batches contain the cells " + iter2str(blue(sect)))
        
        summeta = a.meta.concat(b.meta)
        sumquant = a.quant.concat(b.quant, sample_id) if a.quant is not None else None
        sumtcrs = a.tcrs.concat(b.tcrs, sample_id)
        intersection = { key: value for key, value in sumtcrs._cell2tcrs.items() if key in cell_part.cells() }
        sumtcrs = TCRs(intersection, sample_id)
        sumqc = a.qc
        #assert set(a.vbs.columns) == set(a.vbs.columns), 'Error: different columns'
        assert a.diffexpr == b.diffexpr
        sumvbs = None   # TODO

        
        ab = cls(meta=summeta, quant=sumquant, tcrs=sumtcrs, qc=sumqc, sample_id=sample_id, internal_name=internal_name,
                 diffexpr=a.diffexpr, vbs=sumvbs, color=color, print_stats=(len(other)==0), cell_part=cell_part)
        #ab.get_clones_lazy(onlytag=self.onlytag, diffexpr=self.diffexpr)
        
        if len(other) == 0:
            return ab
        elif len(other) == 1:
            return cls.merge(ab, other[0], tag=sample_id, color=color)
        else:
            cde = cls.merge(*other, tag=sample_id, color=color)
            return cls.merge(ab, cde, tag=sample_id, color=color)
        
    #def find_clones(self):
    #    if not self.clones:
    #        self.get_clones_lazy(onlytag=self.onlytag)
    #    return self.clones
    
    def get_onlytag(self):
        if self.cell_part:
            return '_'.join(self.cell_part.clonotype)
        else:
            return ''
        #return '_'.join(self.cell_part.clonotype) if self.cell_part is not None else ''
    
    def get_important_genes2pvals(self):
        if self.clones is not None and self.clones.best_hypo():
            return self.clones.best_hypo().best_genes
            #g = self.clones.best_hypo().best_genes
            #return g[g < 0.05]
        else:
            return None
    
    def get_important_genes(self):
        genes2pvals = self.get_important_genes2pvals()
        if genes2pvals is not None:
            return genes2pvals.index.tolist()
        else:
            return []
    
    #def get_union_important_genes2pval(self):
    #    return self.get_clones_lazy().union_important_genes
    
    #def get_union_important_genes(self):
    #    return self.get_important_genes2pval().index.tolist()
    
    #def genes2quant(self, genes):
    #    #if self.important_quant.empty:
    #    #genes = self.get_union_important_genes()
    #    if genes is None:
    #        return self.quant.ge.drop(self.quant.ge.index)
    #    elif genes == 'important':
    #        return self.quant.ge[ self.get_important_genes() ]
    #    else:
    #        return self.quant.ge[ genes ]
        
    def __repr__(self):
        #clf = HypoClassifier(self.clones.best_hypo(), kernel='linear', C=1.0)
        #if not clf.scores is None:
        #    acc = "%0.2f (+/- %0.2f)" % (clf.scores.mean(), clf.scores.std() * 2)
        #else:
        #    acc = red('NOT TRAINED')
        res = ''
        if type(self.color) is str:
            tag = xkcd2console(self.tag, self.color)
        else:
            tag = self.tag
        return tag

        res += '{}: {} cells\n'.format(tag, blue(len(self.cells())))
        #if not self.clones is None:
        #    res += self.clones.__repr__()
        return res
    
    def print_stats(self):
        print(" --= batch %s =--" % (self.tag))
        if (self.meta):  self.meta.print_stats()
        if (self.tcrs):  self.tcrs.print_stats()
        if (self.quant): self.quant.print_stats()
        if (self.qc):    self.qc.print_stats()
        #meta.print_stats()
        #masks = pd.DataFrame({
        #    'TCRs': [len(x) for x in meta.recombinants],
        #}, index=meta.index)
        #print (sum(self.TCRs > 0), 'single cells with at least one TCR reconstructed')
        #print (sum(masks.cells == 1), 'single cells out of', len(masks.cells), 'wells')
