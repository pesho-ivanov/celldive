import subprocess# import PIPE, run#, call

from pathlib import Path
from structs.clone_struct import CellPart
from utils import *

class GSEA:
    """ http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    """
    def __init__(self, batch, *, out_dir, geneset):
        if not batch.cell_part:
            warning('No clonality information for batch {} => No GSEA analysis!'.format(batch.tag))
            return
        
        #clonal_cells = [ cell for cell in batch.cells() if batch.cell_part.is_clonal(cell) ]
        #nonclonal_cells = [ cell for cell in batch.cells() if batch.cell_part.is_nonclonal(cell) ]
        
        clonal_cells = list(batch.cell_part.part_dict[CellPart.CLONAL])
        nonclonal_cells = list(batch.cell_part.part_dict[CellPart.NONCLONAL])
        print("GSEA for {} clonal and {} nonclonal cells...".format(red(len(clonal_cells)), blue(len(nonclonal_cells))))
        if len(clonal_cells) < 5 or len(nonclonal_cells) < 5:
            warning('At least 5 clonal and 5 nonclonal cells are needed => No GSEA analysis!')
            return
        
        self.chosen_cells = clonal_cells + nonclonal_cells
        self.subbatch = batch.subbatch(self.chosen_cells)
        assert(set(self.subbatch.cells()) == set(self.chosen_cells))
        self.out_dir = out_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        
        self.gsea_jar = '/home/pesho/libs/gsea-3.0.jar'
        self.collapse_chip = 'gseaftp.broadinstitute.org://pub/gsea/annotations/ENSEMBL_human_transcript.chip' # TODO: download offline
        self.pheno_cls = out_dir / 'pheno.cls'
        self.write_pheno_cls(out_fn=self.pheno_cls)
        
        self.expressions = out_dir / 'expressions.txt'
        self.write_expressions(out_fn=self.expressions)
        
        potential_file = dirs['depend'] / 'msigdb_v6.0' / 'msigdb_v6.0_GMTs' / geneset  # geneset == 'c1.all.v6.0.symbols.gmt'
        if potential_file.is_file():
            self.geneset = potential_file
        else:
            self.geneset = geneset
        
        print('GSEA in {}'.format(blue(out_dir)))
        self.execute()
        
    def write_pheno_cls(self, *, out_fn):
        def unique(seq):
            seen = set()
            seen_add = seen.add
            return [x for x in seq if not (x in seen or seen_add(x))]
        
        clonality_arr = [ self.subbatch.cell_part.what_type_str(cell) for cell in self.chosen_cells ]
        assert all([ self.subbatch.cell_part.is_nonclonal(cell) | self.subbatch.cell_part.is_clonal(cell)
                    for cell in self.chosen_cells ]), 'Not all chosen cells are clonal or nonclonal'
        
        f = open(out_fn, "w") 
        print('Writing to {}'.format(blue(out_fn)))
        f.write("{} 2 1\n".format(len(self.chosen_cells))) 
        unique_arr = unique(clonality_arr)
        assert len(unique_arr) == 2
        f.write("# {}\n".format(' '.join(unique_arr)))
        f.write(' '.join(clonality_arr) + '\n')
        f.close() 
        
    def write_expressions(self, *, out_fn):
        df = self.subbatch.quant.ge.T
        df.insert(loc=0, column='descr', value=0)  # add an empty description column as a second column
        print('Writing to {}'.format(blue(out_fn)))
        df.to_csv(out_fn, sep='\t')
        
    def execute(self):
        cmd = ['java', '-cp', self.gsea_jar, '-Xmx2048m', 'xtools.gsea.Gsea',
               '-res', str(self.expressions),
               '-cls', str(self.pheno_cls),
               '-gmx', str(self.geneset),
               '-collapse true', '-mode Max_probe', '-norm meandiv', '-nperm 1000',
               '-permute phenotype', '-rnd_type no_balance', '-scoring_scheme weighted',
               '-rpt_label my_analysis', '-metric Signal2Noise', '-sort real',
               '-order descending',
               '-chip', str(self.collapse_chip),
               '-create_gcts false',
               '-create_svgs false', '-include_only_symbols true', '-make_sets true',
               '-median false', '-num 100', '-plot_top_x 20', '-rnd_seed timestamp',
               '-save_rnd_lists false', '-set_max 500', '-set_min 15', '-zip_report false',
               '-out', str(self.out_dir),
               '-gui false']
        print('Executing GSEA: ', green(' '.join(cmd)))
        #if subprocess.call(cmd):
        #    print(red("ERROR: GSEA failed!"))

        def filter_iterations(line):
            if not line.startswith('Iteration:'):
                yield line
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print(result.returncode, result.stdout, result.stderr)
        
