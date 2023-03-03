from IPython.display import display
from pprint import pprint
from datetime import datetime
import time
from pathlib import Path
import matplotlib.pyplot as plt
import itertools
import subprocess
import sys
import json
from collections import OrderedDict, defaultdict

from structs.batch import Batch

from vis.batch_stats import BatchStats
from vis.ge_scatter_tcr_graph import ScatterTCR
from vis.bipartite_tcrs import BipartiteTCR
from vis.vis_force_graph import draw_TCR_graph
from vis.TCR_heatmap import TCRheatmap
#from vis.quant_heatmap import QuantHeatmap

from structs.vb_trbv import VbTrbv
from structs.clone_struct import SampleColor, CellPart

from algo.filter_batch import FilterBatch
#from algo.hypo_classifier import HypoClassifier
#from algo.grid_cross_validation import GridCrossValidation
from algo.gene_intersect import GeneIntersect
from algo.clonality_stats import ClonalityStats
from algo.dump_clonality import DumpClonality
from algo.stats_tcr import StatsTCR
#from algo.gsea import GSEA

from multiprocessing_logging import install_mp_handler
from multiprocessing import Pool
from multiprocessing.managers import BaseManager, DictProxy

# Logging to both screen and logfile
class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()

class MyManager(BaseManager):
    pass

MyManager.register('defaultdict', defaultdict, DictProxy)
mgr = MyManager()
mgr.start()

#loggig.basicConfig(...)
install_mp_handler()

#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.INFO) #multiprocessing.SUBDEBUG)

from utils import *

EXT = 'png' #'svg' #'png'

def indict(key, d):
    return isinstance(d, collections.Mapping) and key in d and d[key]

class Run:
    config_dumped_file     = 'config.json'
    sample_description_dump_file = 'sample_descriptions.tsv'
    execution_log_file_tmpl1     = 'execution_{}.log'

    cv_file_tmpl           = 'cv-table_B_{}.html'
    bipartite_file_tmpl    = 'bipartite/bipartite_{}' + '.' + EXT
    scatter_file_tmpl      = 'scatter/scatter_{}' + '.' + EXT
#   quantheatmap_file_tmpl = 'heatmap/quantheatmap_{}' + '.' + EXT
    vb_kit_distr_file_tmpl     = 'vb_kit_distr/vb_kit_distr_{}' + '.' + EXT
    vj_heatmap_file_tmpl   = 'vj_heatmap/vj_heatmap_{}' + '.' + EXT
#    gsea_folder_tmpl       = 'gsea_{}/'
#   genes_file_tmpl        = 'genes/genes_{}.csv'
    #de_dir                 = 'DE/'
    de_config_file_tmpl4    = 'DE_genes/genes_{}/DE_config_{}_{}_{}.json'
    de_file_tmpl4           = 'DE_genes/genes_{}/DE_genes_{}_{}_{}.csv'
    gene_table_file_tmpl2  = 'DE_gene_table/genes_{}/DE_gene_table_{}_{}.csv'
    gene_stats_file_tmpl   = 'DE_stats/genes_{}/DE_gene_stats_{}.csv'
    de_sample_by_samples_table = 'DE_sample_by_samples_tables/genes_{}/DE_samples_table_{}_vs_{}_{}.csv'
    upset_file             = 'upset' + '.' + EXT
    clonality_file_tmpl    = 'clonality/clonality_{}.json'
    meta_json_tmpl         = 'meta/meta_{}.json'
#   meta_stats_file        = 'meta_stats.csv'
    clonality_stats_file   = 'clonality_stats.csv'
#   clonality_stats_counts_file_tmpl = 'clonality/clonality_stats_counts_{}.csv'
#   clonality_stats_cells_file_tmpl = 'clonality/clonality_stats_cells_{}.json'
    batch_stats_file_tmpl  = 'stats/stats_sample_{}' + '.' + EXT
    stats_tcr_file_tmpl    = 'stats/stats_tcr_{}' + '.' + EXT
#   log_file               = 'log.txt'
    last_d                 = 'out_latest'

    def __init__(self, dirs):
        self.dirs = dirs
        self.B = {}  # batches 
        self.samples = []
        
        # cv params
        #self.Cs = []
        #self.kernels = []
        #self.gammas = []
        #self.features = []
        
        self.unsuccessful = mgr.defaultdict(list)

        self.bipartite_files = []
        #self.gcv = None
        
        self.samples_description = None  # pd.df
        
    def set_output_dir(self, *, diffexpr=None, nrows, filt): 
        dt = datetime.now()
        time = dt.strftime("%Y%m%d_%Hh%Mm%Ss")
        parts = []
        parts.append(time)
        #if self.samples:
        #    parts.append('samples={}'.format('+'.join([patient + '_'.join([s])
        #                                                 for patient, s, *_ in self.samples])))
        #if diffexpr:
        #    parts.append('DE.{}'.format(','.join(diffexpr)))
        #if nrows:
        #    parts.append('limit_genes.{}'.format(nrows))
        ##print('filt: {}'.format(blue(str(filt).translate(str.maketrans("'{}:", '    ')).replace(" ", ""))))
        #if filt:
        #    parts.append('filt.{}'.format(str(filt).translate(str.maketrans("'{}:=", '     ')).replace(" ", "")))
        #conf_str = str(self.conf).replace(' ', '').replace('{', '').replace('}', '').replace("'", '').replace('"', '').replace(':', '-').replace('/', '|')
        #conf_str = '_'.join("{!s}={!r}".format(key,val) for (key,val) in self.conf.items())
        #parts.append(conf_str)
        #print('running.py, conf_str: ', conf_str)
        self.folder_tag = '_'.join(parts)
        #parts.append('experiment.{}'.format(self.experiment))

        #conf_str = json.dumps(self.conf)
        #conf_str = str(self.conf)
        
        #dirs['out'].mkdir(parents=True, exist_ok=True)
        self.out_dir = dirs['out'] / '_'.join(parts)
        #self.out_dir = Path('/home/pesho/tmp/gsea.test12')  # TODO: remove
        self.out_dir.mkdir()
        latest_dir = dirs.base / self.last_d
        if Path(latest_dir).exists():
            latest_dir.unlink()
        latest_dir.symlink_to(self.out_dir, target_is_directory=True)
        
        print('Output dir: {}'.format(blue(str(self.out_dir))))
        
        ## redirect output to file
        #oldstdout = sys.stdout
        #sys.stdout = open(self.out_dir/self.log_file, 'w')
        #sys.stdout = oldstdout
        
    #def print_batches(self):
    #    # TODO colors are tagged for console
    #    txt = ''
    #    for b in self.B.values():
    #        print(b)
    #        txt += b.__repr__()
    #        
    #    if self.out_dir:
    #        tag = 'features_B={}.txt'.format('-'.join(self.B.keys()))
    #        fn = self.out_dir / tag
    #        print('Writing cross-validation results to {}'.format(fn))
    #        txt_file = open(fn, "w")
    #        txt_file.write(txt)
    #        txt_file.close()
        
    # TODO remove features
    def draw_scatter(self, batch, *, dim_reduct, perplexity=30, features, interactive):
        #if features=='important':
        #    features = get_important_features()
        #    #features = self.find_clones().hypothesss[0].best_features
        #else:
        #    features = None
            
        # colored by sample, colored by clonality, colored by both
        
        #points2d = None
        #assert 'coloring' in self.conf
        
        #for coloring in SampleColor:
        #for coloring in self.conf['coloring']:
        for coloring in self.conf['scatter']['coloring']:
            assert coloring in SampleColor.__dict__.values()
            if not coloring == SampleColor.SAMPLE and not batch.cell_part:
                continue
            get_color = lambda cell, annot=False: batch.get_color(cell, coloring=coloring, annot=annot)
            f = 'gene_list' if Path(features).exists() else str(features)
            out_fn = str(self.scatter_file_tmpl).format(
                '_'.join([batch.tag, coloring, dim_reduct, str(batch.onlytag), f]))
            scatter_tcr = ScatterTCR(batch, dim_reduct=dim_reduct, perplexity=perplexity,
                                     interactive=interactive, 
                                     out_fn=str(self.out_dir/out_fn), get_color=get_color,
                                    ) #points2d=points2d)
            #points2d = scatter_tcr.points2d
            
#    def draw_quantheatmap(self, batch, *, features):
#        for coloring in self.conf['quantheatmap']['coloring']:
#            assert coloring in SampleColor.__dict__.values()
#            #if not coloring == SampleColor.SAMPLE and not batch.cell_part:
#            #    continue
#            get_color = lambda cell, annot=False: batch.get_color(cell, coloring=coloring, annot=annot)
#            
#            features = str(features)
#            f = features if '/' not in features else features.split('/')[-1]
#            specific = '{}_{}({})_{}_{}'.format(batch.tag, f, len(batch.genes()), batch.get_onlytag(), coloring)
#            out_fn = self.out_dir / str(self.quantheatmap_file_tmpl).format(specific)
#            out_fn.parent.mkdir(parents=True, exist_ok=True)
#            #get_color = lambda cell, annot=False: batch.get_color(cell, coloring='byClone', annot=annot)
#            QuantHeatmap.draw(batch, out_fn=out_fn, get_color=get_color)
            
    def exec_per_batch(self, batch, unsuccessful, samples):
        conf = self.conf
        status = {}
        OK = 'OK'
        ERROR = 'ERROR'
        SKIPPED = 'SKIPPED'

        interactive = False # conf['interactive']
        #features = conf['features']
            
        if True: #indict('meta', conf):
            meta_json_file = self.out_dir / str(self.meta_json_tmpl).format(batch.tag)
            batch.meta.dump(meta_json_file, batch.tag)
            status['dump_meta'] = OK

        if indict('stats', conf) and indict('samples', conf['stats']):
            out_fn = str(self.batch_stats_file_tmpl).format(batch.tag)
            batch_stats = BatchStats(batch)
            batch_stats.draw_all(out_fn=self.out_dir/out_fn, interactive=interactive,
                                 filt_genes = conf['filt_genes'], filt_cells = conf['filt_cells'])
            status['stats_batch'] = OK
            
        if indict('stats', conf): # and indict('tcrs', conf['stats']):
            out_fn = str(self.stats_tcr_file_tmpl).format(batch.tag)
            print('Writing to {}'.format(blue(out_fn)))
            stats_tcr = StatsTCR(batch)
            stats_tcr.plot_edit_distance_histograms(out_fn=self.out_dir/out_fn)
            status['edit_distance_histograms'] = OK
            
        if indict('clonality', conf):
            out_fn = self.out_dir / str(self.clonality_file_tmpl).format(batch.tag)
            kallisto_dir = conf['kallisto_dir_tmpl'].format(batch.tag)
            clonality_dumper = DumpClonality(batch, kallisto_dir)
            clonality_dumper.dump(out_fn)
            status['dump_clonality'] = OK
                
        # TODO: the flow doesn't go after this if
        if indict('bipartite', conf):
            #use_cellcolors = conf['bipartite_cellcolors']
            for use_cellcolors in conf['bipartite']: # ['clonal', 'samples', 'rainbow']:
                #if use_cellcolors == 'clonal' and not batch.cell_part:
                if use_cellcolors == SampleColor.CLONAL and batch.onlytag == 'NoClone':
                    continue
                tag = '+'.join([batch.tag, use_cellcolors, batch.get_onlytag()])
                out_fn = str(self.bipartite_file_tmpl).format(tag)
                self.bipartite_files.append(Path(out_fn))
                BipartiteTCR.draw(batch, out_fn=str(self.out_dir/out_fn), 
                                  use_cellcolors=use_cellcolors, 
                                  interactive=interactive)
                status['bipartite_{}'.format(use_cellcolors)] = OK
            
        if indict('tcrs', conf):
            display(batch.tcrs.tcrs)

        def generate_config_file(cell_part, fn):
            print('Dumping DE config to ', fn)
            groups = list(cell_part.part_dict.keys())
            if len(groups) > 2:
                if CellPart.CLONAL in groups and CellPart.NONCLONAL in groups:
                    groups = [ CellPart.CLONAL, CellPart.NONCLONAL ]
            #assert len(groups) == 2
            if len(groups) == 2:
                if CellPart.EARLY_CLONAL in groups:
                    assert CellPart.LATE_CLONAL in groups
                    groups = [ CellPart.LATE_CLONAL, CellPart.EARLY_CLONAL ]
                if CellPart.EARLY_NONCLONAL in groups:
                    assert CellPart.LATE_NONCLONAL in groups
                    groups = [ CellPart.LATE_NONCLONAL, CellPart.EARLY_NONCLONAL ]
            with open(fn, 'w') as f:
                json.dump(groups, f, indent=4)
            return groups
            
        if indict('de', conf):
            #assert indict('clonality', conf)  # `dumper' needed
            sample_name = batch.tag
            kallisto_dir = conf['kallisto_dir_tmpl'].format(sample_name)
            #clonality_file = '/home/pesho/repos/bioviz/in/clonality/clonality_{}.json'.format(sample_name)
            clonality_file = self.out_dir / str(self.clonality_file_tmpl).format(batch.tag)
            gene_lists = conf['de']['gene_lists'] if 'gene_lists' in conf['de'] else [ 'all' ]

            for gene_list in conf['de']['gene_lists']:
                for de_method in conf['de']['methods']:
                    for clear_groups in conf['de']['clear_groups']:
                        clear_groups_str = 'clear_groups' if clear_groups else 'extended_groups'
                        experiment_str = 'de_{}_{}'.format(de_method, clear_groups_str)
                        gene_list_name = 'all' if gene_list == 'all' else str(Path(gene_list).stem)

                        out_fn = self.out_dir / str(self.de_file_tmpl4).format(gene_list_name, batch.tag, de_method, clear_groups_str)
                        out_dir = out_fn.parent
                        out_dir.mkdir(parents=True, exist_ok=True)

                        de_config_file = self.out_dir / str(self.de_config_file_tmpl4).format(gene_list_name, batch.tag, de_method, clear_groups_str)
                        groups = generate_config_file(batch.cell_part, de_config_file)
                        group_sizes = [ len(batch.cell_part.part_dict[group]) for group in groups ]

                        #X_size, Y_size = batch.cell_part.get_group_sizes(clear_groups=clear_groups)
                        #if min(X_size, Y_size) < conf['de']['min_group_size']:
                        assert indict('min_group_size', conf['de'])
                        if min(group_sizes) < conf['de']['min_group_size']:
                            status[experiment_str] = SKIPPED
                            print('Skipping DE on {} because of small group sizes: {}'.format(sample_name, ', '.join([str(x) for x in group_sizes])))
                            continue

                        #dgelist_rds_fn = '~/work/celldive/tmp/dgelist.rds'
 
                        args = [ 'Rscript', 'DE/de.R', sample_name,  str(clonality_file), de_method, str(out_dir), str(out_fn), str(int(clear_groups)), gene_list, str(de_config_file) ]
                        #print('running.py DE:', ' '.join(args))
                        if subprocess.call(args, stderr=subprocess.STDOUT):
                            print(red("ERROR: The DE/de.R didn't run successfully on sample {} with args {}".format(
                                sample_name, ' '.join(args))))
                            #unsuccessful['de'].append(batch.tag)
                            status[experiment_str] = ERROR
                        else:
                            status[experiment_str] = OK

                        # Visualizing statistics on DE results (using Vidger)
                        #args = [ 'Rscript', 'DE/vis.R', rdslist_rds_fn, sample_name, str(out_dir) ]
                        #experiment_str = 'de_vis_{}_{}'.format(de_method, clear_groups_str)
                        #if subprocess.call(args, stderr=subprocess.STDOUT):
                        #    print(red("ERROR: The DE/vis.R didn't run successfully on sample {} with args {}".format(
                        #        sample_name, ' '.join(args))))
                        #    #unsuccessful['de'].append(batch.tag)
                        #    status[experiment_str] = ERROR
                        #else:
                        #    status[experiment_str] = OK

  #      if indict('mutations', conf):
 #           if indict('cnv', conf['mutations']):
 #               out_dir = self.out_dir / str(self.mutations_cnv_dir_tmpl3).format(batch.tag)
 #               out_dir.mkdir(parents=True, exist_ok=True)
 #               args = [ 'Rscript', 'DE/kallisto2genetable.R', sample_name, kallisto_dir, str(clonality_file), de_method, str(out_dir), str(out_fn), str(int(clear_groups)) ]
#    sample_name, clonality_file, countsFromAbundance, gene.table.fn, clear.groups

 #               if subprocess.call(args, stderr=subprocess.STDOUT):
 #                   print(red("ERROR: The DE/de.R didn't run successfully on sample {} with args {}".format(
 #                       sample_name, ' '.join(args))))
 #                   #unsuccessful['de'].append(batch.tag)
 #                   status[experiment_str] = ERROR
 #               else:
 #                   status[experiment_str] = OK

        if indict('vb_kit_distr', conf):
            out_fn = str(self.vb_kit_distr_file_tmpl).format(batch.tag)
            VbTrbv.draw_all(batch, out_fn=self.out_dir/out_fn, interactive=interactive)
            status['vb_kit_distr'] = OK

        #if indict('vj_heatmap', conf):
        #    out_fn = str(self.vj_heatmap_file_tmpl).format(batch.tag)
        #    TCRheatmap.draw(batch, interactive=interactive, out_fn=self.out_dir/out_fn)
        
        #if indict('gsea', conf):
        #    out_dir = self.out_dir / str(self.gsea_folder_tmpl).format(batch.tag)
        #    GSEA(batch, geneset=conf['gsea']['geneset'], out_dir=out_dir)

        #if indict('stats', conf):
        #    batch.print_stats()
                
        #if indict('forcegraph', conf):
        #    draw_TCR_graph(batch)
        #if indict('quantheatmap', conf):
        #    self.draw_quantheatmap(batch, features=features)
        #if indict('cross_validation', conf):
        #    self.gcv.calc_batch_grid_cv(batch, feature=features)
        #    
        #if indict('scatter', conf):
        #    for dim_reduct in conf['scatter']['dim_reduct']:
        #        self.draw_scatter(batch, dim_reduct=dim_reduct, perplexity=perplexity,
        #                          interactive=interactive, features=features)

        return status

    def do_combine_batches(self, combine_batches):
        ''' combine_batches in the form { ('02_clonal_progression', ('02_SE', 'clonal'), ('02_SL', 'clonal')) }
        '''
        print()
        print('Combining batches...')
        # per patient
        ans_dict = {}

        for pair in combine_batches:
            sample_id, _a, _b = pair
            a_sample, a_part = _a
            b_sample, b_part = _b

            a = self.B[a_sample] if a_sample in self.B else ans_dict[a_sample]
            b = self.B[b_sample] if b_sample in self.B else ans_dict[b_sample]

            internal_name = a.internal_name + '_and_' + b.internal_name
            #colored_tags = [xkcd2console(b.tag if b.tag else "", b.color) for b in batches]
            #print('Combining {}\'s batches: {}'.format(blue(patient), ', '.join(colored_tags)))

            #if partition == 'clonalVsNonclonal':
            #    assert False, 'DEBUG: cell_part | cell_part'
            #    sumcell_part = a.cell_part | b.cell_part
            #elif partition == 'earlyVsLateClonal':
            #    sumcell_part = a.cell_part.earlyVsLateclonalVsClonal(b.cell_part)
            #else:
            #    assert False, "No such partition: " + partition

            if a_part == 'clonal' and b_part == 'clonal':
                sumcell_part = a.cell_part.earlyVsLateClonal(b.cell_part)
            elif a_part == 'single' and b_part == 'single':
                sumcell_part = a.cell_part.earlyVsLateSingle(b.cell_part)
            elif a_part == 'all' and b_part == 'all':
                sumcell_part = a.cell_part.merge(b.cell_part)
            else:
                assert False, "No such partition: " + a_part + " and " + b_part
            
            ans_dict[sample_id] = Batch.merge(a, b, sample_id=sample_id, internal_name=internal_name, cell_part=sumcell_part)

        return ans_dict

    #def combine_batches(self, combine_batches):
    #    print()
    #    print('Combining batches...')
    #    # per patient
    #    ans_dict = {}
    #    patients = {}   # patient -> merged batches
    #    if 'patient' in combine_batches or 'all' in combine_batches:
    #        #print(self.samples)
    #        #for patient in set(list(zip(*self.samples))[:2]):
    #        for patient in { sample[:2] for sample in self.samples }:
    #            batches = []
    #            for key, batch in self.B.items():
    #                print('key: ', key, 'patient: ', patient)
    #                p = key[:2]
    #                if p == patient:
    #                    batches.append(batch)
    #            if len(batches) > 1:
    #                colored_tags = [xkcd2console(b.tag if b.tag else "", b.color) for b in batches]
    #                print('Combining patient {}\'s batches: {}'.format(blue(patient), ', '.join(colored_tags)))
    #                b = Batch.merge(*batches, tag=patient)
    #                patients[patient] = b
    #                if 'patient' in combine_batches:
    #                    ans_dict[patient] = b
    #                    
    #    if 'all' in combine_batches:
    #        # all patients together
    #        if len(patients) > 1:
    #            print('Combining {} patients: {}'.format(blue('ALL'), ', '.join([b.tag for b in patients])))
    #            b = Batch.merge(*patients.values(), tag='ALL')
    #        else:
    #            b = next(iter(patients.values()))
    #        ans_dict['ALL'] = b

    #    if 'pairs' in combine_batches:
    #        # pairs of batches
    #        for kv1, kv2 in itertools.combinations(self.B.items(), 2):
    #            k1, b1 = kv1
    #            k2, b2 = kv2
    #            p1, s1 = k1
    #            p2, s2 = k2
    #            if p1 == p2:
    #                if s1[0] == s2[0]: # or s1[-1] == s2[-1]:
    #                    tag = 'pair:{}+{}'.format(b1.tag, b2.tag)
    #                    for p in ['clonalVsNonclonal']: #, 'clonalVsClonal']:
    #                        t = tag + '_' + p
    #                        print('Combining pairs of batches: {} as {}'.format(blue(p1), ', '.join([b1.tag, b2.tag]), p))
    #                        b12 = Batch.merge(b1, b2, tag=t, partition=p)
    #                        ans_dict[t] = b12

    #    #if 'progression' in combine_batches:

    #        
    #    print('{} new batches produced: {}'.format(
    #        green(len(ans_dict)), blue(', '.join(ans_dict.keys()))))
    #    return ans_dict

    def compile_report(self):
        def notebook_copy(in_file, out_file):
            print('Make a copy of the notebook...')
            args = ['jupyter', 'nbconvert', '--to', 'html', in_file, '--output', out_file]
            if subprocess.call(args, stderr=subprocess.STDOUT):
                print(red("ERROR: The notebook didn\'t copy!"))
                self.unsuccessful['jupyter-nbconvert'].append("notebook_copy")
        
        def add_filenames_to_images(in_files):
            print('Collect {} images and stamp them...'.format(blue(len(in_files))))
            args = ['mogrify', '-font', 'Liberation-Sans', '-fill', 'white',
                    '-undercolor', '"#00000080"', '-pointsize', '30',
                    '-gravity', 'SouthEast', '-annotate', '+10+10', '%t', *in_files]
            #mogrify -font Liberation-Sans -fill white -undercolor '#00000080' -pointsize 30 -gravity SouthEast -annotate +10+10 %t *.png  # TODO png -> EXT
            if subprocess.call(args, stderr=subprocess.STDOUT):
                Warning('Adding filename annotations to images failed!')
                self.unsuccessful['mogrify'].append("add_filenames_to_images")

        def images2pdf(in_files, report_fn):
            print('Merge all images to one PDF [{}]...'.format(blue(str(report_fn.name))))
            args = ['convert', *in_files, "-quality", "300", str(report_fn)]
            if subprocess.call(args, stderr=subprocess.STDOUT):
                print(red("ERROR: PDF report not compiles!"))
                self.unsuccessful['convert'].append("imeges2pdf")
        
        def dir2pdf(d, out_fn):
            in_files = sorted([ str(f) for f in Path(d).glob('**/*' + '.' + EXT) ])
            if not in_files:
                Warning('PDF report not compiled! No image files in {}'.format(blue(str(self.out_dir))))
            else:
                #add_filenames_to_images(in_files)
                images2pdf(in_files, out_fn)
                
                    
        print()
        print(blue("   ------- [PDF] --------"))
        
        print('Dumping config to ', self.config_dumped_file)
        with open(self.out_dir/self.config_dumped_file, 'w') as f:
            json.dump(self.conf, f, indent=4)

        print('Copying sample descriptions to {}'.format(blue(self.sample_description_dump_file)))
        args = ['cp', '-v', self.samples_description_tsv, self.out_dir/self.sample_description_dump_file]
        if subprocess.call(args, stderr=subprocess.STDOUT):
            print(red("ERROR: Copying samples description did not finish."))
            self.unsuccessful['dump'].append("copying sample descriptions")

        #notebook_copy(dirs['src']/self.notebook, self.out_dir/self.notebook)
        
        for d in Path(self.out_dir).iterdir():
            if d.is_dir():
                report_fn = self.out_dir / '{}_{}.pdf'.format(d.name, str(self.folder_tag))
                dir2pdf(d, report_fn)
            
        #report_fn = self.out_dir / '{}_{}.pdf'.format(self.experiment, str(self.folder_tag))
        #dir2pdf(self.out_dir, report_fn)
                
        print(blue("   ------- [DONE] --------"))
        
    #def produce_results_for_batch(self, batch):
    #    print()
    #    print(batch)
    #    if type(self.features) is str:
    #        self.features = [self.features]

    #    for f in self.features:
    #        self.exec_per_batch(batch, graphs_dict=self.conf)

    def run_in_parallel(self, f, items):
        '''Runs `f' in parallel on every item of `items'.
           `f' has to accept two arguments: an item and the dict `unsuccessful' for reporting errors.'''
        processes = self.conf['processes'] if 'processes' in self.conf else None

        sample2res = {}

        if processes == 1:
            for item in items:
                sample2res[item] = f(item, self.unsuccessful, self.B)
        else:
            with Pool(processes) as pool:
                try:
                    handlers = []
                    for item in items:
                        handler = pool.apply_async(f, args = (item, self.unsuccessful, self.B),)
                        handlers.append((item, handler))
                    print('Parallel execution of {} tasks in <= {} processes'.format(len(handlers), pool._processes), flush=True)
                    
                    for item, handler in handlers: 
                        sample2res[item] = handler.get()
                except Exception as e:
                    Warning(e)
                    raise Exception("Parallel execution stopped.") from e #"".join(traceback.format_exception(*sys.exc_info())))
                pool.close() 
                pool.join()

        return sample2res

    def print_final_msg(self, sample2statuses):
        print()
        proc2failedsamples = defaultdict(list)
        error_cnt = 0
        for sample, statuses in sample2statuses.items():
            for proc, status in statuses.items():
                proc2failedsamples[proc]
                if status != 'OK':
                    if status == 'ERROR':
                        status = red(status)
                    proc2failedsamples[proc].append((sample, status))
                    error_cnt += 1

        #print('Samples: {}'.format(', '.join([ sample.__repr__() for sample in self.conf['samples'] ])))
        print('Samples: {}'.format(', '.join([ sample.__repr__() for sample in self.B.keys() ])))

        color_f = red if error_cnt else green
        msg = '{} procesures failed.'.format(error_cnt) if error_cnt else 'All runs are successful.'
        print(color_f(' --- Status: {} --- '.format(msg)))
        for proc, l in proc2failedsamples.items():
            if l:
                line = ', '.join([ '{}({})'.format(sample.__repr__(), status) for sample, status in l ])
                print('   {}: {}'.format(red(proc), line))
            else:
                print('   {}: {}'.format(proc, green('OK on all samples')))
        print(color_f(' --------- '))
        
    def produce_results(self):
        conf = self.conf
        #assert((not 'geneintersect' in conf) or (not 'de' in conf)), \
        #        "If 'de' is specified, the genes are automatically intersected. No 'geneintersect' needed!"

        #if indict('geneintersect', conf):
        #    d = conf['geneintersect'] 
        #    if not Path(d).is_dir():
        #        raise Exception('No dir {}'.format(d))
        #elif indict('de', conf):
        #    conf['geneintersect'] = str((self.out_dir/self.de_file_tmpl4.format(de_method)).parent)

        sample2statuses = self.run_in_parallel(self.exec_per_batch, self.B.values())
        #sample2statuses_str = '\n'.join([ '{}: {}'.format(sample, status) for sample, status in sample2statuses.items() ])
        #print('Sample statuses: {}'.format(sample2statuses_str))
 
        #if 'cross_validation' in conf and conf['cross_validation']:
        #    #out_fn = str(self.out_dir / self.cv_file_tmpl).format('-'.join(self.B.keys()))
        #    out_fn = str(self.out_dir / self.cv_file_tmpl).format('cross_validation')
        #    self.gcv.results2html(out_fn)

        if indict('clonality', conf):
            files_dir = (self.out_dir / self.clonality_file_tmpl.format("example")).parent
            out_file = self.out_dir / self.clonality_stats_file
            if files_dir.exists():
                ClonalityStats(files_dir).dump(out_file)
            else:
                warning('No clonality dir {}'.format(files_dir))

        #if indict('meta', conf):
            #display()
            #files_dir = (self.out_dir / self.meta_json_tmpl.format("example")).parent
            #out_file = self.out_dir / self.meta_stats_file
            #ClonalityStats(files_dir).dump(out_file)
            
            #stats = pd.DataFrame()
            ##print(self.B.items())
            #for patient_sample, batch in self.B.items():
            #    out_counts_fn = str(self.out_dir / self.clonality_stats_counts_file_tmpl).format(batch.tag)
            #    out_cells_fn = str(self.out_dir / self.clonality_stats_cells_file_tmpl).format(batch.tag)
            #    row = ClonalityStats(batch).run(out_counts_fn, out_cells_fn)
            #    stats = stats.append(row)
            #display(stats) 
            #print('Writing to', self.clonality_stats_file)
            #stats.to_csv(self.out_dir / self.clonality_stats_file)

        if indict('compare_gene_list', conf) and indict('de', conf):
            for gene_list in conf['de']['gene_lists']:
                gene_list_name = 'all' if gene_list == 'all' else Path(gene_list).stem #name
                de_dir = self.out_dir / str(Path(self.de_file_tmpl4.format(gene_list_name, 'bala', 'nica', 'bla')).parent)
                de_dir.mkdir(parents=True, exist_ok=True)
                gene_list_copy_fn = de_dir / Path(conf['compare_gene_list']).name
                args = [ 'cp', conf['compare_gene_list'], str(gene_list_copy_fn) ]
                print(args)
                if subprocess.call(args, stderr=subprocess.STDOUT):
                    print(red("ERROR: Cannot copy the file."))

        if indict('de', conf):
            for gene_list in conf['de']['gene_lists']:
                gene_list_name = 'all' if gene_list == 'all' else Path(gene_list).stem #name
                if indict('gene_table_intersect', conf['de']):
                    de_dir = conf['de']['gene_table_intersect']
                else:
                    de_dir = self.out_dir / str(Path(self.de_file_tmpl4.format(gene_list_name, 'bala', 'nica', 'bla')).parent)

                gene_intersect = GeneIntersect(de_dir)

                methods = conf['de']['methods'] if indict('methods', conf['de']) else []
                if indict('clear_groups', conf['de']):
                    clear_groups = [ 'clear_groups' if clear_groups else 'extended_groups' for clear_groups in conf['de']['clear_groups'] ]
                else:
                    clear_groups = []

                #groups = [''] + self.samples + methods + clear_groups  # '' creates a group of all samples
                if 'combine_batches' in conf:
                    methods += ['clonal_progression']
                    methods += ['single_progression']
                    methods += ['progression']
                if len(methods) > 1:
                    methods += ['']

                print('Separate tables for {} to be generated...'.format(', '.join(methods)))

                for fn_should_include in methods:
                    fn2sample = gene_intersect.get_fn2sample(fn_should_include)

                    out_fn_tmpl = self.out_dir / self.gene_table_file_tmpl2.format(gene_list_name, fn_should_include, "{}")
                    gene_intersect.dump_gene_tables(out_fn_tmpl=out_fn_tmpl, fn2sample=fn2sample)

                    stats_fn = self.out_dir / self.gene_stats_file_tmpl.format(gene_list_name, fn_should_include)
                    gene_intersect.dump_stats(out_fn=stats_fn, fn2sample=fn2sample)

                # square tables with samples compared to themselves
                for group in [ [''], methods, ]: #clear_groups ]:
                    if group:
                        #print('All combinations among {}:'.format(', '.join(group)))
                        #for fn_should_include_X, fn_should_include_Y in itertools.combinations_with_replacement(group, 2):
                        print('Sqare table for {}:'.format(', '.join(group)))
                        for fn_should_include_X in group:
                            fn_should_include_Y = fn_should_include_X
                            fn2sample_X = gene_intersect.get_fn2sample(fn_should_include_X)
                            fn2sample_Y = gene_intersect.get_fn2sample(fn_should_include_Y)

                            out_fn_tmpl = self.out_dir / self.de_sample_by_samples_table.format(gene_list_name, fn_should_include_X, fn_should_include_Y, '{}')
                            gene_intersect.dump_experiment_comparison_table(out_fn_tmpl=out_fn_tmpl,
                                    fn2sample_X=fn2sample_X, fn2sample_Y=fn2sample_Y)


            # UpSeT
            #de_method = checkbox['de']['method']
            #de_dir = self.out_dir / str(Path(self.de_file_tmpl4.format('ala', 'bala', 'nica', 'bla')).parent)
            #out_fn = self.out_dir / str(self.upset_file)
            #out_fn.parent.mkdir(parents=True, exist_ok=True)
            #args = ['Rscript', 'DE/upset.R', str(de_dir), str(out_fn)]
            #print(args)
            #if subprocess.call(args, stderr=subprocess.STDOUT):
            #    print(red("ERROR: The DE/upset.R didn't run successfully!"))
            #    self.unsuccessful['upset'].append('geneintersect')
            
        self.compile_report() # TODO: move out_fn as default in GeneIntersect

        self.print_final_msg(sample2statuses)
    
    def read_samples_description_tsv(self, fn):
        self.samples_description = pd.read_csv(fn, sep='\t', index_col='sample_id')
        print('Descriptions of {} samples: {}'.format(self.samples_description.shape[0],
                                                      ', '.join(self.samples_description.index)))
    def read_sample(self, sample_id, unsuccessful, samples):
        conf = self.conf
        #patient = color, clonotype in samples:
        row = self.samples_description.loc[sample_id, :]
        color, relative_dir, internal_name, facs = row['color'], Path(row['directory']), row['internal_sample_name'], row['facs_sorting']
        alpha, beta = row['clonal_alpha'], row['clonal_beta']
        #clonotype = self.clonality_rule.replace('alpha', alpha).replace('beta', beta)
        #clonality_method = self.clonality['method'] if self.clonality else None
        clonotype = OrderedDict([('alpha', alpha), ('beta', beta), ('rule', self.clonality['rule']), ('facs', facs)])
        clonality = conf['clonality'] if 'clonality' in conf else None
        self.B.pop(sample_id, None)
        assert not sample_id in self.B
        normalization = conf['normalize'] if 'normalize' in conf else None
        filt_cells = conf['filt_cells'] if 'filt_cells' in conf else None
        filt_tcr_chains = conf['filt_tcr_chains'] if 'filt_tcr_chains' in conf else None
        filt_genes = conf['filt_genes'] if 'filt_genes' in conf else None
        #conf_de = conf['de'] if 'de' in conf else None
        expression_threshold = filt_genes['expression_threshold'] if filt_genes and 'expression_threshold' in filt_genes else None
        min_cells_expressing_gene = filt_genes['min_cells_expressing_gene'] if filt_genes and 'min_cells_expressing_gene' in filt_genes else None
        full_dir = Path(self.samples_description_tsv).parent / relative_dir
        
        kallisto_dir = conf['kallisto_dir_tmpl'].format(sample_id)
        batch = Batch.from_file(full_dir, sample_id=sample_id, internal_name=internal_name, color=color,
            kallisto_dir=kallisto_dir,
            diffexpr=self.diffexprs[0], input_files=self.input_files,  # TODO: fix [0]
            filt_cells=filt_cells,
            filt_tcr_chains=filt_tcr_chains,
            clonality=clonality,
            clonotype=clonotype,
            #conf_de=conf_de,
            expression_threshold=expression_threshold,
            min_cells_expressing_gene=min_cells_expressing_gene,
            nrows=self.nrows, out_dir=self.out_dir, norm=normalization)
        
        #_b.get_clones_lazy(onlytag=clonotype, diffexpr=self.diffexprs[0])

        #_b = _b.subbatch(features=features)
        # def subbatch(self, cells, *, tag=None, features=None, color=False):
        #_b = _b.sieved(filt)
        #_g = FilterBatch.filter_cells(filt_cells)  ## TODO: uncomment!
        #self.B[sample_id] = batch
        #samples[sample_id] = batch
        #print('inside: ',samples)
        return batch
       
    def read_samples(self, samples):
        print('Loading samples...', flush=True)
        id2sample = self.run_in_parallel(self.read_sample, samples)
        for sample_id, sample in id2sample.items():
            self.B[sample_id] = sample

    def run_experiment(self, *, log_fn):
        sys.stdout = Tee(log_fn, 'w')

        start_time = time.time()
        self.read_samples(self.samples)

        if 'combine_batches' in self.conf:
            combined_batches = self.do_combine_batches(self.conf['combine_batches'])
            self.B.update(combined_batches)
            
        self.produce_results()

        (self.out_dir / 'DONE').touch()   # Magic file OK upon sucess
        execution_time = time.time() - start_time
        print('Execution time: {} sec'.format(int(execution_time)))
        sys.stdout.flush()
        
    def init(self, *, clear=True, conf={}):
        if clear:
            self.__init__(self.dirs)
        
        self.conf = conf
        self.input_files = conf['input_files'] if 'input_files' in conf else {}
        self.diffexprs = conf['diffexprs'] if 'diffexprs' in conf else ['bla']  # TODO: remove
        self.nrows = conf['nrows'] if 'nrows' in conf else 0
        self.interactive = conf['interactive'] if 'interactive' in conf else False
        #self.features = conf['features']  # 'important' or 'all'
        self.experiment = conf['experiment'] if 'experiment' in conf else ''
        if 'notebook' in conf:
            self.notebook = conf['notebook']
        else:
            self.notebook = self.experiment + '.ipynb'
        self.samples = conf['samples'] if indict('samples', conf) else []
        self.samples_description_tsv = conf['samples_description_tsv']
        self.clonality = conf['clonality']
        if 'kernels' in conf:
            self.kernels = conf['kernels']  # SVM kernels
            self.Cs = conf['Cs']  # 'SVM coefs'
            self.gammas = conf['gammas']
        else:
            self.kernels = None
            self.Cs = [1]
            self.gammas = [1]
        if self.nrows == 0:
            self.nrows = None
            
        #self.gcv = GridCrossValidation(clf=HypoClassifier, features=self.features, diffexprs=self.diffexprs,
        #                               Cs=self.Cs, kernels=self.kernels, gammas=self.gammas)
        if 'filt_cells' in conf and 'cells' in conf['filt_cells']:
            filt = conf['filt_cells']['cells']
        else:
            filt = None
                
        self.read_samples_description_tsv(self.samples_description_tsv)
        if type(self.samples) is str and self.samples == 'all':
            self.samples = list(self.samples_description.index)

        self.set_output_dir(nrows=self.nrows, filt=filt)
        execution_log_fn_colorful = self.out_dir / self.execution_log_file_tmpl1.format('colorful')

        self.run_experiment(log_fn=execution_log_fn_colorful)

        # colorful execution log to BW log
        execution_log_fn_bw = self.out_dir / self.execution_log_file_tmpl1.format('bw')
        #sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g"
        bw_f = open(execution_log_fn_bw, 'w')
        args = [ 'sed', '-r', 's/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g', str(execution_log_fn_colorful) ]
        print(args)
        if subprocess.call(args, stderr=subprocess.STDOUT, stdout=bw_f):
            print(red("ERROR: Sed failed to convert the colorful {} to BW {}".format(execution_log_fn_colorful, execution_log_fn_bw)))

