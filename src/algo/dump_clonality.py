import json
from collections import OrderedDict #defaultdict
from pathlib import *
from utils import *

from structs.clone_struct import CellPart

class DumpClonality: 
    def __init__(self, batch, kallisto_dir):
        self.batch = batch
        self.kallisto_dir = kallisto_dir

    def cells2json(self, cells):
        d = OrderedDict() #defaultdict(list)
        #cells = list(cells)
        tcr_cells = self.batch.tcrs.cells() 
        #if len(cells) > 0 and cells[0] not in self.batch.tcrs.cells():
        #    print('requested cells: ', cells[:3])
        #    print('available cells', self.batch.tcrs.cells()[:3])
        for cell in sorted(cells):
            if cell in tcr_cells:
                d[cell] = self.batch.tcrs.cell2struct(cell)
                #d[cell] = list()
                #for seqs in self.batch.tcrs.cell2seqs(cell):
                #    d[cell].append(seqs)
            else:
                d[cell] = {}
            d[cell]['microscope'] = int(self.batch.meta.well2count(cell))
            #f = self.batch.cell_part.cell2file(cell)
            #d[cell]['kallisto_subdir'] = str(Path(self.kallisto_dir) / f) if f else None
            d[cell]['kallisto_subdir'] = self.batch.quant.cell2kallisto(cell)
        return d

    def dump(self, fn):
        part = self.batch.cell_part
        if not part:
            raise Exception('No cell partitioning for {}'.format(self.batch.tag))
        
        cell_part = OrderedDict()
        cell_part['clonotype'] = part.clonotype
        # TODO: compute independent of CellPart
        members = [ attr for attr in dir(CellPart) if not callable(getattr(CellPart, attr)) and not attr.startswith("__") ]
        parts = [ part.part_dict.get(getattr(CellPart, member), {}) for member in members ]
        #print('parts: ', parts)
        all_cells = set().union(*parts)
        #print('all cells: ', all_cells)

        #all_cells = part.part_dict[CellPart.CLONAL] | part.part_dict[CellPart.RELATED] | part.part_dict[CellPart.AMBIGUOUS] | part.part_dict[CellPart.NONCLONAL] | part.part_dict[CellPart.BYSTANDERS]
        #print('dump_clonality: sample {} '.format(self.batch.tag), cell_part)
        clonotype = cell_part['clonotype']
        #clonotype['theta_cells'] = part.theta
        #clonotype['theta_perc'] = float('{:.2f}'.format(100.0 * part.theta / len(all_cells)))
        clonotype['max_cc_size [cells]'] = part.max_tc_size
        clonotype['max_cc_size [%]'] = float('{:.2f}'.format(100.0 * part.max_tc_size / len(all_cells)))
        clonotype['reconstructed'] = len(all_cells)

        for group, cells in part.part_dict.items():
            cell_part[group] = self.cells2json(cells)

        to_json = {
            self.batch.tag: cell_part,
        }

        if Path(fn).suffix != 'json':
            Warning('{} is in json format without a json extension.'.format(fn))

        Path(fn).parent.mkdir(parents=True, exist_ok=True)
        with open(fn, 'w') as f:
            json.dump(to_json, f, indent=4)

        print('Clonality JSON dumped: {}'.format(blue(Path(fn).name)))

