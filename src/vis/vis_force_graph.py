from string import Template
from IPython.core.display import HTML, display
import random
from pathlib import Path
from utils import *

_vis_dir = Path('.') / 'vis'
_force_graph_js = str(_vis_dir / 'force_directed_graph.js')

def init_d3():
    html = []
    html.append('<script src="%s/lib/d3.min.js"></script>' % _vis_dir)
    css = '''
    .node {
      stroke: #fff;
      stroke-width: 1.5px;
    }

    .link {
      stroke: #999;
      stroke-opacity: .6;
    }
    '''
    html.append("<style>" + css + "</style>")
    return '\n'.join(html)

def get_html_graph(graph):
    # http://bl.ocks.org/mbostock/4062045
    JS_text = Template('''
                <div id='maindiv${divnum}'></div>
                <script>
                    $main_text
                </script>
                ''')

    divnum = int(random.uniform(0,9999999999))
    data_dict = { 'divnum': divnum, 'data': graph }
    assert_file_exists(_force_graph_js)
    main_text_template = Template(open(_force_graph_js, 'r').read())
    main_text = main_text_template.safe_substitute(data_dict)

    return JS_text.safe_substitute({'divnum': divnum, 'main_text': main_text})

#def test_construct_graph():
#    display(HTML(viz_graph_d3.init_d3()))
#    
#    n_nodes = 30
#    p_edge = 0.05
#    graph = {"nodes": [], "links": []}
#    for i in range(n_nodes):
#        graph["nodes"].append( {"name": "i" + str(i),
#                                #"group": int(random.uniform(1,11)),
#                               } )
#    for i in range(n_nodes):
#        for j in range(n_nodes):
#            if random.uniform(0,1) < p_edge:
#                graph["links"].append( {"source": i,
#                                        "target": j,
#                                        "value": random.uniform(0.5,3)} )
#    display(HTML(viz_graph_d3.draw_graph(graph)))

def construct_tcr_graph(batch):
    graph = {"nodes": [], "links": []}
    
    cells = batch.tcrs.cells()
    for cell in cells:
        sample = batch.meta.cell2sample(cell)
        if not sample:
            sample = 999
        graph["nodes"].append( {"name": str(cell),
                                'color': batch.meta.get_color_rgb(cell),
                                #"group": sample,
                               } )
        
    for a, b, tcr in batch.tcrs.getEdges():
        if a in cells and b in cells:
            graph["links"].append( {"source": cells.index(a),
                                    "target": cells.index(b),
                                    "value": random.uniform(0.5,3)} )
    
    return graph

_initialized = False

def draw_TCR_graph(batch):
    global _initialized
    if not _initialized:
        display(HTML(init_d3()))
        _initialized = True
        print('initialized!')
    graph = construct_tcr_graph(batch)
    display(HTML(get_html_graph(graph)))