import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from utils import *

class Violin:
    @classmethod
    def draw(cls, genes_df, out_fn=None, interactive=False):
        #if genes_df is None or genes_df.empty:
        #    warning('No genes to be plotted as violins!')
        #    return
        
        
        #x = np.random.poisson(lam =3, size=100)
        #y = np.random.choice(["S{}".format(i+1) for i in range(6)], size=len(x))
        #df = pd.DataFrame({"Scenario":y, "LMP":x})
        
        numpy_arrays = genes_df.values
        
        if numpy_arrays is None:
            warning('No genes for a violin to be drawns.')
            return

        plt.ioff()
        fig, ax = plt.subplots(1, 1)

        ax.violinplot(dataset = numpy_arrays)

        #ax.set_title(batch.label)
        ax.yaxis.grid(True)
        ax.set_xlabel('Genes')
        ax.set_ylabel('Expression')
        
        ax.set_xticks(range(1, len(genes_df.columns)+1))
        ax.set_xticklabels(genes_df.columns) #, rotation=-30) #, ha='center')
        
        plot(plt, interactive, out_fn)
        
    @classmethod
    def draw_example(cls, out_fn=None, interactive=False):
        plt.ioff()
        
        x = np.random.poisson(lam =3, size=100)
        y = np.random.choice(["S{}".format(i+1) for i in range(6)], size=len(x))
        df = pd.DataFrame({"Scenario":y, "LMP":x})

        fig, ax = plt.subplots(1, 1)

        ax.violinplot(dataset = [df[df.Scenario == 'S1']["LMP"],
                                   df[df.Scenario == 'S2']["LMP"],
                                   df[df.Scenario == 'S3']["LMP"],
                                   df[df.Scenario == 'S4']["LMP"],
                                   df[df.Scenario == 'S5']["LMP"],
                                   df[df.Scenario == 'S6']["LMP"] ] )

        ax.set_title('Day Ahead Market')
        ax.yaxis.grid(True)
        ax.set_xlabel('Scenario')
        ax.set_ylabel('LMP ($/MWh)')
        
        plot(plt, interactive, out_fn)