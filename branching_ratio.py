
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',   
    r'\usepackage{sansmath}',  
    r'\sansmath'               
    r'\usepackage{siunitx}',   
    r'\sisetup{detect-all}',   
]  

import argparse, sys, os, itertools, pickle
import numpy as np
from brian2.units import mV, ms, second


from .axes.network import *
from .axes.neuron import *
from .axes.synapse import *
from .axes.parameter_display import *


def branching_ratio_figure(bpath, nsp):

    fig = pl.figure()
    ax_lines, ax_cols = 1,1
    axs = {}
    for x,y in itertools.product(range(ax_lines),range(ax_cols)):
        axs['%d,%d'%(x+1,y+1)] = pl.subplot2grid((ax_lines, ax_cols), (x, y))

    fig.set_size_inches(6*1.6/2,2.*1.6)


    # for _,ax in axs.items():
    #     ax.axis('off')

    tmin1, tmax1 = 0*second, nsp['T1']
    tmin3, tmax3 = nsp['T1']+nsp['T2'], nsp['T1']+nsp['T2']+nsp['T3']

    tmin5 =  nsp['T1']+nsp['T2']+nsp['T3']+nsp['T4']
    
    # raster_plot(axs['1,1'], bpath, nsp, tmin=tmin1, tmax=tmax1)
    raster_plot(axs['1,1'], bpath, nsp, tmin=tmin5, tmax=tmin5+5*second)

    bin_w = 10*ms

    ft = branching_ratio('', bpath, nsp, bin_w)

    axs['1,1'].text(0.6, -0.5, 'mre %.6f' %ft.mre,
                    transform=axs['1,1'].transAxes)

    axs['1,1'].text(0.05,-0.5, 'ATotalMax %.2f' %nsp['ATotalMax'],
                    transform=axs['1,1'].transAxes)
    axs['1,1'].text(0.05,-0.65, 'iATotalMax %.4f' %nsp['iATotalMax'],
                    transform=axs['1,1'].transAxes)


    pl.tight_layout()

    directory = "figures/branching_ratio"
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    pl.savefig(directory+"/{:s}.png".format(bpath[-4:]),
               dpi=300, bbox_inches='tight')

    



    
if __name__ == "__main__":

    
    # return a list of each build (simulations run)
    # e.g. build_dirs = ['builds/0003', 'builds/0007', ...]
    # sorted to ensure expected order
    build_dirs = sorted(['builds/'+pth for pth in next(os.walk("builds/"))[1]])
    
    for bpath in build_dirs:

        try:
            with open(bpath+'/raw/namespace.p', 'rb') as pfile:
                nsp=pickle.load(pfile)
                
            # if bpath[-4:]=='0000':

            branching_ratio_figure(bpath, nsp)

        except FileNotFoundError:
            print(bpath[-4:], "reports: No namespace data. Skipping.")
