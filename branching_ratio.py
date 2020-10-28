
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

import mrestimator as mre
from typing import Dict


def h_offset(exp_fit: mre.FitResult, expoffset_fit: mre.FitResult) -> bool:
    """
    Test the H_offset hypothesis as described in (Wilting&Priesemann 2018; Supplemental Note 5)

    :return: if True, then the data set should be rejected for m estimation
    """

    return 2*expoffset_fit.ssres < exp_fit.ssres


def h_tau(exp_fit: mre.FitResult, expoffset_fit: mre.FitResult) -> bool:
    """
    Test the H_tau hypothesis as described in (Wilting&Priesemann 2018; Supplemental Note 5)

    :return: if True, then the data set should be rejected for m estimation
    """
    tau_exp, tau_offset = exp_fit.tau, expoffset_fit.tau

    return (np.abs(tau_exp - tau_offset) / np.min((tau_exp, tau_offset))) > 2


def h_lin(rk: mre.CoefficientResult, exp_fit: mre.FitResult) -> bool:
    """
    Test the H_lin hypothesis as described in (Wilting&Priesemann 2018; Supplemental Note 5)

    :return: if True, then the data set should be rejected for m estimation
    """
    lin_fit = mre.fit(rk, fitfunc=mre.f_linear)

    return lin_fit.ssres < exp_fit.ssres


def branching_ratio_figure(bpath):
    nsp = None
    try:
        with open(bpath+'/raw/namespace.p', 'rb') as pfile:
            nsp=pickle.load(pfile)
    except FileNotFoundError:
        print(bpath[-4:], "reports: No namespace data. Skipping.")
        return

    size_factor = 2
    rc('font', size=str(8*size_factor))

    fig: pl.Figure = pl.figure()
    ax_lines, ax_cols = 2,2
    axs: Dict[str, pl.Axes] = {}
    for x,y in itertools.product(range(ax_lines),range(ax_cols)):
        axs['%d,%d'%(x+1,y+1)] = pl.subplot2grid((ax_lines, ax_cols), (x, y))

    fig.set_size_inches(ax_cols*size_factor*6*1.6/2,ax_lines*size_factor*2.*1.6)

    tmin5 =  nsp['T1']+nsp['T2']+nsp['T3']+nsp['T4']

    raster_plot(axs['1,1'], bpath, nsp, tmin=tmin5, tmax=tmin5+nsp["T5"])

    bin_w = 4.7*ms  # BiExp EE time to peak

    rk, ft, (counts, bins) = branching_ratio('', bpath, nsp, bin_w)
    expoffset_fit = mre.fit(rk, fitfunc=mre.f_exponential_offset)
    Hs = h_offset(ft, expoffset_fit), h_tau(ft, expoffset_fit), h_lin(rk, ft)
    H_offset_rej, H_tau_rej, H_lin_rej = ["rej" if H else "acc" for H in Hs]

    tbl = axs['2,1'].table(loc="center",
                           cellText=(("$\Delta{}t$", str(bin_w)),
                                     ("mre", f"{ft.mre:0.6f}"),
                                     ("$H_{offset}$", H_offset_rej),
                                     ("$H_{\\tau}$", H_tau_rej,),
                                     ("$H_{lin}$", H_lin_rej),
                                     ("ATotalMax", nsp['ATotalMax']),
                                     ("iATotalMax", nsp['iATotalMax']),
                                     ("$mu_e$", nsp["mu_e"]),
                                     ("$mu_i$", nsp["mu_i"]),
                                     )
                           )
    tbl.scale(1, 0.6*size_factor)
    axs['2,1'].axis('tight')
    axs['2,1'].axis('off')

    bin_scales = [0.25, 0.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    def fit_bin(bin_w, bin_s):
        rk, _, _ = branching_ratio('', bpath, nsp, bin_w*bin_s)
        return mre.fit(rk, fitfunc=mre.f_exponential_offset)

    fits = [fit_bin(bin_w, bin_s).mre for bin_s in bin_scales]
    fits_expected = [ft.mre**bin_s for bin_s in bin_scales]
    bin_scales.insert(2, 1)
    fits.insert(2, ft.mre)
    fits_expected.insert(2, ft.mre)

    markersize = str(int(5 * size_factor))
    axs['1,2'].plot(bin_scales, fits, label="actual", linestyle="None", marker="D", markersize=markersize)
    axs['1,2'].plot(bin_scales, fits_expected, label="expected", linestyle="None", marker="P", markersize=markersize)
    axs['1,2'].legend()
    axs['1,2'].set_title(f"$actual=m(binsize=k\\Delta{{}}t)=m(binsize=\\Delta{{}}t)^k=expected$, $\Delta{{}}t={bin_w}$")
    axs['1,2'].set_xlabel("$k$")
    axs['1,2'].set_ylabel("$\hat{m}$")

    axs['2,2'].plot(bins[1:], counts)
    axs['2,2'].set_xlim(0, bins[-1])
    axs['2,2'].set_title(f"Activity with binsize={bin_w}")

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
        branching_ratio_figure(bpath)
