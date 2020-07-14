from analysis.axes.network import branching_ratio
from analysis.multisim.multisimplotrunner import *

from brian2.units import *
import numpy as np
import mrestimator as mre


class SrvPrbBranching(MultiSimPlotRunner):

    bin_w = 4.7*ms
    
    def __init__(self, name="srvprb_branching", plot_count=(1, 1)):
        super(SrvPrbBranching, self).__init__(name, plot_count)

    def plot(self, directories, nsps, fig, axs):
        ax = axs[0][0]
        fontsize = 18
        revreg, subcrit = "reverberating", "sub-critical"
        colors = {revreg: "green", subcrit: "blue"}
        data, mres = {revreg: [], subcrit: []}, {revreg: [], subcrit: []}
        for dir, nsp in zip(directories, nsps):
            # we assume that only directories with data have been selected
            # that passes the stationary tests
            rk, ft, _ = branching_ratio('', dir, nsp, self.bin_w)
            fit = mre.fit(rk, fitfunc=mre.f_exponential_offset)
            if fit.mre > 0.995:
                print(f"skipping: mre is critical or super-critical: {dir}")
                continue

            group = revreg if fit.mre > 0.9 else subcrit
            data[group].append(self.unpickle(dir, "survival_full_t"))
            mres[group].append(fit.mre)

        assert(len(data[revreg]) > 0 and len(data[subcrit]) > 0)

        # strct plasticity is only applied every 1000ms
        # so bin size needs to be at least 1s, which logspace doesn't
        # do for the range <10e1
        bins = np.logspace(np.log10(10e1),
                           np.log10(np.max([data[revreg][i]['t_split']/second for i in range(0, len(data[revreg]))])+0.1),
                           num=70)  # 82 exhibits jump in blue curve
        bins = np.hstack((np.arange(1, 10, 1), bins))
        centers = (bins[:-1] + bins[1:])/2.

        ax.loglog()

        for label, dfs in data.items():
            color = colors[label]
            counts = [np.histogram(df['full_t'], bins=bins, density=True)[0] for df in dfs]
            avg, std = np.mean(counts, axis=0), np.std(counts, axis=0)
            ax.plot(centers, avg, label=label, color=color)
            ax.fill_between(centers, avg-std, avg+std, facecolor=f"light{color}", linewidth=0)

        group_info = [(label, len(dfs), np.mean(mres[label]), np.std(mres[label])) for label, dfs in data.items()]
        group_info_str = [f"{label}: N={no:2d}, $\hat{{m}}={mean:.2f}\pm{std:.2f}$" for label, no, mean, std in group_info]
        ax.legend(loc="lower left", fontsize=fontsize)
        # ax.text(1.0, -0.25, '\n'.join(group_info_str), transform=ax.transAxes, ha='right')
        ax.set_xlabel("survival time [s]", fontsize=fontsize)
        ax.set_ylabel("probaiblity density", fontsize=fontsize)
        ax.set_title("Survival times of EE synapses", fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)


if __name__ == '__main__':
    SrvPrbBranching().run()