from analysis.axes.network import branching_ratio
from analysis.multisim.multisimplotrunner import *

from brian2.units import *
import numpy as np
import mrestimator as mre


class SrvPrbBranching(MultiSimPlotRunner):

    bin_w = 4.7*ms
    
    def __init__(self):
        super(SrvPrbBranching, self).__init__(
            name="srvprb_branching", plot_count=(1, 1)
        )

    def plot(self, directories, nsps, axs):
        ax = axs[0][0]
        revreg, subcrit = "reverberating regime ($0.9 < m < 0.995$)", "sub-critical ($m <= 0.9$)"
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

        bins = np.logspace(np.log10(1),
                           np.log10(data[revreg][0]['t_split'] / second + 0.1),
                           num=100)
        centers = (bins[:-1] + bins[1:])/2.

        for label, dfs in data.items():
            color = colors[label]
            counts = np.array([np.histogram(df['full_t'], bins=bins, density=True)[0] for df in dfs])
            avg, std = np.mean(counts, axis=0), np.std(counts, axis=0)
            ax.plot(centers, avg, label=label, color=color)
            ax.fill_between(centers, avg-std, avg+std, facecolor=f"light{color}", linewidth=0)

        group_info = [(label, len(dfs), np.mean(mres[label]), np.std(mres[label])) for label, dfs in data.items()]
        group_info_str = [f"{label}: {no:2d}, mean {mean:.2f}, std {std:.2f}" for label, no, mean, std in group_info]
        ax.legend()
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.text(1.0, -0.25, '\n'.join(group_info_str), transform=ax.transAxes, ha='right')
        ax.set_xlabel("survival time [s]")
        ax.set_ylabel("density")
        ax.set_title("EE synapse survival by branching factor")


if __name__ == '__main__':
    SrvPrbBranching().run()