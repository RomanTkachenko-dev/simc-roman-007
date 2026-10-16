"""Plot the original SIMC omega model cross section after acceptance cuts.

sigcm is d2sigma/dt/dphi in microbarn/MeV^2/radian. The active model
has no phi dependence, so integration over phi and unit conversion give
dsigma/dt = 2*pi*1e6*sigcm in microbarn/GeV^2. Preserve its sign.
"""
from pathlib import Path
import os
import sys

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

data = np.load(BASE / 'sample_arrays.npz')
accepted = data['pass_cut']
t = data['t'][accepted]  # This source's branch already stores -t.
cross_section = data['sigcm'][accepted] * (2 * np.pi * 1e6)
assert len(t) == 595
assert np.allclose(cross_section, data['dsigma_dt_ub_GeV2'][accepted])

edges = np.arange(.82, 1.10001, .02)
centers, means = [], []
for low, high in zip(edges[:-1], edges[1:]):
    selected = (t >= low) & (t < high)
    if selected.any():
        centers.append((low + high) / 2)
        means.append(cross_section[selected].mean())

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                     'axes.labelsize': 15, 'axes.titlesize': 15,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
fig, ax = plt.subplots(figsize=(9.2, 6.2))
fig.subplots_adjust(left=.15, right=.97, bottom=.14, top=.89)
ax.scatter(t, cross_section, s=17, color='#185d96', alpha=.48,
           edgecolors='none', label='SIMC model values (595 accepted events)')
ax.plot(centers, means, 'o-', color='black', lw=1.5, ms=4,
        label=r'Mean in each 0.02 GeV$^2$ bin')
ax.axhline(0, color='#777777', lw=1)
ax.set(xlabel=r'$-t$ [GeV$^2$]',
       ylabel=r'$d\sigma/dt$ [$\mu$b/GeV$^2$]',
       title=r'SIMC $\omega$ after acceptance cuts',
       xlim=(.82, 1.10), ylim=(-.07, .006))
ax.xaxis.set_major_locator(MultipleLocator(.05))
ax.xaxis.set_minor_locator(MultipleLocator(.01))
ax.yaxis.set_major_locator(MultipleLocator(.01))
ax.yaxis.set_minor_locator(MultipleLocator(.002))
ax.tick_params(which='both', top=True, right=True)
ax.tick_params(which='major', length=6)
ax.tick_params(which='minor', length=3)
ax.legend(loc='upper right', frameon=False, fontsize=10)
ax.grid(which='major', alpha=.13)
fig.savefig(BASE / 'simc_omega_cross_section_after_cuts.png', dpi=180)
fig.savefig(BASE / 'simc_omega_cross_section_after_cuts.svg')
plt.close(fig)
print('Saved standalone SIMC PNG and SVG.')
