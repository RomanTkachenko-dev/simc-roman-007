"""Diagnostic plots of signed SIMC values; all stored events, no weighting."""
from pathlib import Path
import sys, os
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.labelsize': 12, 'axes.titlesize': 12,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
BLUE, GREEN = '#185d96', '#21845a'

def load(name):
    with uproot.open(BASE / f'{name}.root') as f:
        a = f['h10'].arrays(library='np')
    assert len(a['Weight']) == 1000
    return a

def style(ax, axis):
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-3, 3))
    ax.ticklabel_format(axis=axis, style='sci', scilimits=(-3, 3), useMathText=True)
    getattr(ax, f'{axis}axis').set_major_formatter(fmt)
    ax.tick_params(top=True, right=True)
    ax.grid(alpha=.13)

def hist(ax, x, name, color):
    assert np.isfinite(x).all()
    vmax = float(x.max())
    edges = np.linspace(-.12*vmax, 1.06*vmax, 51)
    values, _ = np.histogram(x, edges)
    assert values.sum() == 1000
    ax.stairs(values, edges, color=color, linewidth=1.5)
    ax.axvline(0, color='#777777', linestyle='--', linewidth=1)
    ax.set(xlabel=f'Stored {name} value', ylabel='Events / bin',
           title=f'{name} distribution', xlim=(edges[0], edges[-1]))
    ax.text(.97, .94, f'Positive: {(x>0).sum()}\nNegative: {(x<0).sum()}',
            transform=ax.transAxes, ha='right', va='top', fontsize=10)
    style(ax, 'x')

def scatter(ax, x, y, xlabel, name, color):
    ax.scatter(x, y, s=9, alpha=.42, color=color, edgecolors='none')
    ax.axhline(0, color='#777777', linestyle='--', linewidth=1)
    ax.set(xlabel=xlabel, ylabel=f'Stored {name} value', title=f'{name} vs kinematics',
           ylim=(-.08*float(y.max()), 1.10*float(y.max())))
    ax.text(.97, .95, f'Negative: {(y<0).sum()} / {len(y)}',
            transform=ax.transAxes, ha='right', va='top', fontsize=10)
    style(ax, 'y')

def save(fig, filename):
    fig.text(.5, .025, 'All 1,000 stored events; no additional analysis cuts. Raw signed branch values.',
             ha='center', fontsize=10, color='#444444')
    for ext in ['png', 'svg']:
        fig.savefig(BASE/f'{filename}.{ext}', dpi=180)
    plt.close(fig)

e = load('test_eep_h_root')
fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
fig.subplots_adjust(left=.075, right=.985, bottom=.20, top=.76, wspace=.34)
fig.suptitle(r'Bill SIMC master (2021): elastic $H(e,e^{\prime}p)$', y=.94, fontsize=16)
scatter(ax[0], e['Q2'], e['sigcc'], r'$Q^2$ [GeV$^2$]', 'sigcc', BLUE)
hist(ax[1], e['sigcc'], 'sigcc', BLUE)
hist(ax[2], e['Weight'], 'Weight', GREEN)
save(fig, 'bill_2021_elastic_sign_check')

p = load('test_pion_h_sign')
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
fig.subplots_adjust(left=.105, right=.98, bottom=.105, top=.85, hspace=.56, wspace=.32)
fig.suptitle(r'Bill SIMC master (2021): $H(e,e^{\prime}\pi^+)n$', y=.95, fontsize=16)
scatter(ax[0,0], p['t'], p['sigcm'], r'$-t$ [GeV$^2$]', 'sigcm', BLUE)
hist(ax[0,1], p['sigcm'], 'sigcm', BLUE)
hist(ax[1,0], p['siglab'], 'siglab', BLUE)
hist(ax[1,1], p['Weight'], 'Weight', GREEN)
save(fig, 'bill_2021_pion_sign_check')
print('Saved elastic and pion sign-check PNG/SVG figures.')
