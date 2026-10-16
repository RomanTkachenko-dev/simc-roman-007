"""Reproduce only Henry's explicitly assumed cross-section central curves.

Does not execute his file-loading/yield-fitting code or invent error bars.
"""
from pathlib import Path
import ast
import hashlib
import json
import os
import sys

OUT = Path(__file__).resolve().parent
SOURCE = Path('D:/EIC/Reference/plot_exp_1over35.py')
sys.path.insert(0, str(OUT.parent / 'omega_cross_section_check/.deps'))
os.environ['MPLCONFIGDIR'] = str(OUT / '.mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tree = ast.parse(SOURCE.read_text(encoding='utf-8-sig'))
keep = []
names = {'W_edges', 't_edges', 'B_slope', 'A_intercepts'}
for node in tree.body:
    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in node.targets):
        keep.append(node)
    if isinstance(node, ast.FunctionDef) and node.name == 'y_theory':
        keep.append(node)
exec(compile(ast.Module(body=keep, type_ignores=[]), str(SOURCE), 'exec'))

t_plot = np.linspace(t_edges[0], t_edges[-1], 300)
t_mid = (t_edges[:-1] + t_edges[1:]) / 2
rows = []
fig, ax = plt.subplots(figsize=(9, 5.7))
colors = ['#0072b2', '#d55e00', '#009e73', '#e69f00']
for i, A in enumerate(A_intercepts):
    y = y_theory(t_plot, A)
    ax.plot(t_plot, y, color=colors[i], linewidth=1.8,
            label=rf'${W_edges[i]:.2f}<W<{W_edges[i+1]:.2f}$ GeV')
    rows.append({'W_GeV': [float(W_edges[i]), float(W_edges[i+1])],
                 'A_ub_per_GeV2': A, 'at_abs_t_0p95_ub_per_GeV2': float(y_theory(.95, A)),
                 't_bin_centers_GeV2': t_mid.tolist(),
                 'assumed_dsigma_dabs_t_ub_per_GeV2': y_theory(t_mid, A).tolist()})
ax.set(xlabel=r'$|t|$ [GeV$^2$]', ylabel=r'$d\sigma/d|t|$ [$\mu$b/GeV$^2$]',
       title=r"Henry's assumed $\omega$ cross-section model", xlim=(t_edges[0], t_edges[-1]), ylim=(0, None))
ax.legend(frameon=False)
ax.text(.03, .04, r'$A\,e^{-5|t|}$ — reference curves, not an extracted cross section',
        transform=ax.transAxes, fontsize=10)
fig.tight_layout()
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'henry_omega_reference.{ext}', dpi=180)
plt.close(fig)
report = {
    'status': 'Assumed reference model reproduced verbatim from Henry supplied plotting script; not a new measurement or updated SIMC model.',
    'source': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'formula': 'dsigma/dabs(t) = A * exp(-B*abs(t))', 'B_GeV_minus2': B_slope,
    'rows': rows,
    'limitations': ['The script assigns cross-section central values with y_theory; fitted omega counts set relative statistical errors only.',
                    'No new error bars: the referenced Analyzed_Pythia_557p034_thin.root is not available in the checked reference/run locations.',
                    'The formula does not specify Q2 or epsilon dependence. Its applicability to our accepted omega events is not established.',
                    'The combined-template-fit amplitudes are not used to normalize these reference curves.',
                    'No cross-section function or omega event weights in SIMC have been replaced.'],
    'needed_for_physical_SIMC_estimate': 'A specified omega model valid at the actual W,Q2,t/u,epsilon and the appropriate generation Jacobian, confirmed by Bill; then rerun and normalize the generated sample.'
}
(OUT / 'henry_omega_reference.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps([{'W': r['W_GeV'], 'dsigma_at_t_0p95': r['at_abs_t_0p95_ub_per_GeV2']} for r in rows], indent=2))
