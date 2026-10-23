"""Check omega after the rho-only Jacobian and shared ntuple-header fixes."""
from pathlib import Path
import ast
import hashlib
import json
import os
import re
import struct
import sys

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

for path, selected in [(BASE.parent / 'sign_check_2021' / 'analyze_sign.py', ['read_bin', 'stats']),
                       (BASE.parent / 'omega_W247_Q235' / 'analyze.py', ['cut'])]:
    module = ast.parse(path.read_text())
    exec(compile(ast.Module(body=[n for n in module.body if isinstance(n, ast.FunctionDef) and n.name in selected], type_ignores=[]), str(path), 'exec'))

name = 'omega_Roman_007'
report = {'scope': 'Original omega input rerun before and after the rho Jacobian and shared NtupleInit fixes. physics_omega.f unchanged; it still uses Jttheta_fx.',
          'runs': {}}
arrays = {}
for variant in ['before', 'after']:
    folder = BASE / variant
    a = arrays[variant] = read_bin(folder / f'{name}.bin')
    selection = cut(a)
    with uproot.recreate(folder / f'{name}.root') as f:
        t = f.mktree('h10', {k: 'float64' for k in a})
        t.extend(a)
        f['schema_note'] = ('Original extra factor header shifts eepx tail labels.' if variant == 'before'
                            else 'Corrected eepx tail labels; omega model and Jacobian unchanged.')
    report['runs'][variant] = {
        'stored': len(selection), 'after_cuts': int(selection.sum()),
        'attempts': int(re.search(r'Last Event\s+(\d+)', (folder / 'simc.log').read_text())[1]),
        'columns': len(a), 'input_sha256': hashlib.sha256((folder / f'{name}.inp').read_bytes()).hexdigest(),
        'all_stored': {k: stats(a[k]) for k in ['sigcm', 'siglab', 'Weight']},
        'after_cuts_stats': {k: stats(a[k][selection]) for k in ['sigcm', 'siglab', 'Weight']},
    }
b, a = arrays['before'], arrays['after']
assert report['runs']['before']['input_sha256'] == report['runs']['after']['input_sha256']
assert len(a) == 65 and len(b) == 66 and 'factor' not in a
for k in list(a)[:58]:
    assert np.array_equal(a[k], b[k]), k
tail = {'phimq': 'factor', 'minus_u': 'phimq', 'phicmi': 'minus_u',
        'tprimei': 'phicmi', 'wcmi': 'tprimei', 'epsiloni': 'wcmi', 'ui': 'epsiloni'}
for new_key, old_key in tail.items():
    assert np.array_equal(a[new_key], b[old_key]), new_key
assert np.array_equal(cut(a), cut(b))
report['verification'] = {'cross_sections_weights_and_kinematics_exactly_unchanged': True,
                          'corrected_tail_labels_verified': True}
(BASE / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')

fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
for ax, key in zip(axes, ['sigcm', 'siglab', 'Weight']):
    values = a[key]
    lo, hi = min(0., values.min()), max(0., values.max())
    bins = np.linspace(lo - .02 * (hi - lo), hi + .02 * (hi - lo), 46)
    ax.hist(a[key], bins=bins, histtype='step', linewidth=2, color='#125ba6', label='After rho/header fixes')
    ax.hist(b[key], bins=bins, histtype='step', linewidth=1.5, color='#d45a16', linestyle='--', label='Original')
    ax.axvline(0, color='.45', linewidth=.8)
    ax.set_xlabel(key + '\n[native SIMC units]')
    ax.set_ylabel('Stored events / bin')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax.grid(alpha=.16)
    ax.legend(frameon=False, fontsize=9)
fig.suptitle(r"$\omega$: H(e,e'p)$\omega$", fontsize=16)
fig.tight_layout(rect=(0, 0, 1, .95))
for ext in ['png', 'svg']:
    fig.savefig(BASE / f'omega_after_rho_fixes.{ext}', dpi=180)
plt.close(fig)
print(json.dumps(report, indent=2))
