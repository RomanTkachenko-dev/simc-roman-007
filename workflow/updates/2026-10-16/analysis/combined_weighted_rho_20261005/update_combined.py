"""Replace the old count-filled rho template with corrected event weights.

Reproduce the previous fixed-shape, nonnegative seven-amplitude fit first.
All output scales remain diagnostic and must not be interpreted as cross sections.
"""
from pathlib import Path
import ast
import hashlib
import itertools
import json
import os
import sys

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(OUT / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from uproot.writing.identify import to_TH1x, to_TAxis

BASE = Path('/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN')
PRIOR = BASE / 'meeting_20260916_template_fit'
RHO = OUT.parent / 'rho_jacobian_fix_20261005/after/rho_Roman_007.root'
CHANNELS = ['omega', 'phi', 'rho', 'eta', 'etap', 'xphsp', 'xphsp3']
WEIGHTED = {'omega', 'rho', 'xphsp', 'xphsp3'}
COLORS = ['#0000ff', '#e07000', '#00b000', '#ff8000', '#a040ff', '#707070', '#00afc1']
LABELS = [r'$\omega$: Weight (provisional)', r'$\phi$: counts', r'$\rho$: Weight',
          r'$\eta$: counts', r"$\eta'$: counts", r'$2\pi$: Weight', r'$3\pi$: Weight']
module = ast.parse((OUT.parent / 'omega_W247_Q235/analyze.py').read_text())
exec(compile(ast.Module(body=[n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'cut'], type_ignores=[]), 'existing_cuts', 'exec'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def solve_nnls(design, target):
    # Seven parameters allow exhaustive enumeration of all 128 possible
    # active sets. Each feasible least-squares solution is a candidate.
    best = None
    for bits in itertools.product([False, True], repeat=design.shape[1]):
        active = np.array(bits)
        x = np.zeros(design.shape[1])
        if active.any():
            trial, _, _, _ = np.linalg.lstsq(design[:, active], target, rcond=None)
            if np.any(trial < 0):
                continue
            x[active] = trial
        residual = design @ x - target
        score = float(residual @ residual)
        if best is None or score < best[0]:
            best = score, x
    score, x = best
    gradient = design.T @ (design @ x - target)
    assert np.max(np.abs(gradient[x > 0]), initial=0) < 1e-5
    assert np.all(gradient[x == 0] >= -1e-5)
    return x, score

old_report = json.loads((PRIOR / 'summary.json').read_text())
weight_meta = json.loads((BASE / 'meeting_20260916_weight_applied/summary.json').read_text())
original_meta = json.loads((BASE / 'meeting_20260916_overlay/summary.json').read_text())
sources = [PRIOR / 'results.root', PRIOR / 'summary.json', RHO,
           BASE / 'meeting_20260916_weight_applied/summary.json']
hashes = {str(p): sha(p) for p in sources}
raw, variances = {}, {}
with uproot.open(PRIOR / 'results.root') as f:
    edges = f['pythia_counts_after_cuts'].axis().edges()
    counts = f['pythia_counts_after_cuts'].values()
    pyflow = f['pythia_counts_after_cuts'].values(flow=True)
    display_py = f['pythia_display_unchanged'].values()
    prior_sum = f['fitted_sum_display'].values()
    for ch in CHANNELS:
        raw[ch] = f[ch + '_input_after_cuts'].values(flow=True)
        variances[ch] = f[ch + '_input_after_cuts'].variances(flow=True)
assert np.array_equal(edges, np.linspace(0, 2, 101))
display = weight_meta['common_display_scale_retained']
np.testing.assert_allclose(display_py, counts * display, rtol=1e-12)

def refill(a, weighted):
    mask = cut(a)
    x = a['missmass'][mask]
    w = a['Weight'][mask].astype(float) if weighted else np.ones(mask.sum())
    # ROOT bins: values equal to xmax belong to overflow.
    index = np.searchsorted(edges, x, side='right')
    h = np.bincount(index, weights=w, minlength=102)
    v = np.bincount(index, weights=w*w, minlength=102)
    return h, v, int(mask.sum())

# Independently reproduce all old raw templates from their source events.
for ch in CHANNELS:
    path = BASE / 'runs' / (ch + '_Roman_007') / (ch + '_Roman_007.root')
    hashes[str(path)] = sha(path)
    with uproot.open(path) as f:
        data = f['h10'].arrays(library='np')
    h, v, naccepted = refill(data, ch in WEIGHTED - {'rho'})
    np.testing.assert_allclose(h, raw[ch], rtol=1e-10, atol=1e-28)
    np.testing.assert_allclose(v, variances[ch], rtol=1e-10, atol=1e-28)
    assert naccepted == original_meta['channels'][ch]['after_cuts']

total_counts = counts.sum()
sigma = np.sqrt(np.maximum(counts, 1.))
old_matrix = np.column_stack([raw[ch][1:-1] / raw[ch][1:-1].sum() for ch in CHANNELS])
old_design = old_matrix * total_counts / sigma[:, None]
target = counts / sigma
old_coef, old_chi2 = solve_nnls(old_design, target)
expected_coef = np.array([old_report['step2']['coefficients'][ch]['fraction_of_Pythia_count_total'] for ch in CHANNELS])
np.testing.assert_allclose(old_coef, expected_coef, rtol=1e-9, atol=1e-12)
np.testing.assert_allclose(old_chi2, old_report['step2']['chi2_Pythia_variance_only'], rtol=1e-12)
np.testing.assert_allclose(old_matrix @ old_coef * total_counts * display, prior_sum, rtol=1e-10, atol=1e-12)

with uproot.open(RHO) as f:
    rho = f['h10'].arrays(library='np')
rho_cut = cut(rho)
assert np.all(rho['Weight'][rho_cut] > 0)
rho_counts, _, accepted_rho = refill(rho, False)
np.testing.assert_array_equal(rho_counts, raw['rho'])
old_rho_raw = raw['rho'].copy()
raw['rho'], variances['rho'], _ = refill(rho, True)
matrix = np.column_stack([raw[ch][1:-1] / raw[ch][1:-1].sum() for ch in CHANNELS])
design = matrix * total_counts / sigma[:, None]
coef, chi2 = solve_nnls(design, target)
# Solver closure using seven strictly positive known coefficients.
truth = np.array([.12, .03, .15, .01, .02, .27, .4])
recovered, _ = solve_nnls(design, design @ truth)
np.testing.assert_allclose(recovered, truth, atol=1e-12)
parts = matrix * coef[None, :] * total_counts * display
summed = parts.sum(axis=1)
residual = display_py - summed
np.testing.assert_allclose(summed, matrix @ coef * total_counts * display, rtol=1e-12)
active = coef > 0
covariance = np.zeros((7, 7))
covariance[np.ix_(active, active)] = np.linalg.inv(design[:, active].T @ design[:, active])
coefficient_report = {}
for j, ch in enumerate(CHANNELS):
    area = raw[ch][1:-1].sum()
    coefficient_report[ch] = {
        'fill': 'Weight' if ch in WEIGHTED else 'counts',
        'template_visible_integral': float(area),
        'fraction_of_Pythia_total': float(coef[j]),
        'Pythia_count_equivalent': float(coef[j] * total_counts),
        'previous_Pythia_count_equivalent': float(old_coef[j] * total_counts),
        'scale_raw_to_Pythia_counts': float(coef[j] * total_counts / area),
        'at_zero_boundary': bool(not active[j]),
        'conditional_count_error_fixed_templates': float(np.sqrt(covariance[j, j]) * total_counts) if active[j] else None,
    }
report = {
    'status': 'Provisional seven-channel template fit; not an absolute yield or cross-section extraction.',
    'changed_input': 'rho count template replaced by corrected rho event-Weight template; remaining six raw templates unchanged.',
    'rho_source': str(RHO), 'rho_accepted_events': accepted_rho,
    'rho_sum_weight_after_cuts': float(raw['rho'].sum()),
    'cuts': weight_meta['cuts'], 'range_GeV': [0, 2], 'bins': 100,
    'fit': 'Nonnegative seven free amplitudes, fixed histogram shapes, observed-count variance max(n,1).',
    'solver': 'Exhaustive 128 active-set least-squares enumeration; original NNLS/ROOT result reproduced; KKT and synthetic closure checked.',
    'Pythia_visible_counts': float(total_counts), 'Pythia_underflow': float(pyflow[0]), 'Pythia_overflow': float(pyflow[-1]),
    'display_scale_unchanged': display, 'previous_chi2': old_chi2, 'updated_chi2': chi2,
    'nominal_ndf': 93, 'updated_chi2_over_nominal_ndf': chi2/93,
    'coefficients': coefficient_report,
    'caveats': ['omega retains its previous Weight-filled shape, but its underlying sigcm is negative after cuts; omega component remains provisional.',
                'phi, eta and eta-prime remain count-filled templates.',
                'Fixed-template fit ignores finite Monte Carlo template uncertainties; the fit is not a validated background subtraction.',
                weight_meta['pythia_variable_caveat'],
                'Fitted amplitudes are Pythia count equivalents, not luminosity-normalized physical yields or cross sections.',
                'No updated omega model has been specified yet.'],
    'sources_sha256': hashes,
}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
fig, (ax, resax) = plt.subplots(2, 1, figsize=(12, 8.2), sharex=True,
                              gridspec_kw={'height_ratios': [3.8, 1], 'hspace': .06})
for j, (color, label) in enumerate(zip(COLORS, LABELS)):
    ax.stairs(parts[:, j], edges, color=color, linewidth=1.2, label=label)
ax.stairs(summed, edges, color='red', linewidth=1.8, label='Template fit (provisional)')
ax.stairs(display_py, edges, color='black', linewidth=1.8, label='Pythia Background')
ax.set_ylabel('Pythia-referenced ordinate [a.u.]')
ax.set_ylim(0, max(summed.max(), display_py.max()) * 1.12)
ax.legend(loc='upper right', frameon=False, fontsize=10)
resax.stairs(residual, edges, color='black', linewidth=1.2)
resax.axhline(0, color='.4', linestyle=':', linewidth=.8)
resax.set_ylabel('Pythia - fit\n[a.u.]')
resax.set_xlabel('Missing mass [GeV]')
resax.set_xlim(0, 2)
for panel in [ax, resax]:
    panel.minorticks_on()
    panel.tick_params(which='both', top=True, right=True)
fig.subplots_adjust(left=.10, right=.98, top=.98, bottom=.09)
for ext in ['png', 'svg']:
    fig.savefig(OUT / f'combined_signal_weighted_rho.{ext}', dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 5))
ax.stairs(old_matrix[:, 2], edges, label='Previous rho: counts', color='#777777', linestyle='--', linewidth=1.7)
ax.stairs(matrix[:, 2], edges, label='Corrected rho: Weight', color='#00a050', linewidth=1.7)
ax.set(xlabel='Missing mass [GeV]', ylabel='Fraction per bin', xlim=(0, 2))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / 'rho_template_shape_comparison.png', dpi=180)
plt.close(fig)

def root_hist(name, contents, variance, entries):
    centers = (edges[:-1] + edges[1:]) / 2
    return to_TH1x(name, '', contents.astype('float64'), float(entries),
                   float(contents[1:-1].sum()), float(variance[1:-1].sum()),
                   float(contents[1:-1] @ centers), float(contents[1:-1] @ centers**2),
                   variance.astype('float64'), to_TAxis('xaxis', 'Missing mass [GeV]', 100, 0., 2.))

with uproot.recreate(OUT / 'results.root') as f:
    for j, ch in enumerate(CHANNELS):
        f[ch + '_input_after_cuts'] = root_hist(ch, raw[ch], variances[ch], original_meta['channels'][ch]['after_cuts'])
        f[ch + '_fitted_display'] = (parts[:, j], edges)
    f['fitted_sum_display'] = (summed, edges)
    f['pythia_display_unchanged'] = (display_py, edges)
    f['pythia_counts_after_cuts'] = root_hist('Pythia', pyflow, pyflow, pyflow.sum())
    f['residual_display'] = (residual, edges)
    f['previous_fitted_sum_display'] = (prior_sum, edges)
    f['metadata_json'] = json.dumps(report)
with uproot.open(OUT / 'results.root') as f:
    np.testing.assert_allclose(sum(f[ch + '_fitted_display'].values() for ch in CHANNELS), f['fitted_sum_display'].values(), rtol=1e-12)
    np.testing.assert_array_equal(f['rho_input_after_cuts'].values(flow=True), raw['rho'])
    np.testing.assert_array_equal(f['rho_input_after_cuts'].variances(flow=True), variances['rho'])
for p, digest in hashes.items():
    assert sha(Path(p)) == digest
(OUT / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'previous_chi2': old_chi2, 'new_chi2': chi2,
                  'rho_accepted': accepted_rho, 'coefficients': coefficient_report}, indent=2))
