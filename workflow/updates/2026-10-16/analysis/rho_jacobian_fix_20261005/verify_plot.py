"""Audit the two requested source edits and the paired rho reruns."""
from pathlib import Path
import ast
import difflib
import hashlib
import json
import os
import re
import struct
import sys
import tarfile

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = ast.parse((BASE.parent / 'sign_check_2021' / 'analyze_sign.py').read_text())
exec(compile(ast.Module(body=[n for n in parser.body if isinstance(n, ast.FunctionDef)], type_ignores=[]), 'validated_ntuple_reader', 'exec'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def cuts(a):
    return ((abs(a['hsdelta']) < 8) & (abs(a['hsxptar']) < .080)
            & (abs(a['hsyptar']) < .035) & (a['ssdelta'] > -10)
            & (a['ssdelta'] < 22) & (abs(a['ssxptar']) < .040)
            & (abs(a['ssyptar']) < .024))

source = Path('//wsl.localhost/Ubuntu/opt/simc/roman007-rho-jacobian-fix-20261005')
archive = Path('D:/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/source_b6025c4.tar.gz')
changed, checked, patches = [], 0, []
with tarfile.open(archive) as tf:
    for m in tf:
        if not m.isfile() or not m.name.endswith(('.f', '.inc')):
            continue
        old = tf.extractfile(m).read()
        new = (source / m.name).read_bytes()
        checked += 1
        if old != new:
            changed.append(m.name)
            old_text, new_text = old.decode(), new.decode()
            if m.name == 'physics_rho_recoil.f':
                expected = old_text.replace('peep_rho=Jttheta_fx*gamma_T*Soding*ntup%sigcm',
                                            'peep_rho=Jttheta*gamma_T*Soding*ntup%sigcm')
            elif m.name == 'NtupleInit.f':
                expected, count = re.subn(r"(?m)^\s*m = m\+1\n[^\n]*NtupleTag\(m\) = 'factor'\s*! 55\n", '', old_text)
                assert count == 1
            else:
                raise AssertionError(m.name)
            assert new_text == expected, m.name
            patches.extend(difflib.unified_diff(old_text.splitlines(True), new_text.splitlines(True),
                                               fromfile='a/' + m.name, tofile='b/' + m.name))
assert set(changed) == {'physics_rho_recoil.f', 'NtupleInit.f'}
(BASE / 'source_changes' / 'requested_fixes.patch').write_text(''.join(patches))

name = 'rho_Roman_007'
original_input = Path(f'D:/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/{name}/{name}.inp')
report = {'base_commit': 'b6025c462ef33eeb620bdf4b399cd602c92e89bf',
          'fixed_source_directory': str(source), 'source_files_checked': checked,
          'changed_source_files': changed, 'executable_sha256': sha(source / 'simc'),
          'input_sha256': sha(original_input), 'input_identical': True,
          'analysis_cuts': 'abs(hsdelta)<8, abs(hsxptar)<0.080, abs(hsyptar)<0.035, -10<ssdelta<22, abs(ssxptar)<0.040, abs(ssyptar)<0.024',
          'runs': {}}
arrays = {}
for variant in ['before', 'after']:
    folder = BASE / variant
    assert (folder / f'{name}.inp').read_bytes() == original_input.read_bytes()
    a = arrays[variant] = read_bin(folder / f'{name}.bin')
    assert len(a['Weight']) == 1000
    mask = cuts(a)
    with uproot.recreate(folder / f'{name}.root') as f:
        tree = f.mktree('h10', {k: 'float64' for k in a})
        tree.extend(a)
        f['schema_note'] = ('Original extra factor header shifts the eepx tail columns.' if variant == 'before'
                            else 'Removed extra factor header: eepx tail labels align with results_write.f.')
    with uproot.open(folder / f'{name}.root') as f:
        for k in a:
            assert np.array_equal(f['h10'][k].array(library='np'), a[k], equal_nan=True), k
    log = (folder / 'simc.log').read_text()
    report['runs'][variant] = {
        'stored_events': len(mask), 'after_cuts': int(mask.sum()),
        'generation_attempts': int(re.search(r'Last Event\s+(\d+)', log)[1]),
        'columns': len(a),
        'all_stored': {k: stats(a[k]) for k in ['sigcm', 'siglab', 'Weight']},
        'after_cuts_stats': {k: stats(a[k][mask]) for k in ['sigcm', 'siglab', 'Weight']},
    }

b, a = arrays['before'], arrays['after']
assert 'factor' in b and 'factor' not in a
assert len(b) == 66 and len(a) == 65
stable = list(b)[:57]
for k in stable:
    if k not in ['siglab', 'Weight']:
        assert np.array_equal(a[k], b[k]), k
assert np.array_equal(cuts(a), cuts(b))
# Header deletion shifts labels back onto the existing data slots.
tail = {'thetamq': 'thetamq', 'phimq': 'factor', 'minus_u': 'phimq',
        'phicmi': 'minus_u', 'tprimei': 'phicmi', 'wcmi': 'tprimei',
        'epsiloni': 'wcmi', 'ui': 'epsiloni'}
for new_key, old_key in tail.items():
    assert np.array_equal(a[new_key], b[old_key]), (new_key, old_key)
assert np.all((a['epsiloni'] > 0) & (a['epsiloni'] < 1))

# Independent Jacobian ratio reconstructed from vertex-level variables.
# This input has a stationary hydrogen target and no energy loss or smearing.
mp = .93827231
proton_p = 1.1 * (1 + a['hsdeltai'] / 100)
proton_e = np.sqrt(proton_p**2 + mp**2)
q2, w, u = a['Q2i'], a['wcmi'], -a['ui']
nu = (w**2 + q2 - mp**2) / (2 * mp)
q_p_costheta = (u + q2 - mp**2) / 2 + nu * proton_e
mass = np.sqrt(w**2 - 2 * mp * proton_e + u + q2)
expected_ratio = (mp + nu - proton_e * q_p_costheta / proton_p**2) / mass
actual_ratio = a['siglab'] / b['siglab']
assert np.allclose(actual_ratio, expected_ratio, rtol=2e-8, atol=1e-10)
assert np.allclose(a['Weight'] / b['Weight'], actual_ratio, rtol=1e-12, atol=1e-12)
report['verification'] = {
    'same_initial_random_state': (BASE / 'before' / f'{name}_start_random_state.dat').read_bytes() == (BASE / 'after' / f'{name}_start_random_state.dat').read_bytes(),
    'same_events_and_sigcm': True, 'tail_column_mapping_verified': tail,
    'jacobian_ratio_max_absolute_difference': float(abs(actual_ratio - expected_ratio).max()),
    'weight_ratio_matches_cross_section_ratio': True,
}
assert report['verification']['same_initial_random_state']
(BASE / 'verification_summary.json').write_text(json.dumps(report, indent=2) + '\n')

fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
for ax, key in zip(axes, ['sigcm', 'siglab', 'Weight']):
    values = np.concatenate([b[key], a[key]])
    lo, hi = min(0., values.min()), max(0., values.max())
    bins = np.linspace(lo - .02 * (hi - lo), hi + .02 * (hi - lo), 46)
    for variant, color, style, label in [('before', '#d45a16', '--', 'Before fix'), ('after', '#125ba6', '-', 'After fix')]:
        ax.hist(arrays[variant][key], bins=bins, histtype='step', linewidth=1.8,
                color=color, linestyle=style, label=label)
    ax.axvline(0, color='.45', linewidth=.8)
    ax.set_xlabel(key + '\n[native SIMC units]')
    ax.set_ylabel('Stored events / bin')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax.grid(alpha=.16)
    ax.legend(frameon=False)
fig.suptitle(r"Recoil $\rho^0$: H(e,e'p)$\rho^0$", fontsize=16)
fig.tight_layout(rect=(0, 0, 1, .95))
for ext in ['png', 'svg']:
    fig.savefig(BASE / f'rho_before_after_fix.{ext}', dpi=180)
plt.close(fig)
print(json.dumps(report, indent=2))
