"""Compare identical-input SIMC runs, retaining signed native output values."""
from pathlib import Path
import ast
import hashlib
import json
import os
import re
import struct
import sys
import tarfile

BASE = Path(__file__).resolve().parent
OLD = BASE.parent / 'sign_check_2021'
PLOTS_ONLY = '--plots-only' in sys.argv
sys.path.insert(0, str(BASE.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Reuse the previously validated binary reader and statistics functions,
# without running that script's output-writing top-level statements.
tree = ast.parse((OLD / 'analyze_sign.py').read_text())
functions = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)], type_ignores=[])
exec(compile(functions, str(OLD / 'analyze_sign.py'), 'exec'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

report = {
    'versions': {'2021': '8b43f9d0dd3be977df93b35abc528a389ca70d5c',
                 'Alicia_Acpostuma': 'b6025c462ef33eeb620bdf4b399cd602c92e89bf'},
    'selection': 'All 1000 stored events per run; no additional analysis cuts.',
    'values': 'Signed native ntuple values; no event weighting, normalization, absolute value, or unit conversion.',
    'inputs': {}, 'tests': {},
}
for name in ['test_eep_h_root.inp', 'extra_hydrogen_root.inp', 'test_pion_h_sign.inp', 'extra_nucpi_sign.inp']:
    assert (OLD / name).read_bytes() == (BASE / name).read_bytes()
    report['inputs'][name] = {'byte_identical': True, 'sha256': sha(BASE / name)}

source = Path('//wsl.localhost/Ubuntu/opt/simc/roman007-b6025c4')
archive = Path('D:/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/source_b6025c4.tar.gz')
checked, differences = 0, []
with tarfile.open(archive) as tf:
    for m in tf:
        if m.isfile() and m.name.endswith(('.f', '.inc')):
            local = source / m.name
            if tf.extractfile(m).read() != local.read_bytes():
                differences.append(m.name)
            checked += 1
report['source_audit'] = {'checked_fortran_and_include_files': checked,
                          'differences_from_archive': differences,
                          'source_archive_sha256': sha(archive), 'executable_sha256': sha(source / 'simc')}
assert not differences, differences

for name, keys, title in [
    ('test_eep_h_root', ['sigcc', 'Weight'], "Elastic H(e,e'p)"),
    ('test_pion_h_sign', ['sigcm', 'siglab', 'Weight'], r"Pion H(e,e'$\pi^+$)n"),
]:
    arrays = {'2021': read_bin(OLD / f'{name}.bin'), 'Alicia': read_bin(BASE / f'{name}.bin')}
    print(name, 'Alicia columns:', list(arrays['Alicia']))
    for data in arrays.values():
        assert len(data['Weight']) == 1000
    branch_map = {k: k for k in keys}
    if name == 'test_eep_h_root':
        # Both results_write.f versions store main%sigcc, but the 2026
        # version routes doing_eep into the meson-style ntuple layout:
        # 2021: slot 43, sigcc; 2026: slot 45, siglab.
        branch_map['sigcc'] = 'siglab'
    assert all(k in arrays['2021'] and branch_map[k] in arrays['Alicia'] for k in keys)
    if not PLOTS_ONLY:
        with uproot.recreate(BASE / f'{name}.root') as f:
            t = f.mktree('h10', {k: 'float64' for k in arrays['Alicia']})
            t.extend(arrays['Alicia'])
    with uproot.open(BASE / f'{name}.root') as f:
        for k in keys:
            assert np.array_equal(f['h10'][branch_map[k]].array(library='np'), arrays['Alicia'][branch_map[k]])
    if name == 'test_eep_h_root':
        arrays['Alicia']['sigcc'] = arrays['Alicia']['siglab']
    test = {'reaction': title, 'runs': {}, 'alicia_native_branch_mapping': branch_map,
            'mean_ratios_alicia_over_2021': {},
            'initial_random_states_byte_identical': (OLD / f'{name}_start_random_state.dat').read_bytes() == (BASE / f'{name}_start_random_state.dat').read_bytes()}
    for label, folder in [('2021', OLD), ('Alicia', BASE)]:
        log = (folder / f'{name}.log').read_text()
        test['runs'][label] = {'stored_events': 1000,
            'generation_attempts': int(re.search(r'Last Event\s+(\d+)', log)[1]),
            'signs': {k: stats(arrays[label][k]) for k in keys},
            'binary_sha256': sha(folder / f'{name}.bin'),
            'initial_random_state_sha256': sha(folder / f'{name}_start_random_state.dat')}
    for k in keys:
        test['mean_ratios_alicia_over_2021'][k] = float(arrays['Alicia'][k].mean() / arrays['2021'][k].mean())
    common = arrays['2021'].keys() & arrays['Alicia'].keys()
    test['same_index_exactly_equal_branches'] = sorted(k for k in common if np.array_equal(arrays['2021'][k], arrays['Alicia'][k]))
    report['tests'][name] = test
    fig, axes = plt.subplots(1, len(keys), figsize=(6 * len(keys), 4.8))
    for ax, key in zip(np.atleast_1d(axes), keys):
        values = np.concatenate([arrays[label][key] for label in arrays])
        lo, hi = min(0., float(values.min())), float(values.max())
        bins = np.linspace(lo - .025 * (hi - lo), hi + .025 * (hi - lo), 46)
        for label, color, style in [('2021', '#125ba6', '-'), ('Alicia', '#d45a16', '--')]:
            ax.hist(arrays[label][key], bins=bins, histtype='step', linewidth=1.8, color=color, linestyle=style, label=label)
        ax.axvline(0, color='0.4', lw=.8)
        display_key = 'Cross section' if name == 'test_eep_h_root' and key == 'sigcc' else key
        ax.set_xlabel(f'{display_key}\n[native SIMC units]')
        ax.set_ylabel('Stored events / bin')
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        ax.legend(frameon=False)
        ax.grid(alpha=.18)
    fig.suptitle(title, fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, .95))
    outname = 'elastic_comparison' if name == 'test_eep_h_root' else 'pion_comparison'
    fig.savefig(BASE / f'{outname}.png', dpi=180)
    fig.savefig(BASE / f'{outname}.svg')
    plt.close(fig)

report['scope'] = 'Elastic and pion comparisons only. These tests do not validate the omega model or identify the reported omega sign bug.'
if not PLOTS_ONLY:
    (BASE / 'comparison_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
