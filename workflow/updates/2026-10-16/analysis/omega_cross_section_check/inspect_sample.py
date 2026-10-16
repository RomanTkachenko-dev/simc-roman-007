"""Read the archived omega events; preserve signs and original run files."""
from pathlib import Path
import sys, json, hashlib, math

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / '.deps'))
import numpy as np
import uproot

RUN = Path('D:/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/omega_Roman_007')
ROOT = RUN / 'omega_Roman_007.root'
with uproot.open(ROOT) as f:
    a = f['h10'].arrays(library='np')
    a = {k: np.asarray(v, dtype=float) for k, v in a.items()}

cut = ((abs(a['hsdelta']) < 8) & (abs(a['hsxptar']) < .080)
       & (abs(a['hsyptar']) < .035) & (a['ssdelta'] > -10)
       & (a['ssdelta'] < 22) & (abs(a['ssxptar']) < .040)
       & (abs(a['ssyptar']) < .024))
assert cut.sum() == 595
assert math.isclose(a['Weight'][cut].sum(), 5.12814492630298e-5, rel_tol=1e-12)
a['dsigma_dt_ub_GeV2'] = a['sigcm'] * (2 * np.pi * 1e6)
mp = .93827231
# The frozen NtupleInit.f includes an extra 'factor' label after thetamq;
# results_write.f does not write that extra field. Match column positions,
# not misleading labels, in this archived output. Source files are unchanged.
a['model_minus_u_GeV2'] = a['epsiloni']
a['model_W_GeV'] = a['tprimei']
a['model_epsilon'] = a['wcmi']
a['omega_theta_cm_deg'] = 180 - np.degrees(a['thetacmi'])
u = abs(a['model_minus_u_GeV2'])
q2 = a['Q2i']
w = a['model_W_GeV']
eps = a['model_epsilon']
a['sigma_T_model'] = (20.16527 - 4.2124 * u) / q2**2 / (w**2 - mp**2)**2
a['sigma_L_model'] = (12.2546 - 29.8629 * u) / q2**4 / (w**2 - mp**2)**2
a['model_from_branches'] = a['sigma_T_model'] + eps * a['sigma_L_model']
a['siglab_over_sigcm'] = a['siglab'] / a['sigcm']

def stats(x, weights):
    return {'min': float(x.min()), 'max': float(x.max()), 'mean': float(x.mean()),
            'median': float(np.median(x)), 'weighted_mean': float(np.average(x, weights=weights)),
            'negative': int((x < 0).sum()), 'positive': int((x > 0).sum()),
            'q05': float(np.quantile(x, .05)), 'q95': float(np.quantile(x, .95))}

keys = ['Q2','Q2i','W','Wi','wcmi','t','ti','ui','minus_u','epsilon','epsiloni',
        'thetacmi','sigcm','siglab','Weight','dsigma_dt_ub_GeV2',
        'model_minus_u_GeV2','model_W_GeV','model_epsilon','omega_theta_cm_deg',
        'sigma_T_model','sigma_L_model','model_from_branches','siglab_over_sigcm']
report = {'source_file': str(ROOT), 'sha256': hashlib.sha256(ROOT.read_bytes()).hexdigest(),
          'stored_events': len(cut), 'after_cuts': int(cut.sum()),
          'negative_sigcm_before_cuts': int((a['sigcm'] < 0).sum()),
          'cut': 'abs(hsdelta)<8 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035 && -10<ssdelta<22 && abs(ssxptar)<0.040 && abs(ssyptar)<0.024',
          'conversion': 'dsigma_dt_ub_GeV2 = 2*pi*1e6*sigcm; phi independent active model; no abs()',
          'schema_warning': {'reason': 'Extra factor name at NtupleInit.f:178-179 has no matching entry at results_write.f:183-190',
                             'actual_model_W_GeV': 'stored label tprimei',
                             'actual_model_epsilon': 'stored label wcmi',
                             'actual_model_minus_u_GeV2': 'stored label epsiloni',
                             'unaffected': ['sigcm','siglab','Weight','Q2','Q2i','W','Wi','t','ti','thetacmi']},
          'branches': {k: stats(a[k][cut], a['Weight'][cut]) for k in keys},
          'model_reproduction': {'max_abs_difference': float(np.max(abs(a['model_from_branches'] - a['dsigma_dt_ub_GeV2']))),
            'max_relative_difference': float(np.max(abs((a['model_from_branches'] - a['dsigma_dt_ub_GeV2'])/a['dsigma_dt_ub_GeV2'])))}}
assert np.allclose(a['model_from_branches'], a['dsigma_dt_ub_GeV2'], rtol=2e-4, atol=2e-7), report['model_reproduction']
(BASE/'sample_summary.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
np.savez(BASE/'sample_arrays.npz', **a, pass_cut=cut)
np.savetxt(BASE/'accepted_events.csv', np.column_stack([a[k][cut] for k in keys]),
           delimiter=',',header=','.join(keys),comments='')
print(json.dumps(report, indent=2))
