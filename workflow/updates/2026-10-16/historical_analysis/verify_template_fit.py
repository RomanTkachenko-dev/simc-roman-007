"""Read-only round-trip checks for the generated ROOT artifact."""
from pathlib import Path
import json
import math
import ROOT
base = Path(__file__).resolve().parent/'meeting_20260916_template_fit'
report = json.loads((base/'summary.json').read_text())
f = ROOT.TFile.Open(str(base/'results.root'), 'READ')
assert f and not f.IsZombie() and not f.TestBit(ROOT.TFile.kRecovered)
channels = ['omega','phi','rho','eta','etap','xphsp','xphsp3']
for suffix, total in [('_phase_boost_step1','phase_boost_sum'),('_fitted_display','fitted_sum_display')]:
    for i in range(102):
        expected = math.fsum(f.Get(ch+suffix).GetBinContent(i) for ch in channels)
        actual = f.Get(total).GetBinContent(i)
        assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-25)
assert f.Get('ROOT_fit_result').Status() == 0
assert json.loads(f.Get('metadata').GetTitle()) == report
for stem in ['step1_phase_boost_linear','step1_phase_boost_log',
             'step2_seven_channel_fit','step2_seven_channel_fit_log']:
    assert f.Get('c_'+stem)
    for ext in ['png','pdf']: assert (base/(stem+'.'+ext)).stat().st_size > 1000
f.Close()
print('PASS: ROOT round-trip, all-bin component sums, fit result, metadata, canvases, PNG/PDF outputs.')
