"""Manual phase-space rescaling, linear ordinate, no fit or new normalization."""
from pathlib import Path
import json
import math
import ROOT
import apply_requested_weights as plotting

BASE = Path(__file__).resolve().parent
PRIOR = BASE/'meeting_20260916_weight_applied'
OUT = BASE/'manual_phase_scaling_3e11_7e12'
FACTORS = {'xphsp': 3e11, 'xphsp3': 7e12}


def main():
    OUT.mkdir(exist_ok=True)
    target = OUT/'manual_phase_scales.root'
    assert not target.exists(), 'Previous results will not be overwritten.'
    source_path = PRIOR/'weight_applied_overlay.root'
    meta = json.loads((PRIOR/'summary.json').read_text())
    sources = {str(p): plotting.sha(p) for p in [source_path, PRIOR/'summary.json', Path(__file__), BASE/'apply_requested_weights.py']}
    source = ROOT.TFile.Open(str(source_path), 'READ')
    assert source and not source.IsZombie()
    parts = []
    for ch in plotting.CHANNELS:
        old = plotting.get(source, ch+'_partial_weighting_reference_after', ch+'_reference')
        h = old.Clone(ch+'_manual_scaled')
        factor = FACTORS.get(ch, 1.)
        h.Scale(factor)
        for i in range(102):
            plotting.close(h.GetBinContent(i), old.GetBinContent(i)*factor)
            plotting.close(h.GetBinError(i), old.GetBinError(i)*factor)
        h.SetLineColor(ROOT.TColor.GetColor(plotting.COLORS[ch]))
        parts.append(h)
    py = plotting.get(source, 'pythia_background_unchanged', 'pythia_background_unchanged')
    plotting.equal(py, source.Get('pythia_background_unchanged'))
    source.Close()
    total = parts[0].Clone('manual_diagnostic_sum'); total.Reset('ICES')
    for h in parts: total.Add(h)
    for i in range(102):
        plotting.close(total.GetBinContent(i), math.fsum(h.GetBinContent(i) for h in parts))
        plotting.close(total.GetBinError(i)**2, math.fsum(h.GetBinError(i)**2 for h in parts))
    plotting.LABELS = ['#omega: Weight', '#phi / 10: counts', '#rho #times 2: counts',
        '#eta / 22.5: counts', "#eta' / 15: counts",
        '2#pi: Weight #times 4 #times 3#upoint10^{11}', '3#pi: Weight #times 7#upoint10^{12}']
    report = {
        'fit_performed': False, 'logarithmic_axis': False, 'new_normalization': False,
        'additional_multipliers_relative_to_previous_weight_applied_plot': FACTORS,
        'meaning': 'Replace the trial 1e10 with 3e11 (2pi) and 7e12 (3pi); do not multiply on top of 1e10.',
        'previous_2pi_times_4_retained': True,
        'cuts': meta['cuts'], 'weighted_channels': meta['weighted_channels'],
        'all_other_channels_and_Pythia_unchanged': True,
        'warning': meta['warning'], 'variable_caveat': meta['pythia_variable_caveat'],
        'sources_sha256': sources,
        'verification': 'Contents and errors scaled bin-by-bin; sum and variance checked including flows; ROOT reopened; sources unchanged.'}
    dest = ROOT.TFile(str(target), 'RECREATE')
    for h in [*parts, total, py]: h.Write()
    plotting.draw(parts, total, py, OUT, 'manual_phase_scales_linear', False)
    ROOT.TNamed('metadata', json.dumps(report)).Write(); dest.Close()
    check = ROOT.TFile.Open(str(target), 'READ')
    assert check and not check.IsZombie() and not check.TestBit(ROOT.TFile.kRecovered)
    for h in [*parts, total, py]: plotting.equal(h, check.Get(h.GetName()))
    assert json.loads(check.Get('metadata').GetTitle()) == report
    check.Close()
    for path, digest in sources.items(): assert plotting.sha(Path(path)) == digest
    (OUT/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS: manual linear plot; no fit; unchanged other channels/Pythia; scale/error/sum/source/ROOT checks.')
    print(OUT/'manual_phase_scales_linear.png')


if __name__ == '__main__': main()
