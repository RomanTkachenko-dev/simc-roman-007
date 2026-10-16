"""Match the total SIMC omega-region peak height to Pythia Background."""
from pathlib import Path
import json
import math
import ROOT
from compare_pythia_scaled_sum import CHANNELS, render, sha

BASE = Path(__file__).resolve().parent
OUT = BASE / 'pythia_vs_simc_omega_peak'
WINDOW = (0.70, 0.86)


def peak(hist):
    bins = [i for i in range(1, hist.GetNbinsX()+1)
            if WINDOW[0] <= hist.GetBinCenter(i) < WINDOW[1]]
    i = max(bins, key=hist.GetBinContent)
    return {'bin': i, 'low_GeV': hist.GetBinLowEdge(i),
            'high_GeV': hist.GetBinLowEdge(i+1), 'height': hist.GetBinContent(i)}


def close(a, b):
    assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-11), (a, b)


def main():
    OUT.mkdir(exist_ok=True)
    sp = BASE / 'scaled_missing_mass/scaled_missing_mass.root'
    pp = BASE.parent / 'BACKGROUND_ROOT_ANALYSIS/Pythia_Background/histograms.root'
    sm = json.loads((sp.parent/'summary.json').read_text())
    pm = json.loads((pp.parent/'summary.json').read_text())
    assert sm['cuts'] == pm['cut'], 'Acceptance cuts differ'
    hashes = {str(p): sha(p) for p in (sp, pp)}
    sf = ROOT.TFile.Open(str(sp), 'READ')
    pf = ROOT.TFile.Open(str(pp), 'READ')
    assert sf and not sf.IsZombie() and pf and not pf.IsZombie()
    total = sf.Get('sum_scaled_missing_mass').Clone('simc_sum_before_peak_matching')
    pythia = pf.Get('mx_raw_after').Clone('pythia_background_counts')
    components = [sf.Get(ch+'_scaled_missing_mass').Clone(ch+'_before_peak_matching') for ch in CHANNELS]
    sf.Close(); pf.Close()
    for h in [total, pythia, *components]:
        assert h.GetNbinsX() == 100
        for i in range(1, 102):close(h.GetBinLowEdge(i), (i-1)*.02)
        for i in range(102):assert math.isfinite(h.GetBinContent(i)) and h.GetBinContent(i) >= 0
    for i in range(102):close(total.GetBinContent(i), sum(h.GetBinContent(i) for h in components))
    simc_peak, pythia_peak = peak(total), peak(pythia)
    assert simc_peak['height'] > 0 and pythia_peak['height'] > 0
    factor = pythia_peak['height']/simc_peak['height']
    matched = total.Clone('simc_sum_peak_matched_counts'); matched.Scale(factor)
    mc = [h.Clone(ch+'_peak_matched_counts') for ch, h in zip(CHANNELS, components)]
    for h in mc:h.Scale(factor)
    close(peak(matched)['height'], pythia_peak['height'])

    # A single DISPLAY factor is applied to BOTH spectra and all components.
    # Pythia has area one; the matched SIMC area is deliberately not forced to one.
    display_factor = 1./pythia.Integral('width')
    dp = pythia.Clone('pythia_reference_density'); dp.Scale(display_factor)
    dt = matched.Clone('simc_sum_peak_matched_density'); dt.Scale(display_factor)
    dc = [h.Clone(ch+'_peak_matched_density') for ch, h in zip(CHANNELS, mc)]
    for h in dc:h.Scale(display_factor)
    for i in range(102):
        close(matched.GetBinContent(i), sum(h.GetBinContent(i) for h in mc))
        close(dt.GetBinContent(i), sum(h.GetBinContent(i) for h in dc))
        for original, raw, density in zip([total, *components], [matched, *mc], [dt, *dc]):
            close(raw.GetBinContent(i), original.GetBinContent(i)*factor)
            close(raw.GetBinError(i), original.GetBinError(i)*factor)
            close(density.GetBinContent(i), original.GetBinContent(i)*factor*display_factor)
            close(density.GetBinError(i), original.GetBinError(i)*factor*display_factor)
    close(peak(dt)['height'], peak(dp)['height'])
    close(dp.Integral('width'), 1.)
    labels = ['#omega', '#phi / 10', '#rho #times 2', '#eta / 15', "#eta' / 15", 'X-phase-space #times 4']
    assert sm['factors'] == dict(zip(CHANNELS, [1., .1, 2., 1./15., 1./15., 4.]))
    report = {
        'method': 'Maximum bin height of the TOTAL SIMC spectrum in the omega window matched to the maximum Pythia Background bin height in the same window. No fit or background subtraction.',
        'window_GeV': WINDOW, 'bin_width_GeV': .02,
        'simc_peak_before': simc_peak, 'pythia_peak': pythia_peak,
        'common_SIMC_factor_relative_to_scaled_counts': factor,
        'SIMC_factor_relative_to_previous_unit_area_density': factor*display_factor*total.Integral('width'),
        'unchanged_relative_channel_factors': sm['factors'],
        'common_display_factor_for_both_spectra': display_factor,
        'density_peak_height_SIMC': peak(dt)['height'],
        'density_peak_height_Pythia': peak(dp)['height'],
        'density_area_SIMC': dt.Integral('width'),
        'density_area_Pythia': dp.Integral('width'),
        'cuts': sm['cuts'], 'sources_sha256': hashes,
        'pythia_original_file': pm['file'], 'pythia_tree': pm['tree'],
        'pythia_histogram': 'mx_raw_after', 'simc_original_histogram': 'missmass_raw_after',
        'variable_assumption': 'Pythia mx is provisionally treated as reconstructed missing mass in GeV. Its writer/definition has not been supplied or independently verified.',
        'cautions': ['Empirical peak-height matching, not luminosity or cross-section normalization.',
                     'The matched SIMC peak includes all six channels, not the isolated omega component.',
                     'Raw-count histograms with user channel scales; event Weight is not used.',
                     'Errors scaled conditional on the chosen factor; uncertainty of the estimated normalization factor is not propagated.',
                     'Underflow/overflow excluded from peak search and display denominator; their stored contents are scaled consistently.',
                     'Pythia is a separate comparison spectrum and is not added to the SIMC sum.'],
        'verification': 'Identical binning/cuts; equal omega-region maximum heights; all component bins and errors share one factor; component sums agree; input hashes unchanged.'
    }
    output = ROOT.TFile(str(OUT/'omega_peak_matched.root'), 'RECREATE')
    for h in [total, pythia, *components, matched, *mc, dp, dt, *dc]:h.Write()
    render(dp, dt, dc, labels, OUT, 'omega_peak_matched_channels', True, True,
           y_title='Pythia-referenced density [GeV^{-1}]')
    render(dp, dt, dc, labels, OUT, 'omega_peak_matched_totals', True,
           y_title='Pythia-referenced density [GeV^{-1}]')
    render(pythia, matched, mc, labels, OUT, 'omega_peak_matched_counts', False, True)
    ROOT.TNamed('analysis_metadata', json.dumps(report)).Write()
    output.Close()
    assert all(sha(Path(p)) == value for p, value in hashes.items()), 'Input changed'
    (OUT/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':main()
