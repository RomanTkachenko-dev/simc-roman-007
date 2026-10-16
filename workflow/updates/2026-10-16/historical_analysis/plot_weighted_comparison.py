"""Signed event Weight after acceptance cuts; retain requested display scales."""
from pathlib import Path
import json
import math
import ROOT
from compare_pythia_scaled_sum import CHANNELS, sha
from normalize_omega_peak import BASE, peak
from plot_scaled_missing_mass import COLORS

OUT = BASE/'pythia_vs_simc_weighted'
SCALES = dict(zip(CHANNELS, [1., .1, 2., 1./22.5, 1./15., 4.]))
LABELS = ['#omega', '#phi / 10', '#rho #times 2', '#eta / 22.5', "#eta' / 15", 'X-phase-space #times 4']


def check(a, b):
    assert math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-28), (a, b)


def draw(total, components, name, ytitle, pythia=None):
    c = ROOT.TCanvas('canvas_'+name, '', 1450, 880)
    c.SetLeftMargin(.12); c.SetRightMargin(.04); c.SetBottomMargin(.12); c.SetTopMargin(.045); c.SetTicks(1, 1)
    frame = total.Clone('frame_'+name); frame.Reset('ICES'); frame.SetTitle('')
    frame.GetXaxis().SetTitle('Missing mass [GeV]'); frame.GetYaxis().SetTitle(ytitle)
    for axis in [frame.GetXaxis(), frame.GetYaxis()]:axis.SetTitleSize(.04)
    frame.GetYaxis().SetTitleOffset(1.4)
    histograms = [total, *components]+([pythia] if pythia else [])
    ymin = min(0., *(h.GetMinimum() for h in histograms))
    ymax = max(0., *(h.GetMaximum() for h in histograms))
    span = ymax-ymin
    frame.SetMinimum(ymin-.06*span if ymin < 0 else 0.)
    frame.SetMaximum(ymax+.16*span); frame.Draw('AXIS')
    zero = ROOT.TLine(0., 0., 2., 0.); zero.SetLineColor(ROOT.kGray+1); zero.SetLineStyle(3); zero.Draw()
    for h in components:h.Draw('HIST SAME')
    total.Draw('HIST SAME')
    if pythia:pythia.SetLineColor(ROOT.kBlack); pythia.SetLineWidth(3); pythia.Draw('HIST SAME')
    legend = ROOT.TLegend(.73, .48, .947, .93)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(.026)
    for h, label in zip(components, LABELS):legend.AddEntry(h, label, 'l')
    legend.AddEntry(total, 'Weighted SIMC sum', 'l')
    if pythia:legend.AddEntry(pythia, 'Pythia Background', 'l')
    legend.Draw(); c.RedrawAxis()
    for ext in ['png', 'pdf']:c.SaveAs(str(OUT/(name+'.'+ext)))
    c.Write(); c.Close()


def main():
    OUT.mkdir(exist_ok=True)
    pm = json.loads((BASE.parent/'BACKGROUND_ROOT_ANALYSIS/Pythia_Background/summary.json').read_text())
    cuts = pm['cut']
    report = {'cuts': cuts, 'selection_expression': 'Weight*('+cuts+')',
              'channel_scales': SCALES, 'channels': {},
              'weight_handling': 'Signed Weight for every accepted SIMC event; no abs, clipping, or removal of negative weights.',
              'pythia': {'original_file': pm['file'], 'tree': pm['tree'], 'histogram': 'mx_raw_after', 'event_weight': 'Not available; raw counts retained.'},
              'warning': 'Diagnostic signed-weight comparison, not validated physical yields. Negative weights need interpretation/validation. Channel factors are user-selected, not luminosity factors.',
              'variable_assumption': 'Pythia mx is provisionally treated as reconstructed missing mass in GeV; physical equivalence to SIMC missmass is unverified.',
              'error_handling': 'Sumw2; independent-channel variances added. Peak-normalization-factor uncertainty not propagated.'}
    source_hashes = {}; originals = []; components = []
    for ch, color in zip(CHANNELS, COLORS):
        path = BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root')
        source_hashes[str(path)] = sha(path)
        f = ROOT.TFile.Open(str(path), 'READ'); assert f and not f.IsZombie()
        tree = f.Get('h10'); assert tree and tree.GetBranch('Weight') and tree.GetBranch('missmass')
        selection = ROOT.TTreeFormula('cut_'+ch, cuts, tree)
        h = ROOT.TH1D(ch+'_signed_weight_after_cuts', '', 100, 0, 2); h.Sumw2()
        weights = []; negative_all = 0
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            w, m = float(tree.Weight), float(tree.missmass)
            assert math.isfinite(w) and math.isfinite(m)
            negative_all += w < 0
            if selection.EvalInstance():
                h.Fill(m, w); weights.append(w)
        assert len(weights) == tree.GetEntries(cuts)
        check(h.Integral(0,101), math.fsum(weights))
        check(sum(h.GetBinError(i)**2 for i in range(102)), math.fsum(w*w for w in weights))
        saved_path = path.parent/'plots/histograms.root'
        saved_meta = json.loads((saved_path.parent/'summary.json').read_text())
        assert saved_meta['cut'] == cuts
        saved = ROOT.TFile.Open(str(saved_path), 'READ')
        previous = saved.Get('missmass_weighted_after')
        for i in range(102):
            check(h.GetBinContent(i), previous.GetBinContent(i))
            check(h.GetBinError(i), previous.GetBinError(i))
        saved.Close()
        report['channels'][ch] = {'file': str(path), 'all_events': int(tree.GetEntries()),
            'accepted': len(weights), 'negative_all': negative_all,
            'negative_accepted': sum(w<0 for w in weights), 'positive_accepted': sum(w>0 for w in weights),
            'zero_accepted': sum(w==0 for w in weights), 'sum_weight_after_cuts': math.fsum(weights),
            'sum_weight_visible': h.Integral(), 'scale': SCALES[ch]}
        del selection
        f.Close()
        hc = h.Clone(ch+'_weighted_scaled'); hc.Scale(SCALES[ch])
        hc.SetLineColor(color); hc.SetLineWidth(2); hc.SetLineStyle(1)
        originals.append(h); components.append(hc)
    total = components[0].Clone('weighted_simc_sum'); total.Reset('ICES')
    for h in components:total.Add(h)
    total.SetLineColor(ROOT.TColor.GetColor('#FF0000')); total.SetLineWidth(3)
    for i in range(102):
        check(total.GetBinContent(i), math.fsum(h.GetBinContent(i) for h in components))
        check(total.GetBinError(i)**2, math.fsum(h.GetBinError(i)**2 for h in components))
    pp = BASE.parent/'BACKGROUND_ROOT_ANALYSIS/Pythia_Background/histograms.root'
    source_hashes[str(pp)] = sha(pp)
    pf = ROOT.TFile.Open(str(pp), 'READ'); pythia = pf.Get('mx_raw_after').Clone('pythia_raw_after'); pf.Close()
    assert pythia.GetNbinsX() == 100
    for i in range(1,102):check(pythia.GetBinLowEdge(i), total.GetBinLowEdge(i))
    sp, ppeak = peak(total), peak(pythia)
    report['weighted_sum_omega_region_peak'] = sp; report['pythia_omega_region_peak'] = ppeak
    output = ROOT.TFile(str(OUT/'weighted_comparison.root'), 'RECREATE')
    for h in [*originals, *components, total, pythia]:h.Write()
    draw(total, components, 'signed_weighted_simc', 'Scaled sum of Weight / 0.02 GeV')
    if sp['height'] > 0:
        factor = ppeak['height']/sp['height']
        display = 1./pythia.Integral('width')
        dp = pythia.Clone('pythia_reference_density'); dp.Scale(display)
        dt = total.Clone('weighted_simc_peak_matched_density'); dt.Scale(factor*display)
        dc = [h.Clone(ch+'_weighted_peak_matched_density') for ch,h in zip(CHANNELS,components)]
        for original, h in zip(components,dc):
            h.Scale(factor*display)
            for i in range(102):
                check(h.GetBinContent(i), original.GetBinContent(i)*factor*display)
                check(h.GetBinError(i), original.GetBinError(i)*factor*display)
        check(peak(dt)['height'], peak(dp)['height'])
        for i in range(102):check(dt.GetBinContent(i), math.fsum(h.GetBinContent(i) for h in dc))
        report['normalization'] = {'method': 'Same previous prescription: maxima of TOTAL spectra in [0.70,0.86) GeV; recomputed for signed weights, not reused from raw counts.',
            'common_SIMC_factor': factor, 'common_display_factor': display,
            'SIMC_density_area': dt.Integral('width'), 'Pythia_density_area': dp.Integral('width')}
        for h in [dp,dt,*dc]:h.Write()
        draw(dt,dc,'weighted_omega_matched_channels','Pythia-referenced density [GeV^{-1}]',dp)
    else:
        report['normalization'] = 'Not applied: no positive weighted total peak in omega window; a sign-flipping normalization would be misleading.'
    report['source_sha256'] = source_hashes
    report['verification'] = 'Rebuilt from event trees; cut counts, signed sums, Sumw2 and every bin agree with saved weighted histograms; channel sums checked; source hashes unchanged.'
    assert all(sha(Path(p)) == digest for p,digest in source_hashes.items())
    ROOT.TNamed('metadata',json.dumps(report)).Write(); output.Close()
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':main()
