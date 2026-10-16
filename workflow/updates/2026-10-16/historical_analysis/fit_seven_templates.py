"""Bill's two requests: boost weighted phase spaces, then fit seven fixed shapes.

PROVISIONAL: omega/2pi/3pi use Weight, other four channels remain counts.
The fit is nonnegative, fixed-template, observed-variance least squares.
It is not an absolute cross-section measurement or a validated yield extraction.
"""
from pathlib import Path
import argparse
import json
import math
import hashlib
import numpy as np
from scipy.optimize import nnls
import ROOT

BASE = Path(__file__).resolve().parent
CHANNELS = ['omega', 'phi', 'rho', 'eta', 'etap', 'xphsp', 'xphsp3']
TITLES = ['#omega', '#phi', '#rho', '#eta', "#eta'", '2#pi legacy', '3#pi']
COLORS = ['#0000FF', '#E07000', '#00B000', '#FF8000', '#A040FF', '#707070', '#00AFC1']
WEIGHTED = {'omega', 'xphsp', 'xphsp3'}
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)
ROOT.TH1.SetDefaultSumw2(True)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def get(f, key, name):
    obj = f.Get(key)
    assert obj and obj.InheritsFrom('TH1'), key
    h = obj.Clone(name)
    h.SetDirectory(0)
    assert h.GetNbinsX() == 100
    assert h.GetXaxis().GetXmin() == 0 and h.GetXaxis().GetXmax() == 2
    return h


def vals(h):
    return np.array([h.GetBinContent(i) for i in range(1, 101)])


def add(parts, name):
    total = parts[0].Clone(name)
    total.Reset('ICES')
    for h in parts:
        total.Add(h)
    np.testing.assert_allclose(vals(total), np.sum([vals(h) for h in parts], axis=0), rtol=1e-12)
    return total


def draw(parts, total, py, labels, out, name, log=False, fitted=False):
    c = ROOT.TCanvas('c_' + name, '', 1500, 950 if fitted else 900)
    if fitted:
        top = ROOT.TPad(name+'_top', '', 0, .29, 1, 1)
        bottom = ROOT.TPad(name+'_bottom', '', 0, 0, 1, .29)
        top.Draw(); bottom.Draw(); top.cd()
    else:
        top = c
    top.SetLeftMargin(.12); top.SetRightMargin(.035)
    top.SetTopMargin(.045); top.SetBottomMargin(.025 if fitted else .13)
    top.SetTicks(1, 1); top.SetLogy(log)
    frame = total.Clone('frame_' + name); frame.Reset('ICES'); frame.SetTitle('')
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    frame.GetYaxis().SetTitle('Pythia-referenced ordinate [a.u.]' if fitted else 'Diagnostic ordinate [a.u.]')
    frame.GetYaxis().SetTitleOffset(1.3)
    if fitted:
        frame.GetXaxis().SetLabelSize(0); frame.GetXaxis().SetTitleSize(0)
    positives = [v for h in [*parts, total, py] for v in vals(h) if v > 0]
    frame.SetMinimum(max(1e-19, 10**(math.floor(math.log10(min(positives)))-1)) if log else 0)
    frame.SetMaximum(max(total.GetMaximum(), py.GetMaximum()) * (20 if log else 1.18))
    frame.Draw('AXIS')
    leg = ROOT.TLegend(.735, .48 if fitted else .45, .955, .93)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(.024)
    for h, label, color in zip(parts, labels, COLORS):
        h.SetLineColor(ROOT.TColor.GetColor(color)); h.SetLineWidth(2)
        h.Draw('HIST SAME'); leg.AddEntry(h, label, 'l')
    total.SetLineColor(ROOT.kRed); total.SetLineWidth(3); total.Draw('HIST SAME')
    py.SetLineColor(ROOT.kBlack); py.SetLineWidth(3); py.Draw('HIST SAME')
    leg.AddEntry(total, 'Template fit (provisional)' if fitted else 'Mixed sum (diagnostic)', 'l')
    leg.AddEntry(py, 'Pythia Background', 'l'); leg.Draw(); top.RedrawAxis()
    if fitted:
        bottom.cd(); bottom.SetLeftMargin(.12); bottom.SetRightMargin(.035)
        bottom.SetTopMargin(.02); bottom.SetBottomMargin(.34); bottom.SetTicks(1, 1)
        resid = py.Clone('residual_'+name); resid.Add(total, -1)
        # Residual content only: no misleading independent sum of fitted-model errors.
        for i in range(102): resid.SetBinError(i, 0)
        resid.SetLineWidth(2)
        resid.GetYaxis().SetTitle('Pythia - fit [a.u.]')
        resid.GetXaxis().SetTitle('Missing mass [GeV]')
        resid.GetXaxis().SetLabelSize(.095); resid.GetXaxis().SetTitleSize(.11)
        resid.GetYaxis().SetLabelSize(.08); resid.GetYaxis().SetTitleSize(.075)
        resid.GetYaxis().SetTitleOffset(.65); resid.GetYaxis().SetNdivisions(505)
        lim = max(abs(v) for v in vals(resid)) * 1.18
        resid.SetMinimum(-lim); resid.SetMaximum(lim); resid.Draw('HIST')
        zero = ROOT.TLine(0, 0, 2, 0); zero.SetLineStyle(2); zero.Draw()
    c.cd()
    for ext in ['png', 'pdf']: c.SaveAs(str(out/(name+'.'+ext)))
    c.Write(); c.Close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=BASE/'meeting_20260916_template_fit')
    args = ap.parse_args(); out = args.output
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'results.root').exists(), 'Use a new output directory; old results are preserved.'
    prior = BASE/'meeting_20260916_weight_applied'
    meta = json.loads((prior/'summary.json').read_text())
    original = json.loads((BASE/'meeting_20260916_overlay/summary.json').read_text())
    pybase = BASE.parent/'BACKGROUND_ROOT_ANALYSIS/Pythia_Background'
    pymeta = json.loads((pybase/'summary.json').read_text())
    assert meta['cuts'] == pymeta['cut']
    hashes = dict(meta['sources_sha256'])
    for p in [prior/'summary.json', prior/'weight_applied_overlay.root',
              pybase/'histograms.root', Path(__file__), BASE/'fit_seven_templates.C']:
        hashes[str(p)] = sha(p)
    for p, digest in hashes.items(): assert sha(p) == digest, ('Changed input', p)
    f = ROOT.TFile.Open(str(prior/'weight_applied_overlay.root'), 'READ')
    raw = [get(f, ch+'_event_contribution_after_cuts', ch+'_input_after_cuts') for ch in CHANNELS]
    previous = [get(f, ch+'_partial_weighting_reference_after', ch+'_previous') for ch in CHANNELS]
    pydisplay = get(f, 'pythia_background_unchanged', 'pythia_display_unchanged'); f.Close()
    f = ROOT.TFile.Open(str(pybase/'histograms.root'), 'READ')
    pyraw = get(f, 'mx_raw_after', 'pythia_counts_after_cuts'); f.Close()
    display = meta['common_display_scale_retained']
    np.testing.assert_allclose(vals(pydisplay), vals(pyraw)*display, rtol=1e-12)
    assert np.all(vals(pyraw) == np.round(vals(pyraw)))
    for h in raw: assert np.all(vals(h) >= 0) and h.Integral() > 0
    # Independently refill every template from its original SIMC events.
    for ch, h in zip(CHANNELS, raw):
        tf = ROOT.TFile.Open(str(BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root')), 'READ')
        tree = tf.Get('h10'); cut = ROOT.TTreeFormula('cut_'+ch, meta['cuts'], tree)
        check = h.Clone('check_'+ch); check.Reset('ICES'); accepted = 0
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            if cut.EvalInstance():
                check.Fill(float(tree.missmass), float(tree.Weight) if ch in WEIGHTED else 1.)
                accepted += 1
        assert accepted == meta['channel_checks'][ch]['accepted_events']
        for i in range(102):
            assert math.isclose(check.GetBinContent(i), h.GetBinContent(i), rel_tol=1e-10, abs_tol=1e-28)
            assert math.isclose(check.GetBinError(i), h.GetBinError(i), rel_tol=1e-10, abs_tol=1e-28)
        del cut; tf.Close()
    fout = ROOT.TFile(str(out/'results.root'), 'RECREATE')
    # Step 1: exactly x1e10 extra for BOTH phase spaces; everything else untouched.
    boosted = []
    for ch, old in zip(CHANNELS, previous):
        h = old.Clone(ch+'_phase_boost_step1')
        factor = 1e10 if ch in {'xphsp', 'xphsp3'} else 1
        h.Scale(factor)
        np.testing.assert_allclose(vals(h), vals(old)*factor, rtol=1e-12)
        boosted.append(h)
    boostsum = add(boosted, 'phase_boost_sum')
    labels1 = ['#omega: Weight', '#phi / 10: counts', '#rho #times 2: counts',
               '#eta / 22.5: counts', "#eta' / 15: counts",
               '2#pi: Weight #times 4 #times 10^{10}', '3#pi: Weight #times 10^{10}']
    for log in [False, True]:
        draw(boosted, boostsum, pydisplay, labels1, out,
             'step1_phase_boost_'+('log' if log else 'linear'), log=log)
    # Step 2: unit-sum templates are only a numerically stable reparameterization.
    # No old relative scales, omega peak matching, or Step-1 boost enters the fit.
    areas = np.array([h.Integral() for h in raw])
    templates = []
    for ch, h, area in zip(CHANNELS, raw, areas):
        t = h.Clone(ch+'_unit_sum_template'); t.Scale(1/area); templates.append(t)
    matrix = np.array([vals(t) for t in templates]).T
    counts = vals(pyraw); n = float(counts.sum())
    # Include ALL 100 bins in [0,2), including bins not covered by MC templates.
    # Poisson count variance is n_i; empty observed bins use variance=1.
    sigma = np.sqrt(np.maximum(counts, 1.))
    design = matrix*n/sigma[:, None]; target = counts/sigma
    coefficients, norm = nnls(design, target)
    gradient = design.T @ (design@coefficients-target)
    assert all(g >= -1e-4 for g in gradient[coefficients < 1e-10])
    active = coefficients > 1e-10
    # Exact synthetic closure of the nonnegative linear-template solver.
    truth = np.array([.12, .03, .15, .01, .02, .27, .4])
    recovered, _ = nnls(design, design@truth)
    np.testing.assert_allclose(recovered, truth, atol=1e-9)
    ROOT.gROOT.ProcessLine('.L '+str(BASE/'fit_seven_templates.C'))
    for j, t in enumerate(templates): ROOT.set_seven_template(j, t)
    fitdata = pyraw.Clone('fit_data_fraction_per_bin'); fitdata.Scale(1/n)
    for i in range(1, 101): fitdata.SetBinError(i, float(sigma[i-1]/n))
    fun = ROOT.TF1('seven_channel_fit_function', ROOT.seven_template_total, 0., 2., 7)
    for j, ch in enumerate(CHANNELS):
        fun.SetParName(j, ch); fun.SetParameter(j, max(float(coefficients[j]), 1e-6))
        fun.SetParLimits(j, 0., 10.)
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
    result = fitdata.Fit(fun, 'SQNR0')
    root_coefs = np.array([fun.GetParameter(j) for j in range(7)])
    root_chi2 = float(result.Chi2()); chi2 = float(norm**2)
    assert int(result) == 0, ('ROOT minimization failed', int(result))
    assert abs(root_chi2-chi2) < max(.01, chi2*1e-7), (root_chi2, chi2)
    np.testing.assert_allclose(matrix@root_coefs, matrix@coefficients, rtol=1e-3, atol=1e-6)
    # Use exact NNLS optimum (zero at a boundary is exact), cross-checked by ROOT.
    fitted = []
    for ch, h, a in zip(CHANNELS, templates, coefficients):
        part = h.Clone(ch+'_fitted_display'); part.Scale(float(a*n*display)); fitted.append(part)
    fitsum = add(fitted, 'fitted_sum_display')
    labels2 = [title+(': Weight' if ch in WEIGHTED else ': counts') for title, ch in zip(TITLES, CHANNELS)]
    draw(fitted, fitsum, pydisplay, labels2, out, 'step2_seven_channel_fit', fitted=True)
    draw(fitted, fitsum, pydisplay, labels2, out, 'step2_seven_channel_fit_log', log=True, fitted=True)
    residual = pyraw.Clone('residual_counts'); pull = pyraw.Clone('residual_in_Pythia_sigma')
    prediction = matrix@coefficients*n
    for i in range(1, 101):
        residual.SetBinContent(i, float(counts[i-1]-prediction[i-1])); residual.SetBinError(i, float(sigma[i-1]))
        pull.SetBinContent(i, float((counts[i-1]-prediction[i-1])/sigma[i-1])); pull.SetBinError(i, 0.)
    covariance = np.zeros((7, 7))
    covariance[np.ix_(active, active)] = np.linalg.inv(design[:, active].T@design[:, active])
    coef_report = {}
    for j, ch in enumerate(CHANNELS):
        oldfactor = meta['channel_checks'][ch]['channel_scale']*meta['common_SIMC_scale_retained']
        coef_report[ch] = {
            'fill_weight': 'Weight' if ch in WEIGHTED else '1 (unweighted provisional template)',
            'template_integral': float(areas[j]), 'fraction_of_Pythia_count_total': float(coefficients[j]),
            'fitted_count_equivalent': float(coefficients[j]*n),
            'scale_per_original_histogram_to_Pythia_counts': float(coefficients[j]*n/areas[j]),
            'multiplier_relative_to_previous_display_curve': float(coefficients[j]*n/areas[j]/oldfactor),
            'at_nonnegative_boundary': bool(not active[j]),
            'conditional_count_error_fixed_templates': float(math.sqrt(covariance[j,j])*n) if active[j] else None,
            'root_fraction': float(root_coefs[j])}
    uncovered = (matrix.sum(axis=1) == 0) & (counts > 0)
    report = {
        'cuts': meta['cuts'], 'root_version': ROOT.gROOT.GetVersion(),
        'step1': {'additional_phase_factors': {'xphsp': 1e10, 'xphsp3': 1e10},
                  'other_channels_and_Pythia_unchanged': True, 'old_scales_retained': meta['channel_checks']},
        'step2': {'model': 'sum of seven fixed histogram shapes with nonnegative free amplitudes',
                  'source_example': 'https://github.com/billlee77/omega_analysis/blob/master/fitting_functions.C',
                  'method': 'observed-count-variance least squares; scipy NNLS independently checked by ROOT Minuit2',
                  'fit_range_GeV': [0, 2], 'bins': 100, 'bin_width_GeV': .02,
                  'variance': 'Pythia raw count per bin; 1 for empty Pythia bins',
                  'old_scales_used_in_fit': False, 'shape_smoothing_or_interpolation': False,
                  'chi2_Pythia_variance_only': chi2, 'ndf_nominal_100_minus_7': 93,
                  'chi2_over_nominal_ndf': chi2/93,
                  'active_parameters': int(active.sum()), 'ROOT_status': int(result),
                  'ROOT_covariance_status': int(result.CovMatrixStatus()), 'ROOT_chi2': root_chi2,
                  'matrix_rank': int(np.linalg.matrix_rank(design)),
                  'matrix_condition_number': float(np.linalg.cond(design)),
                  'Pythia_in_range_counts': n, 'Pythia_underflow': pyraw.GetBinContent(0),
                  'Pythia_overflow': pyraw.GetBinContent(101),
                  'uncovered_bins': [{'bin': int(i+1), 'center_GeV': pyraw.GetBinCenter(int(i+1)),
                                      'Pythia_count': float(counts[i])} for i in np.where(uncovered)[0]],
                  'uncovered_counts': float(counts[uncovered].sum()),
                  'coefficients': coef_report,
                  'conditional_fraction_covariance': covariance.tolist()},
        'caveats': [
            'Only omega, legacy 2pi and 3pi are Weight-filled. phi/rho/eta/etap are count templates pending weight correction.',
            'Pythia is supplied Background tree after the same cuts, with no Weight branch.',
            meta['pythia_variable_caveat'],
            '1000 stored SIMC events per channel: fixed-template fit ignores finite MC shape uncertainty. Conditional errors are not total uncertainties.',
            'A converged optimizer does not mean the template model describes Pythia. Inspect residuals and chi2.',
            'Amplitudes are diagnostic count equivalents, not physical yields, branching fractions or cross sections.',
            'Flows excluded from fit; all 100 visible bins retained, including uncovered bins. No hidden smoothing or added background.',
            'Parameters at zero boundary have no symmetric error reported; Gaussian covariance is conditional on active set.',
            'Old absolute display multiplier retained only to show the same black curve; not used to determine fit amplitudes.'
        ],
        'sources_sha256': hashes,
        'verification': ['All seven templates independently refilled from source events with cuts and appropriate fill weight.',
                         'Step1: only phase spaces multiplied by 1e10; black curve unchanged.',
                         'NNLS synthetic seven-component closure passed; independent ROOT minimization agrees.',
                         'Histogram sums checked bin by bin; inputs hash-verified before/after; previous outputs preserved.']}
    for h in [*raw, *templates, *boosted, boostsum, *fitted, fitsum, pyraw, pydisplay, residual, pull, fitdata]: h.Write()
    result.Get().Write('ROOT_fit_result'); fun.Write()
    ROOT.TNamed('metadata', json.dumps(report)).Write(); fout.Close()
    for p, digest in hashes.items(): assert sha(p) == digest, ('Changed input', p)
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    rows = ['channel\tfill\tfitted_count_equivalent\tscale_to_original_histogram\tmultiplier_vs_previous_curve\tboundary']
    for ch, d in coef_report.items():
        rows.append(f"{ch}\t{d['fill_weight']}\t{d['fitted_count_equivalent']:.10g}\t{d['scale_per_original_histogram_to_Pythia_counts']:.10g}\t{d['multiplier_relative_to_previous_display_curve']:.10g}\t{d['at_nonnegative_boundary']}")
    (out/'fit_coefficients.txt').write_text('\n'.join(rows)+'\n')
    print(json.dumps({'output':str(out), 'fit':report['step2']}, indent=2))


if __name__ == '__main__':
    main()
