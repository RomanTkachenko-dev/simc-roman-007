"""Reduce eta alone by another 1.5, retaining the established omega scale."""
import json
import ROOT
from normalize_omega_peak import BASE, peak, close
from compare_pythia_scaled_sum import CHANNELS, render, sha

SOURCE = BASE / 'pythia_vs_simc_omega_peak'
OUT = BASE / 'pythia_vs_simc_omega_peak_eta_reduced'
COLOR = '#FF0000'


def main():
    OUT.mkdir(exist_ok=True)
    source = SOURCE/'omega_peak_matched.root'
    digest = sha(source)
    report = json.loads((SOURCE/'summary.json').read_text())
    f = ROOT.TFile.Open(str(source), 'READ')
    assert f and not f.IsZombie()
    objects = []
    groups = []
    for suffix, pythia_name in [('counts', 'pythia_background_counts'),
                                ('density', 'pythia_reference_density')]:
        pythia = f.Get(pythia_name).Clone(pythia_name)
        old_sum = f.Get('simc_sum_peak_matched_'+suffix)
        total = old_sum.Clone('simc_sum_peak_matched_'+suffix)
        total.Reset('ICES')
        components = []
        for channel in CHANNELS:
            name = channel+'_peak_matched_'+suffix
            # The prior file also has a canvas named omega_peak_matched_counts;
            # cycle 1 is the histogram, the later cycle is the canvas.
            original = f.Get(name+';1')
            assert original.InheritsFrom('TH1'), name
            h = original.Clone(name)
            scale = 1./1.5 if channel == 'eta' else 1.
            h.Scale(scale)
            for i in range(102):
                close(h.GetBinContent(i), original.GetBinContent(i)*scale)
                close(h.GetBinError(i), original.GetBinError(i)*scale)
            components.append(h)
            total.Add(h)
        total.SetLineColor(ROOT.TColor.GetColor(COLOR))
        for i in range(102):
            close(total.GetBinContent(i), sum(h.GetBinContent(i) for h in components))
            close(total.GetBinError(i)**2, sum(h.GetBinError(i)**2 for h in components))
            close(total.GetBinContent(i), old_sum.GetBinContent(i)-f.Get('eta_peak_matched_'+suffix).GetBinContent(i)/3.)
        report['updated_'+suffix] = {'simc_peak': peak(total), 'pythia_peak': peak(pythia),
                                    'simc_integral_width': total.Integral('width')}
        groups.append((pythia, total, components))
        objects.extend([pythia, total, *components])
    f.Close()
    report['revision'] = 'Eta alone reduced by another factor 1.5; sum recomputed and drawn bright red. Previous common omega normalization and all other channels retained exactly.'
    report['prior_relative_channel_factors'] = report.pop('unchanged_relative_channel_factors')
    report['relative_channel_factors'] = dict(report['prior_relative_channel_factors'])
    report['relative_channel_factors']['eta'] /= 1.5
    report['sum_color'] = COLOR
    report['revision_source_sha256'] = digest
    report['density_area_SIMC'] = groups[1][1].Integral('width')
    report['density_peak_height_SIMC'] = peak(groups[1][1])['height']
    report['verification'] = 'Every eta bin/error divided by 1.5; other components unchanged; sum contents and variances checked in all bins; previous common normalization retained; source ROOT unchanged.'
    # Historical unit-area conversion is no longer a uniform factor after this eta-only edit.
    report['prior_SIMC_factor_relative_to_previous_unit_area_density'] = report.pop('SIMC_factor_relative_to_previous_unit_area_density')
    labels = ['#omega', '#phi / 10', '#rho #times 2', '#eta / 22.5', "#eta' / 15", 'X-phase-space #times 4']
    output = ROOT.TFile(str(OUT/'omega_peak_matched.root'), 'RECREATE')
    for h in objects:h.Write()
    pythia, total, components = groups[1]
    for detailed, name in [(True, 'omega_peak_matched_channels'), (False, 'omega_peak_matched_totals')]:
        render(pythia, total, components, labels, OUT, name, True, detailed,
               y_title='Pythia-referenced density [GeV^{-1}]', sum_color=COLOR)
    render(*groups[0], labels, OUT, 'omega_peak_matched_count_comparison', False, True, sum_color=COLOR)
    ROOT.TNamed('analysis_metadata', json.dumps(report)).Write()
    output.Close()
    assert sha(source) == digest
    (OUT/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['relative_channel_factors', 'updated_density', 'verification']}, indent=2))


if __name__ == '__main__':main()
