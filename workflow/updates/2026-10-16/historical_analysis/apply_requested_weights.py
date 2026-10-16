"""Apply original event Weight to omega and the legacy/3pi phase-space curves.

This is deliberately a PARTIALLY weighted diagnostic, not a physical total.
All other curves and all historical scale factors remain unchanged.
"""
from pathlib import Path
import json, math, argparse
import ROOT
from compare_pythia_scaled_sum import sha

BASE=Path(__file__).resolve().parent
PRIOR=BASE/'meeting_20260916_overlay'
CHANNELS=['omega','phi','rho','eta','etap','xphsp','xphsp3']
WEIGHTED={'omega','xphsp','xphsp3'}
COLORS={'omega':'#0000FF','phi':'#E07000','rho':'#00B000','eta':'#FF8000',
        'etap':'#A040FF','xphsp':'#707070','xphsp3':'#00AFC1'}
LABELS=['#omega: Weight','#phi / 10: counts','#rho #times 2: counts',
        '#eta / 22.5: counts',"#eta' / 15: counts",'2#pi legacy #times 4: Weight','3#pi #times 1: Weight']
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)
ROOT.TH1.SetDefaultSumw2(True)

def close(a,b):
    assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-28),(a,b)

def equal(a,b):
    assert a.GetNbinsX()==b.GetNbinsX()
    for i in range(a.GetNbinsX()+2):
        close(a.GetBinContent(i),b.GetBinContent(i))
        close(a.GetBinError(i)**2,b.GetBinError(i)**2)

def get(f,name,clone):
    h=f.Get(name);assert h and h.InheritsFrom('TH1'),name
    h=h.Clone(clone);h.SetDirectory(0);return h

def draw(parts,total,py,out,name,log):
    c=ROOT.TCanvas('canvas_'+name,'',1500,920)
    c.SetLeftMargin(.13);c.SetRightMargin(.035);c.SetBottomMargin(.13);c.SetTopMargin(.045)
    c.SetTicks(1,1);c.SetLogy(log)
    frame=total.Clone('frame_'+name);frame.Reset('ICES');frame.SetTitle('')
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    frame.GetYaxis().SetTitle('Diagnostic ordinate [a.u.]')
    frame.GetYaxis().SetTitleOffset(1.45)
    values=[h.GetBinContent(i) for h in [*parts,total,py] for i in range(1,101) if h.GetBinContent(i)>0]
    frame.SetMinimum(10.**(math.floor(math.log10(min(values)))-1) if log else 0.)
    frame.SetMaximum(max(total.GetMaximum(),py.GetMaximum())*(20 if log else 1.18))
    frame.Draw('AXIS')
    leg=ROOT.TLegend(.72,.46,.955,.89);leg.SetBorderSize(0);leg.SetFillStyle(0)
    leg.SetTextFont(42);leg.SetTextSize(.022)
    for h,label in zip(parts,LABELS):
        h.SetLineWidth(2);h.Draw('HIST SAME');leg.AddEntry(h,label,'l')
    total.SetLineColor(ROOT.kRed);total.SetLineWidth(3);total.Draw('HIST SAME')
    py.SetLineColor(ROOT.kBlack);py.SetLineWidth(3);py.Draw('HIST SAME')
    leg.AddEntry(total,'Mixed sum (diagnostic)','l');leg.AddEntry(py,'Pythia Background','l');leg.Draw()
    c.RedrawAxis()
    for ext in ['png','pdf']:c.SaveAs(str(out/(name+'.'+ext)))
    c.Write();c.Close()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=BASE/'meeting_20260916_weight_applied')
    ap.add_argument('--refresh-own-output',action='store_true',help='Regenerate only this script\'s already verified output directory.')
    args=ap.parse_args();out=args.output
    out.mkdir(parents=True,exist_ok=True)
    if (out/'weight_applied_overlay.root').exists():
        assert args.refresh_own_output,'Refusing to overwrite existing result'
        existing=json.loads((out/'summary.json').read_text())
        assert existing['weighted_channels']==sorted(WEIGHTED)
        assert out.resolve()==(BASE/'meeting_20260916_weight_applied').resolve()
    meta=json.loads((PRIOR/'summary.json').read_text())
    paths=[PRIOR/'summary.json',PRIOR/'overlay_results.root']
    hashes={str(p):sha(p) for p in paths}
    for p,digest in meta['source_hashes'].items():
        assert sha(Path(p))==digest,('Changed input',p)
        hashes[p]=digest
    source=ROOT.TFile.Open(str(PRIOR/'overlay_results.root'),'READ');assert source and not source.IsZombie()
    cut=meta['cuts'];scales=meta['relative_scales_retained']
    common=meta['previous_scale_overlay']['common_SIMC_factor_frozen']
    display=meta['previous_scale_overlay']['common_display_factor_frozen']
    parts=[];unscaled=[];checks={}
    for ch in CHANNELS:
        old=get(source,ch+'_previous_reference_after','old_'+ch)
        h=get(source,ch+('_weighted_after' if ch in WEIGHTED else '_raw_after'),ch+'_event_contribution_after_cuts')
        if ch in WEIGHTED:
            path=BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root')
            tf=ROOT.TFile.Open(str(path),'READ');tree=tf.Get('h10')
            formula=ROOT.TTreeFormula('cuts_'+ch,cut,tree)
            independent=ROOT.TH1D(ch+'_independent_weight_check','',100,0,2)
            values=[]
            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                if formula.EvalInstance():
                    w=float(tree.Weight);assert math.isfinite(w) and w>0
                    independent.Fill(float(tree.missmass),w);values.append(w)
            equal(independent,h)
            assert len(values)==meta['channels'][ch]['after_cuts']
            close(h.Integral(0,101),math.fsum(values))
            close(sum(h.GetBinError(i)**2 for i in range(102)),math.fsum(w*w for w in values))
            del formula;tf.Close()
        unscaled.append(h)
        scaled=h.Clone(ch+'_partial_weighting_reference_after')
        scaled.Scale(scales[ch]*common*display)
        scaled.SetLineColor(ROOT.TColor.GetColor(COLORS[ch]))
        for i in range(102):
            close(scaled.GetBinContent(i),h.GetBinContent(i)*scales[ch]*common*display)
            close(scaled.GetBinError(i),h.GetBinError(i)*scales[ch]*common*display)
        if ch not in WEIGHTED:equal(scaled,old)
        checks[ch]={'fill_weight':'original Weight' if ch in WEIGHTED else '1 (unchanged)',
            'accepted_events':meta['channels'][ch]['after_cuts'],
            'channel_scale':scales[ch],'unscaled_integral_all_bins':h.Integral(0,101),
            'visible_integral_prior':old.Integral(),'visible_integral_now':scaled.Integral(),
            'visible_integral_ratio_to_prior':scaled.Integral()/old.Integral()}
        parts.append(scaled)
    py=get(source,'pythia_previous_reference_after','pythia_background_unchanged')
    oldtotal=get(source,'sum_seven_previous_reference_after','previous_unweighted_diagnostic_sum')
    total=parts[0].Clone('partially_weighted_diagnostic_sum');total.Reset('ICES')
    for h in parts:total.Add(h)
    for i in range(102):
        close(total.GetBinContent(i),math.fsum(h.GetBinContent(i) for h in parts))
        close(total.GetBinError(i)**2,math.fsum(h.GetBinError(i)**2 for h in parts))
    source.Close()
    report={'cuts':cut,'weighted_channels':sorted(WEIGHTED),'channel_checks':checks,
        'common_SIMC_scale_retained':common,'common_display_scale_retained':display,
        'pythia_unchanged':True,'new_fit_or_normalization':False,
        'warning':'PARTIALLY weighted numerical diagnostic only. Four channels remain counts; three use original Weight. Mixed sum has no physical cross-section/yield interpretation. Small weighted magnitudes are retained literally, not rescaled to the old counts or omega peak.',
        'phase_space_note':'xphsp is the existing legacy mode 1, referred to as two-pion phase space by Bill; not a new full exclusive two-pion generator.',
        'pythia_variable_caveat':meta['variable_caveat'],
        'sources_sha256':hashes,
        'verification':'Requested three weighted histograms independently refilled from events, including under/overflow and Sumw2; other four curves exactly unchanged; old scale factors retained; summed bins and variances checked; Pythia and all source files unchanged.'}
    fout=ROOT.TFile(str(out/'weight_applied_overlay.root'),'RECREATE')
    for h in [*unscaled,*parts,total,oldtotal,py]:h.Write()
    draw(parts,total,py,out,'overlay_weight_applied_linear',False)
    draw(parts,total,py,out,'overlay_weight_applied_log',True)
    for path,digest in hashes.items():assert sha(Path(path))==digest,('Changed source',path)
    ROOT.TNamed('metadata',json.dumps(report)).Write();fout.Close()
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'output':str(out),'checks':checks,'verification':report['verification']},indent=2))

if __name__=='__main__':main()
