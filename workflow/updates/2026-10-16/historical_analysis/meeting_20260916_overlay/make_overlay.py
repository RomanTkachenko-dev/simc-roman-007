"""Bill meeting follow-up: frozen-scale overlays and separate weight checks.

Read existing samples only. No simulation, fit, new normalization, or sign repair.
"""
from pathlib import Path
import json, math, hashlib, shutil
import ROOT

BASE=Path('/mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN')
OUT=BASE/'meeting_20260916_overlay'
CHANNELS=['omega','phi','rho','eta','etap','xphsp','xphsp3']
TITLES=['#omega','#phi','#rho','#eta',"#eta'",'X-phase-space (legacy)','3#pi']
LABELS=['#omega','#phi / 10','#rho #times 2','#eta / 22.5',"#eta' / 15",'X-phase-space #times 4','3#pi #times 1 (test)']
COLORS=['#0000FF','#E07000','#00B000','#FF8000','#A040FF','#707070','#00AFC1']
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)
ROOT.TH1.SetDefaultSumw2(True)

def sha(p):
    with Path(p).open('rb') as f:
        h=hashlib.sha256()
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def close(a,b):
    assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-24),(a,b)

def gethist(path,key,name):
    f=ROOT.TFile.Open(str(path),'READ');assert f and not f.IsZombie(),str(path)
    h=f.Get(key);assert h and h.InheritsFrom('TH1'),key
    result=h.Clone(name);result.SetDirectory(0);f.Close();return result

def hist_equal(a,b):
    assert a.GetNbinsX()==b.GetNbinsX()
    for i in range(a.GetNbinsX()+2):
        close(a.GetBinContent(i),b.GetBinContent(i))
        close(a.GetBinError(i)**2,b.GetBinError(i)**2)

def save(canvas,name):
    for ext in ['png','pdf']:canvas.SaveAs(str(OUT/(name+'.'+ext)))
    canvas.Write();canvas.Close()

def overlay(pythia,parts,total,name,mode,log=False):
    c=ROOT.TCanvas('canvas_'+name,'',1450,880)
    c.SetLeftMargin(.13);c.SetRightMargin(.04);c.SetBottomMargin(.13);c.SetTopMargin(.07);c.SetTicks(1,1)
    c.SetLogy(log)
    frame=total.Clone('frame_'+name);frame.Reset('ICES');frame.SetTitle('')
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    frame.GetYaxis().SetTitle('Pythia-referenced density [GeV^{-1}]' if mode=='previous_reference' else 'Counts / 0.02 GeV (SIMC relative scales)')
    frame.GetYaxis().SetTitleOffset(1.45)
    peak=max(pythia.GetMaximum(),total.GetMaximum())
    frame.SetMinimum(.02 if log else 0);frame.SetMaximum(peak*(8 if log else 1.2));frame.Draw('AXIS')
    legend=ROOT.TLegend(.72,.44,.946,.91)
    legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.025)
    for h,label,color in zip(parts,LABELS,COLORS):
        h.SetLineColor(ROOT.TColor.GetColor(color));h.SetLineWidth(2);h.Draw('HIST SAME');legend.AddEntry(h,label,'l')
    total.SetLineColor(ROOT.TColor.GetColor('#FF0000'));total.SetLineWidth(3);total.Draw('HIST SAME')
    pythia.SetLineColor(ROOT.kBlack);pythia.SetLineWidth(3);pythia.Draw('HIST SAME')
    legend.AddEntry(total,'SIMC sum (7 channels)','l');legend.AddEntry(pythia,'Pythia Background','l');legend.Draw()
    note=ROOT.TLatex();note.SetNDC();note.SetTextFont(42);note.SetTextSize(.023)
    caption='Acceptance cuts; no event Weight; previous reference scale retained' if mode=='previous_reference' else 'Acceptance cuts; no event Weight; no Pythia normalization'
    note.DrawLatex(.13,.96,caption)
    c.RedrawAxis();save(c,name)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'overlay_results.root').exists(),'Refusing to overwrite completed output'
    seven=BASE/'seven_channel_sum_10600'
    old=BASE/'pythia_vs_simc_omega_peak_eta_reduced'
    pydir=BASE.parent/'BACKGROUND_ROOT_ANALYSIS/Pythia_Background'
    sm=json.loads((seven/'summary.json').read_text())
    om=json.loads((old/'summary.json').read_text())
    pm=json.loads((pydir/'summary.json').read_text())
    cut=sm['cuts_for_missing_mass'];assert cut==om['cuts']==pm['cut']
    assert not pm['has_weight']
    sources=[seven/'seven_channel_sum.root',seven/'summary.json',old/'omega_peak_matched.root',old/'summary.json',pydir/'histograms.root',pydir/'summary.json']
    sources += [BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root') for ch in CHANNELS]
    sources += [BASE/'runs'/(ch+'_Roman_007')/'plots/histograms.root' for ch in CHANNELS]
    hashes={str(p):sha(p) for p in sources}
    for p,digest in sm['sources'].items():assert sha(p)==digest,('Source changed',p)
    raw=[];weighted=[];stats={};all3=[]
    for ch in CHANNELS:
        p=BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root')
        f=ROOT.TFile.Open(str(p),'READ');t=f.Get('h10');assert t
        selection=ROOT.TTreeFormula('selection_'+ch,cut,t)
        h=ROOT.TH1D(ch+'_raw_after','',100,0,2)
        hw=ROOT.TH1D(ch+'_weighted_after','',100,0,2)
        accepted=[];allw=[]
        for i in range(t.GetEntries()):
            t.GetEntry(i);w=float(t.Weight);m=float(t.missmass)
            assert math.isfinite(w) and math.isfinite(m)
            allw.append(w)
            if selection.EvalInstance():h.Fill(m);hw.Fill(m,w);accepted.append(w)
        if ch=='xphsp3':all3=allw
        assert len(allw)==sm['channels'][ch]['before']['entries']
        assert len(accepted)==sm['channels'][ch]['after']['entries']
        histpath=p.parent/'plots/histograms.root'
        hist_equal(h,gethist(histpath,'missmass_raw_after','check_raw_'+ch))
        hist_equal(hw,gethist(histpath,'missmass_weighted_after','check_weighted_'+ch))
        close(hw.Integral(0,101),math.fsum(accepted))
        close(sum(hw.GetBinError(i)**2 for i in range(102)),math.fsum(w*w for w in accepted))
        raw.append(h);weighted.append(hw)
        stats[ch]={'before_cuts':len(allw),'after_cuts':len(accepted),'negative_before':sum(w<0 for w in allw),'negative_after':sum(w<0 for w in accepted),'min_weight':min(allw),'max_weight':max(allw),'sum_weight_after':math.fsum(accepted)}
        del selection;f.Close()
    parts=[]
    for ch,h in zip(CHANNELS,raw):
        hs=h.Clone(ch+'_relative_scaled_after');hs.Scale(sm['scales'][ch]);parts.append(hs)
    total=parts[0].Clone('sum_seven_relative_scaled_after');total.Reset('ICES')
    for h in parts:total.Add(h)
    hist_equal(total,gethist(seven/'seven_channel_sum.root','sum_seven_after','previous_seven_sum'))
    py=gethist(pydir/'histograms.root','mx_raw_after','pythia_counts_after')
    for h in [py,total,*parts]:
        assert h.GetNbinsX()==100
        for i in range(1,102):close(h.GetXaxis().GetBinLowEdge(i),total.GetXaxis().GetBinLowEdge(i))
    common=om['common_SIMC_factor_relative_to_scaled_counts']
    display=om['common_display_factor_for_both_spectra']
    ref_parts=[]
    for ch,h in zip(CHANNELS,parts):
        ref=h.Clone(ch+'_previous_reference_after');ref.Scale(common*display);ref_parts.append(ref)
        if ch!='xphsp3':
            hist_equal(ref,gethist(old/'omega_peak_matched.root',ch+'_peak_matched_density;1','old_reference_'+ch))
    ref_sum=ref_parts[0].Clone('sum_seven_previous_reference_after');ref_sum.Reset('ICES')
    for h in ref_parts:ref_sum.Add(h)
    ref_py=py.Clone('pythia_previous_reference_after');ref_py.Scale(display)
    hist_equal(ref_py,gethist(old/'omega_peak_matched.root','pythia_reference_density;1','old_pythia'))
    for group,hs in [(total,parts),(ref_sum,ref_parts)]:
        for i in range(102):
            close(group.GetBinContent(i),math.fsum(h.GetBinContent(i) for h in hs))
            close(group.GetBinError(i)**2,math.fsum(h.GetBinError(i)**2 for h in hs))
    wh=gethist(seven/'seven_channel_sum.root','xphsp3_original_Weight_all_events','three_pi_Weight_no_analysis_cuts')
    check=wh.Clone('three_pi_weight_check');check.Reset('ICES')
    for w in all3:check.Fill(w)
    hist_equal(wh,check);assert wh.Integral()==1000 and min(all3)>0
    output=ROOT.TFile(str(OUT/'overlay_results.root'),'RECREATE')
    for h in [*raw,*weighted,*parts,total,py,*ref_parts,ref_sum,ref_py,wh]:h.Write()
    overlay(py,parts,total,'overlay_no_pythia_normalization_linear','relative_only')
    overlay(py,parts,total,'overlay_no_pythia_normalization_log','relative_only',True)
    overlay(ref_py,ref_parts,ref_sum,'overlay_previous_scale_no_refit','previous_reference')
    # Reuse the established Weight diagnostic formatting, but freshly verify events.
    import sys
    sys.path.insert(0,str(BASE))
    from plot_weight_distributions import panel, describe
    c=ROOT.TCanvas('canvas_three_pi_weight','',1050,750)
    refs=panel(c,wh,'3#pi: original Weight, no analysis cuts',describe(all3))
    save(c,'three_pi_Weight_no_analysis_cuts')
    # Apply original Weight only to the explicitly requested positive-weight samples.
    c=ROOT.TCanvas('canvas_positive_channels_weighted','',1650,900);c.Divide(3,2,.003,.003)
    keep=[]
    for col,index in enumerate([0,5,6],1):
        ch=CHANNELS[index];assert stats[ch]['negative_before']==0
        for row,h0 in enumerate([raw[index],weighted[index]]):
            pad=c.cd(col+3*row);pad.SetLeftMargin(.17);pad.SetBottomMargin(.15);pad.SetTopMargin(.12);pad.SetTicks(1,1)
            h=h0.Clone(ch+'_diagnostic_'+str(row));h.SetLineColor(ROOT.kRed if row else ROOT.kBlue)
            h.SetLineWidth(2);h.SetMinimum(0);h.SetMaximum(h.GetMaximum()*1.25)
            h.GetXaxis().SetTitle('Missing mass [GeV]');h.GetYaxis().SetTitle('Sum of Weight / 0.02 GeV' if row else 'Events / 0.02 GeV')
            h.GetYaxis().SetTitleOffset(1.6);h.Draw('HIST')
            title=ROOT.TLatex();title.SetNDC();title.SetTextFont(42);title.SetTextSize(.045)
            title.DrawLatex(.17,.94,TITLES[index]+(' - with Weight' if row else ' - no Weight'))
            keep.extend([h,title])
    save(c,'omega_legacy_3pi_weighted_after_cuts')
    report={'source_hashes':hashes,'cuts':cut,'channels':stats,'relative_scales_retained':sm['scales'],
        'raw_overlay':'Event Weight is NOT applied to any of the seven components. Existing relative channel factors retained; 3pi coefficient remains provisional 1. Pythia is a separate curve, never included in the SIMC sum.',
        'no_pythia_normalization':'Pythia counts unchanged; only previously established relative SIMC factors applied. Log version exposes the different sample sizes; no unit-area scaling.',
        'previous_scale_overlay':{'common_SIMC_factor_frozen':common,'common_display_factor_frozen':display,'new_fit_or_peak_matching':False,'description':'All six previous component curves and Pythia exactly match the previous reference-density result. New 3pi has coefficient 1 and the same frozen common SIMC factor. The seven-channel sum is not re-matched to the omega peak.'},
        'weight_diagnostics':'3pi original Weight distribution uses all stored events and unit-count fills; no analysis cuts. Separate omega, legacy phase space and 3pi missing-mass histograms apply original Weight after the common cuts, without arbitrary channel factors or Pythia normalization.',
        'variable_caveat':'Pythia mx treated provisionally as reconstructed missing mass in GeV. Equivalence to SIMC missmass and correspondence to Henry analyzed/thin file remain unverified.',
        'physics_caveat':'Diagnostic comparisons, not cross sections or yields. Raw 3pi does not include its Weight-carried phase-space factor. Negative weights unchanged. No hybrid weighted/unweighted total, fit, luminosity conversion, or source simulation rerun.',
        'verification':'Seven source trees refilled; raw and weighted bins and variances match stored histograms; cuts/binning agree; raw seven-channel sum unchanged; historical six-component and Pythia reference curves unchanged; source hashes unchanged.'}
    assert all(sha(p)==digest for p,digest in hashes.items())
    ROOT.TNamed('analysis_metadata',json.dumps(report)).Write();output.Close()
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    shutil.copy2(__file__,OUT/'make_overlay.py')
    print(json.dumps({'output':str(OUT),'three_pi':stats['xphsp3'],'old_scale':common,'checks':report['verification']},indent=2))

if __name__=='__main__':main()
