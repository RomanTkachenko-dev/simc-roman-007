"""Compare Pythia Background mx with the current scaled SIMC mass sum."""
from pathlib import Path
import argparse, hashlib, json, math
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
CHANNELS=['omega','phi','rho','eta','etap','xphsp']

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def render(pythia,total,components,labels,destination,name,normalized,detailed=False,y_title=None,sum_color='#00AFC1'):
    c=ROOT.TCanvas(name,'',1450,880)
    c.SetLeftMargin(.12);c.SetRightMargin(.04);c.SetBottomMargin(.12);c.SetTopMargin(.045);c.SetTicks(1,1)
    frame=total.Clone(name+'_frame');frame.Reset('ICES');frame.SetTitle('')
    # Axis units and equality of mx with reconstructed mass are provisional;
    # the supplied writer is unavailable. This assumption is recorded in metadata.
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    frame.GetYaxis().SetTitle('Normalized density [GeV^{-1}]' if normalized else 'Entries / 0.02 GeV (SIMC scaled)')
    if y_title:frame.GetYaxis().SetTitle(y_title)
    frame.GetXaxis().SetTitleSize(.04);frame.GetYaxis().SetTitleSize(.04);frame.GetYaxis().SetTitleOffset(1.4)
    frame.SetMinimum(0);frame.SetMaximum(max(pythia.GetMaximum(),total.GetMaximum())*1.16);frame.Draw('AXIS')
    if detailed:
        for h in components:h.SetLineWidth(2);h.Draw('HIST SAME')
    total.SetLineColor(ROOT.TColor.GetColor(sum_color));total.SetLineWidth(3);total.SetLineStyle(1)
    pythia.SetLineColor(ROOT.kBlack);pythia.SetLineWidth(3);pythia.SetLineStyle(1)
    total.Draw('HIST SAME');pythia.Draw('HIST SAME')
    legend=ROOT.TLegend(.73,.48 if detailed else .78,.947,.93)
    legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.026)
    if detailed:
        for h,label in zip(components,labels):legend.AddEntry(h,label,'l')
    legend.AddEntry(total,'Scaled SIMC sum','l');legend.AddEntry(pythia,'Pythia Background','l');legend.Draw()
    c.RedrawAxis()
    for ext in ['png','pdf']:c.SaveAs(str(destination/(name+'.'+ext)))
    c.Write();c.Close()

def main():
    base=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--simc-dir',type=Path,default=base/'scaled_missing_mass')
    p.add_argument('--pythia-dir',type=Path,default=base.parent/'BACKGROUND_ROOT_ANALYSIS/Pythia_Background')
    p.add_argument('--output',type=Path,default=base/'pythia_vs_scaled_simc')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    simc_meta=json.loads((args.simc_dir/'summary.json').read_text())
    pythia_meta=json.loads((args.pythia_dir/'summary.json').read_text())
    assert simc_meta['cuts']==pythia_meta['cut'],'Cuts differ'
    simc_path=args.simc_dir/'scaled_missing_mass.root';pythia_path=args.pythia_dir/'histograms.root'
    digests={str(path):sha(path) for path in [simc_path,pythia_path]}
    sf=ROOT.TFile.Open(str(simc_path),'READ');pf=ROOT.TFile.Open(str(pythia_path),'READ')
    if not sf or sf.IsZombie() or not pf or pf.IsZombie():raise RuntimeError('Cannot open input ROOT files')
    total=sf.Get('sum_scaled_missing_mass').Clone('simc_sum_unnormalized')
    pythia=pf.Get('mx_raw_after').Clone('pythia_background_unnormalized')
    components=[sf.Get(c+'_scaled_missing_mass').Clone(c+'_unnormalized') for c in CHANNELS]
    sf.Close();pf.Close()
    for h in [total,pythia,*components]:
        assert h.GetNbinsX()==100 and h.GetXaxis().GetXmin()==0 and h.GetXaxis().GetXmax()==2
        assert all(math.isfinite(h.GetBinContent(i)) and h.GetBinContent(i)>=0 for i in range(102))
    for i in range(102):
        assert math.isclose(total.GetBinContent(i),sum(h.GetBinContent(i) for h in components),rel_tol=1e-12,abs_tol=1e-12)
    simc_area=total.Integral('width');pythia_area=pythia.Integral('width')
    assert simc_area>0 and pythia_area>0
    nt=total.Clone('simc_sum_unit_area');np=pythia.Clone('pythia_background_unit_area')
    nt.Scale(1./simc_area);np.Scale(1./pythia_area)
    nc=[h.Clone(CHANNELS[i]+'_common_normalization') for i,h in enumerate(components)]
    for h in nc:h.Scale(1./simc_area)
    assert math.isclose(nt.Integral('width'),1.,abs_tol=1e-12)
    assert math.isclose(np.Integral('width'),1.,abs_tol=1e-12)
    for i in range(102):
        assert math.isclose(nt.GetBinContent(i),sum(h.GetBinContent(i) for h in nc),rel_tol=1e-12,abs_tol=1e-12)
        for h,orig,factor in [(nt,total,1./simc_area),(np,pythia,1./pythia_area)]:
            assert math.isclose(h.GetBinError(i),orig.GetBinError(i)*factor,rel_tol=1e-12,abs_tol=1e-12)
    factors=simc_meta['factors']
    titles=['#omega','#phi','#rho','#eta',"#eta'",'X-phase-space']
    labels=[label if factors[ch]==1 else f'{label} #times {factors[ch]:.4g}' for ch,label in zip(CHANNELS,titles)]
    report={'stage':'After SIMC / supplied ROOT, after matching acceptance cuts','cuts':simc_meta['cuts'],
        'sources':digests,'pythia_input':{'original_file':pythia_meta['file'],'tree':pythia_meta['tree'],'branch':'mx'},
        'simc_input':{'branch':'missmass','channel_factors':factors},
        'variable_assumption':'For this provisional visual comparison, supplied Pythia mx is treated as a reconstructed missing mass in GeV. The writer/definition was not provided; physical equivalence is not independently verified.',
        'normalization':'Each total separately has integral(density dM)=1 over [0,2) GeV. Underflow and overflow excluded from denominators. All six SIMC components share the SAME SIMC normalization factor; their relative scale factors are preserved.',
        'physics_warning':'Shape comparison only, not physical signal/background yields; no fits or cross-section normalization. Pythia is not added to the six-channel SIMC sum.',
        'binning':{'bins':100,'range_GeV':[0,2],'width_GeV':.02},
        'simc':{'raw_in_range':total.Integral(),'raw_underflow':total.GetBinContent(0),'raw_overflow':total.GetBinContent(101),'density_multiplier':1./simc_area,'normalized_area':nt.Integral('width')},
        'pythia':{'raw_in_range':pythia.Integral(),'raw_underflow':pythia.GetBinContent(0),'raw_overflow':pythia.GetBinContent(101),'density_multiplier':1./pythia_area,'normalized_area':np.Integral('width')},
        'verification':'Matching cuts and bins; both total areas=1; sum of normalized components matches normalized SIMC total in every bin; statistical errors scaled; sources unchanged.'}
    out=ROOT.TFile(str(args.output/'pythia_vs_scaled_simc.root'),'RECREATE')
    for h in [total,pythia,*components,nt,np,*nc]:h.Write()
    render(pythia,total,components,labels,args.output,'comparison_unnormalized_after',False)
    render(np,nt,nc,labels,args.output,'comparison_normalized_after',True)
    render(np,nt,nc,labels,args.output,'comparison_normalized_channels_after',True,True)
    ROOT.TNamed('analysis_metadata',json.dumps(report)).Write();out.Close()
    assert all(sha(Path(path))==value for path,value in digests.items()),'Source file changed'
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'simc':report['simc'],'pythia':report['pythia'],'verification':report['verification']},indent=2))

if __name__=='__main__':main()
