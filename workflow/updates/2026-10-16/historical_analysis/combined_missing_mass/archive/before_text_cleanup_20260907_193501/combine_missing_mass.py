"""Overlay and sum the exact saved post-SIMC histograms; no physics rescaling."""
from pathlib import Path
import argparse, hashlib, json, math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)
ROOT.TH1.AddDirectory(False)
CHANNELS = ['omega', 'phi', 'rho', 'eta', 'etap', 'xphsp']
LABELS = ['#omega', '#phi', '#rho', '#eta', "#eta'", 'X-phase-space']
COLORS = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kViolet+1, ROOT.kGray+2]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def canvas_frame(name, title, state, mode):
    c=ROOT.TCanvas(name,title,1450,880)
    c.SetLeftMargin(.12); c.SetRightMargin(.29); c.SetBottomMargin(.12); c.SetTopMargin(.17)
    c.SetTicks(1,1)
    text=ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42)
    return c,text

def headers(text,title,state,mode):
    text.SetTextSize(.032); text.DrawLatex(.12,.955,title)
    text.SetTextSize(.023)
    stage='After requested acceptance cuts' if state=='after' else 'Before additional analysis cuts'
    text.DrawLatex(.12,.914,'Reconstructed after SIMC | '+stage+' | PROVISIONAL')
    note='Test-sample counts only; no relative channel normalization' if mode=='raw' else 'Signed Weight retained; no relative channel normalization; diagnostic only'
    text.SetTextSize(.021); text.DrawLatex(.12,.876,note)

def axes(h, low, high, mode):
    h.SetTitle('')
    h.GetXaxis().SetTitle('Reconstructed missing mass [GeV]')
    h.GetYaxis().SetTitle('Test-sample events / 0.02 GeV' if mode=='raw' else 'Sum of signed Weight / 0.02 GeV')
    h.GetXaxis().SetTitleSize(.04); h.GetYaxis().SetTitleSize(.04)
    h.GetYaxis().SetTitleOffset(1.4)
    h.SetMinimum(low); h.SetMaximum(high)

def draw(hists,total,path,state,mode,kind,counts):
    name=f'missing_mass_{kind}_{mode}_{state}'
    title='All six channels: missing-mass overlay' if kind=='overlay' else 'All six channels: missing-mass sum'
    c,text=canvas_frame(name,title,state,mode)
    stacked=kind=='sum' and mode=='raw'
    shown=hists if kind=='overlay' else [total]
    lo=min(0.,min(h.GetMinimum() for h in shown))
    hi=max(0.,max(h.GetMaximum() for h in shown))
    span=hi-lo or 1.
    lo=lo-.12*span if lo<0 else 0.
    hi=hi+.16*span
    frame=total.Clone(name+'_frame'); frame.Reset('ICES')
    axes(frame,lo,hi,mode); frame.Draw('AXIS')
    legend=ROOT.TLegend(.745,.48,.98,.82)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(.026)
    stack=None
    if stacked:
        stack=ROOT.THStack(name+'_stack','')
        for index in [5,2,0,3,4,1]:
            h=hists[index]
            h.SetFillStyle(1001); h.SetFillColorAlpha(COLORS[index],.65); h.SetLineWidth(1)
            stack.Add(h)
        stack.Draw('HIST SAME')
    elif kind=='overlay':
        for h in hists: h.Draw('HIST SAME')
    if kind=='sum':
        total.SetLineColor(ROOT.kBlack); total.SetLineWidth(3); total.Draw('HIST SAME')
    if kind=='overlay' or stacked:
        for h,label,count in zip(hists,LABELS,counts):
            legend.AddEntry(h,f'{label}  (N={count})','f' if stacked else 'l')
    if kind=='sum': legend.AddEntry(total,'Sum of six test samples','l')
    legend.Draw()
    text.SetTextSize(.023)
    text.DrawLatex(.745,.41,'No cross-section scaling')
    text.DrawLatex(.745,.375,'No luminosity scaling')
    text.DrawLatex(.745,.34,'X-phase-space #neq verified 2#pi')
    if mode=='weighted':
        text.SetTextColor(ROOT.kRed+1)
        text.DrawLatex(.745,.285,'Negative weights present')
        text.DrawLatex(.745,.25,'Not a physical total yield')
        text.SetTextColor(ROOT.kBlack)
    elif kind=='sum':
        text.DrawLatex(.745,.285,f'Total entries: {sum(counts)}')
        text.DrawLatex(.745,.25,'Arbitrary test-sample mixture')
        text.DrawLatex(.745,.215,f'Displayed: {int(total.Integral())}; M<0: {int(total.GetBinContent(0))}')
    zero=ROOT.TLine(0,0,2,0); zero.SetLineColor(ROOT.kGray+1); zero.SetLineStyle(2)
    if lo<0: zero.Draw()
    c.RedrawAxis()
    headers(text,title,state,mode)
    for extension in ['png','pdf']:c.SaveAs(str(path/(name+'.'+extension)))
    c.Write(name)
    c.Close()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,default=Path(__file__).resolve().parent/'runs')
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'combined_missing_mass')
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    report={'normalization':'Scale=1 for each test sample. No cross-section/luminosity normalization; not a physical combined yield.',
        'weight_policy':'Signed Weight preserved; no abs(Weight), clipping or unit-area renormalization.',
        'binning':{'bins':100,'min_GeV':0,'max_GeV':2,'width_GeV':.02},'sources':{},'results':{}}
    sources=[]; summaries=[]; cuts=[]
    for channel in CHANNELS:
        folder=args.runs/(channel+'_Roman_007')/'plots'
        source=folder/'histograms.root'
        f=ROOT.TFile.Open(str(source),'READ')
        if not f or f.IsZombie(): raise RuntimeError('Cannot open '+str(source))
        sources.append(f)
        summary=json.loads((folder/'summary.json').read_text())
        summaries.append(summary); cuts.append(summary['cut'])
        report['sources'][channel]={'histograms':str(source),'sha256':digest(source),'scale':1.,
            'entries_before':summary['before'],'entries_after':summary['after']}
    assert len(set(cuts))==1,'Input histograms used different cuts'
    report['cut']=cuts[0]
    output=ROOT.TFile(str(args.output/'combined_missing_mass.root'),'RECREATE')
    for state in ['before','after']:
        for mode in ['raw','weighted']:
            histograms=[]
            for source,channel,color in zip(sources,CHANNELS,COLORS):
                original=source.Get(f'missmass_{mode}_{state}')
                assert original and original.GetNbinsX()==100
                assert original.GetXaxis().GetXmin()==0 and original.GetXaxis().GetXmax()==2
                h=original.Clone(f'{channel}_missmass_{mode}_{state}')
                h.SetLineColor(color); h.SetLineWidth(3); h.SetFillStyle(0)
                assert all(math.isfinite(h.GetBinContent(i)) for i in range(102))
                histograms.append(h)
            total=histograms[0].Clone(f'sum_missmass_{mode}_{state}')
            total.Reset('ICES')
            for h in histograms: total.Add(h)
            for i in range(102):
                assert math.isclose(total.GetBinContent(i),sum(h.GetBinContent(i) for h in histograms),rel_tol=1e-12,abs_tol=1e-20)
                assert math.isclose(total.GetBinError(i)**2,sum(h.GetBinError(i)**2 for h in histograms),rel_tol=1e-12,abs_tol=1e-25)
            counts=[s[state] for s in summaries]
            if mode=='raw': assert total.Integral(0,101)==sum(counts)
            report['results'][f'{mode}_{state}']={'entries':sum(counts),'in_range':total.Integral(),
                'underflow':total.GetBinContent(0),'overflow':total.GetBinContent(101),
                'including_flow':total.Integral(0,101),'sum_bins_and_variances_verified':True}
            output.cd()
            for h in histograms:h.Write()
            total.Write()
            for kind in ['overlay','sum']:draw(histograms,total,args.output,state,mode,kind,counts)
    output.cd()
    ROOT.TNamed('cuts',cuts[0]).Write()
    ROOT.TNamed('normalization_warning',report['normalization']).Write()
    output.Close()
    for f in sources:f.Close()
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))

if __name__=='__main__': main()
