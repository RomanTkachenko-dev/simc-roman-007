"""Plot original event Weight as the x variable, using unit-count histogram fills."""
import json
import math
from pathlib import Path
import ROOT
from compare_pythia_scaled_sum import CHANNELS, sha
from normalize_omega_peak import BASE

OUT = BASE/'weight_distributions_after_cuts'
TITLES = ['#omega', '#phi', '#rho', '#eta', "#eta'", 'X-phase-space']
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)


def describe(values):
    return {'entries':len(values), 'negative':sum(w<0 for w in values),
            'positive':sum(w>0 for w in values), 'zero':sum(w==0 for w in values),
            'min':min(values), 'max':max(values), 'mean':math.fsum(values)/len(values),
            'sum_weights':math.fsum(values)}


def panel(pad, hist, label, counts):
    pad.cd(); pad.SetLeftMargin(.14); pad.SetRightMargin(.065)
    pad.SetBottomMargin(.17); pad.SetTopMargin(.14); pad.SetTicks(1,1)
    # Put the power of ten in the axis title, avoiding ROOT's clipped exponent
    # at the right edge of each small pad. Saved data histogram stays unmodified.
    exponent=math.floor(math.log10(max(abs(counts['min']),abs(counts['max']))))
    factor=10.**exponent
    hist=hist.Clone(hist.GetName()+'_display_'+pad.GetName())
    hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin()/factor,hist.GetXaxis().GetXmax()/factor)
    hist.GetXaxis().SetNoExponent(True)
    hist.SetLineColor(ROOT.kBlue+1); hist.SetLineWidth(2); hist.SetMinimum(0)
    hist.SetMaximum(hist.GetMaximum()*1.18)
    hist.GetXaxis().SetTitle(f'Weight [#times10^{{{exponent}}}]'); hist.GetYaxis().SetTitle('Events / bin')
    for ax in [hist.GetXaxis(),hist.GetYaxis()]:
        ax.SetTitleSize(.05); ax.SetLabelSize(.041); ax.SetMaxDigits(3)
    hist.GetXaxis().SetTitleOffset(1.12); hist.GetYaxis().SetTitleOffset(1.22)
    hist.Draw('HIST')
    line = ROOT.TLine(0,0,0,hist.GetMaximum()); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+2); line.Draw()
    title=ROOT.TLatex(); title.SetNDC(); title.SetTextFont(42); title.SetTextAlign(22); title.SetTextSize(.064)
    title.DrawLatex(.54,.94,label)
    note=ROOT.TLatex(); note.SetNDC(); note.SetTextFont(42); note.SetTextAlign(31); note.SetTextSize(.035)
    note.DrawLatex(.92,.805,f"N = {counts['entries']}; Weight < 0: {counts['negative']}")
    pad.RedrawAxis()
    return [hist,line,title,note]


def main():
    OUT.mkdir(exist_ok=True)
    previous=json.loads((BASE/'pythia_vs_simc_weighted/summary.json').read_text())
    cuts=previous['cuts']
    report={'cuts':cuts, 'expression':'Weight', 'selection':cuts,
            'histogram_fill_weight':1, 'channel_scale_factors_applied':False,
            'normalization':'None: x is original signed Weight, y is number of events. No abs, clipping or channel scaling.',
            'binning':'80 bins per channel; a channel-specific range includes all accepted values and zero with padding.',
            'channels':{}, 'sources_sha256':{}}
    histograms=[]
    for channel in CHANNELS:
        path=BASE/'runs'/(channel+'_Roman_007')/(channel+'_Roman_007.root')
        report['sources_sha256'][str(path)]=sha(path)
        f=ROOT.TFile.Open(str(path),'READ'); assert f and not f.IsZombie()
        t=f.Get('h10'); assert t and t.GetBranch('Weight')
        select=ROOT.TTreeFormula('selection_'+channel,cuts,t)
        all_values=[]; values=[]
        for i in range(t.GetEntries()):
            t.GetEntry(i); w=float(t.Weight); assert math.isfinite(w)
            all_values.append(w)
            if select.EvalInstance():values.append(w)
        assert len(values)==t.GetEntries(cuts)
        counts=describe(values)
        assert counts['entries']==previous['channels'][channel]['accepted']
        assert counts['negative']==previous['channels'][channel]['negative_accepted']
        low=min(0.,min(values)); high=max(0.,max(values)); span=high-low
        assert span>0
        low-=.04*span; high+=.04*span
        h=ROOT.TH1D(channel+'_Weight_after_cuts','',80,low,high); h.Sumw2()
        for w in values:h.Fill(w)
        assert h.GetEntries()==len(values) and h.Integral()==len(values)
        assert h.GetBinContent(0)==0 and h.GetBinContent(81)==0
        assert math.isclose(sum(h.GetBinError(i)**2 for i in range(1,81)),len(values),rel_tol=1e-12)
        report['channels'][channel]={'file':str(path),'tree':'h10','before_cuts':describe(all_values),
                                     'after_cuts':counts,'histogram_range':[low,high]}
        histograms.append(h)
        del select
        f.Close()
    out=ROOT.TFile(str(OUT/'weight_distributions.root'),'RECREATE')
    canvas=ROOT.TCanvas('all_six_weight_distributions','',1800,1100); canvas.Divide(3,2,.001,.001)
    keep=[]
    for i,(channel,h,title) in enumerate(zip(CHANNELS,histograms,TITLES),1):
        keep.extend(panel(canvas.cd(i),h,title,report['channels'][channel]['after_cuts']))
    for ext in ['png','pdf']:canvas.SaveAs(str(OUT/('all_six_Weight_after_cuts.'+ext)))
    canvas.Write(); canvas.Close()
    for channel,h,title in zip(CHANNELS,histograms,TITLES):
        c=ROOT.TCanvas('canvas_'+channel+'_Weight','',1050,750)
        refs=panel(c,h,title,report['channels'][channel]['after_cuts'])
        for ext in ['png','pdf']:c.SaveAs(str(OUT/(channel+'_Weight_after_cuts.'+ext)))
        h.Write(); c.Write(); c.Close()
    report['verification']='All six trees read; accepted and negative counts match independent prior audit; histogram counts and variances verified; zero underflow/overflow; sources unchanged.'
    assert all(sha(Path(p))==digest for p,digest in report['sources_sha256'].items())
    ROOT.TNamed('metadata',json.dumps(report)).Write(); out.Close()
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({ch:info['after_cuts'] for ch,info in report['channels'].items()},indent=2))


if __name__=='__main__':main()
