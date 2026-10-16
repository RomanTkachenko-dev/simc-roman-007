"""User-requested scaling of raw reconstructed mass plots, with an overlaid sum."""
from pathlib import Path
import argparse, hashlib, json, math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
CHANNELS=['omega','phi','rho','eta','etap','xphsp']
SCALES=dict(omega=1.,phi=.1,rho=2.,eta=.1,etap=.1,xphsp=2.)
LABELS=['#omega', '#phi #times 0.1', '#rho #times 2', '#eta #times 0.1', "#eta' #times 0.1", 'X-phase-space #times 2']
COLORS=[ROOT.kBlue+1,ROOT.kRed+1,ROOT.kGreen+2,ROOT.kOrange+7,ROOT.kViolet+1,ROOT.kGray+2]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs',type=Path,default=Path(__file__).resolve().parent/'runs')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'scaled_missing_mass')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    report={'factors':SCALES,'input_histogram':'missmass_raw_after','stage':'After SIMC and requested acceptance cuts',
        'normalization':'User-requested display factors applied to raw test-sample histograms; not luminosity or cross-section normalization.',
        'weights':'Original per-event Weight is not used in this raw-count plot. No event data modified.',
        'sources':{}}
    histograms=[]; cuts=[]
    for channel,color in zip(CHANNELS,COLORS):
        folder=args.runs/(channel+'_Roman_007')/'plots'
        path=folder/'histograms.root'
        before_sha=sha(path)
        source=ROOT.TFile.Open(str(path),'READ')
        if not source or source.IsZombie():raise RuntimeError('Cannot open '+str(path))
        original=source.Get('missmass_raw_after')
        assert original and original.GetNbinsX()==100
        assert original.GetXaxis().GetXmin()==0 and original.GetXaxis().GetXmax()==2
        h=original.Clone(channel+'_scaled_missing_mass')
        factor=SCALES[channel]
        h.Scale(factor)
        for i in range(102):
            assert math.isclose(h.GetBinContent(i),original.GetBinContent(i)*factor,rel_tol=1e-12,abs_tol=1e-12)
            assert math.isclose(h.GetBinError(i),original.GetBinError(i)*abs(factor),rel_tol=1e-12,abs_tol=1e-12)
        h.SetLineColor(color); h.SetLineWidth(3); h.SetFillStyle(0)
        summary=json.loads((folder/'summary.json').read_text())
        cuts.append(summary['cut'])
        report['sources'][channel]={'path':str(path),'sha256':before_sha,'raw_entries_after_cuts':summary['after'],
            'factor':factor,'scaled_integral_in_range':h.Integral(),'scaled_integral_including_flow':h.Integral(0,101),
            'scaled_underflow':h.GetBinContent(0),'scaled_overflow':h.GetBinContent(101)}
        source.Close()
        assert sha(path)==before_sha,'Source file changed'
        histograms.append(h)
    assert len(set(cuts))==1
    report['cuts']=cuts[0]
    total=histograms[0].Clone('sum_scaled_missing_mass'); total.Reset('ICES')
    for h in histograms:total.Add(h)
    for i in range(102):
        assert math.isclose(total.GetBinContent(i),sum(h.GetBinContent(i) for h in histograms),rel_tol=1e-12,abs_tol=1e-12)
        assert math.isclose(total.GetBinError(i)**2,sum(h.GetBinError(i)**2 for h in histograms),rel_tol=1e-12,abs_tol=1e-12)
    report['sum']={'in_range':total.Integral(),'including_flow':total.Integral(0,101),
        'underflow':total.GetBinContent(0),'overflow':total.GetBinContent(101)}
    report['verification']='All bins, errors, sum bins and sum variances verified; six source files unchanged.'
    out=ROOT.TFile(str(args.output/'scaled_missing_mass.root'),'RECREATE')
    for h in histograms:h.Write()
    total.Write()
    canvas=ROOT.TCanvas('scaled_overlay_and_sum','',1450,880)
    canvas.SetLeftMargin(.12); canvas.SetRightMargin(.04); canvas.SetBottomMargin(.12); canvas.SetTopMargin(.045)
    canvas.SetTicks(1,1)
    frame=total.Clone('frame'); frame.Reset('ICES'); frame.SetTitle('')
    frame.GetXaxis().SetTitle('Reconstructed missing mass [GeV]')
    frame.GetYaxis().SetTitle('Scaled test-sample events / 0.02 GeV')
    frame.GetXaxis().SetTitleSize(.04); frame.GetYaxis().SetTitleSize(.04); frame.GetYaxis().SetTitleOffset(1.4)
    frame.SetMinimum(0); frame.SetMaximum(total.GetMaximum()*1.16); frame.Draw('AXIS')
    for h in histograms:h.Draw('HIST SAME')
    total.SetLineColor(ROOT.kBlack); total.SetLineWidth(3); total.Draw('HIST SAME')
    legend=ROOT.TLegend(.735,.535,.945,.93)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(.026)
    for h,label in zip(histograms,LABELS):legend.AddEntry(h,label,'l')
    legend.AddEntry(total,'Sum','l'); legend.Draw()
    canvas.RedrawAxis()
    for ext in ['png','pdf']:canvas.SaveAs(str(args.output/('missing_mass_scaled_overlay_sum_after.'+ext)))
    canvas.Write()
    ROOT.TNamed('cuts',cuts[0]).Write()
    ROOT.TNamed('scale_factors',json.dumps(SCALES)).Write()
    ROOT.TNamed('normalization_warning',report['normalization']).Write()
    out.Close();canvas.Close()
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'factors':SCALES,'sum':report['sum'],'verification':report['verification']},indent=2))

if __name__=='__main__':main()
