"""Analyze matched-setting 3pi and extend the provisional raw-count SIMC sum."""
from pathlib import Path
import argparse
import json
import math
import ROOT
from analyze_root import analyze, CUT
from plot_weight_distributions import panel, describe
from compare_pythia_scaled_sum import CHANNELS, sha
from plot_scaled_missing_mass import COLORS

BASE=Path(__file__).resolve().parent
NEW='xphsp3'
TITLES=['#omega','#phi','#rho','#eta',"#eta'",'X-phase-space (legacy)','X-phase-space 3#pi']


def close(a,b):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-25),(a,b)


def sum_canvas(dest,stage,hs,total,six,labels,comparison=False):
    name=('sum_six_vs_seven_' if comparison else 'seven_channel_sum_')+stage
    c=ROOT.TCanvas('canvas_'+name,'',1450,880)
    c.SetLeftMargin(.13);c.SetRightMargin(.04);c.SetBottomMargin(.13);c.SetTopMargin(.045);c.SetTicks(1,1)
    frame=total.Clone('frame_'+name);frame.Reset('ICES');frame.SetTitle('')
    frame.GetXaxis().SetTitle('Reconstructed missing mass [GeV]')
    frame.GetYaxis().SetTitle('Scaled test-sample events / 0.02 GeV')
    frame.GetXaxis().SetTitleSize(.04);frame.GetYaxis().SetTitleSize(.04);frame.GetYaxis().SetTitleOffset(1.5)
    frame.SetMinimum(0);frame.SetMaximum(total.GetMaximum()*1.2);frame.Draw('AXIS')
    leg=ROOT.TLegend(.70,.72 if comparison else .49,.945,.93)
    leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextFont(42);leg.SetTextSize(.026)
    if comparison:
        six.SetLineColor(ROOT.kGray+2);six.SetLineWidth(3);six.Draw('HIST SAME')
        hs[-1].Draw('HIST SAME');leg.AddEntry(six,'Previous 6-channel sum','l');leg.AddEntry(hs[-1],labels[-1],'l')
    else:
        for h,label in zip(hs,labels):h.Draw('HIST SAME');leg.AddEntry(h,label,'l')
    total.SetLineColor(ROOT.TColor.GetColor('#FF0000'));total.SetLineWidth(3);total.Draw('HIST SAME')
    leg.AddEntry(total,'Sum (7 channels)','l');leg.Draw();c.RedrawAxis()
    for ext in ['png','pdf']:c.SaveAs(str(dest/(name+'.'+ext)))
    c.Write();c.Close()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--three-pion-scale',type=float,default=1.)
    ap.add_argument('--output',type=Path,default=BASE/'seven_channel_sum_10600')
    args=ap.parse_args();assert math.isfinite(args.three_pion_scale) and args.three_pion_scale>=0
    dest=args.output;dest.mkdir(exist_ok=True,parents=True)
    channels=CHANNELS+[NEW]
    input_files=[BASE/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.inp') for ch in channels]
    apparatus=[p.read_text().split('begin parm kinematics_main',1)[1] for p in input_files]
    assert len(set(apparatus))==1,'Apparatus/simulation settings differ'
    assert 'doing_Xphasespace = 3' in input_files[-1].read_text()
    old_meta=json.loads((BASE/'pythia_vs_simc_omega_peak_eta_reduced/summary.json').read_text())
    assert old_meta['cuts']==CUT
    scales=dict(old_meta['relative_channel_factors']);scales[NEW]=args.three_pion_scale
    new_dir=BASE/'runs/xphsp3_Roman_007';new_path=new_dir/'xphsp3_Roman_007.root'
    new_summary=analyze(new_path,'h10',new_dir/'plots','3pi SIMC, 10.60 GeV',True)
    assert new_summary['before']==1000
    ROOT.gStyle.SetOptTitle(0);ROOT.gStyle.SetOptStat(0);ROOT.TH1.AddDirectory(False)
    report={'cuts_for_missing_mass':CUT,'weight_distribution_selection':'All 1000 stored events per channel, no additional analysis cuts.',
        'kinematics':{'Ebeam_GeV':10.6,'electron_arm':'SHMS','electron_p_GeV':6.7,'electron_theta_deg':13.,'proton_arm':'HMS','proton_p_GeV':1.1,'proton_theta_deg':32.,
                      'radiation':False,'energy_loss':False,'smearing':False},
        'input_verification':'All seven files identical from kinematics_main through end of file; new mode is doing_Xphasespace=3.',
        'sum_method':'Sum of raw missing-mass counts with the previously requested channel factors. No event Weight and no Pythia scaling. Produced before and after analysis cuts.',
        'scales':scales,'three_pion_scale_status':'Provisional display choice, not a physical normalization or supervisor-approved yield.',
        'physics_warning':'Weights of several original channels remain under review. The raw-count sum is diagnostic only; it does NOT include the new three-body phase-space factor carried by Weight. The standalone new-channel weighted histograms retain that factor. Relative channel yields and possible background double counting are not established.',
        'sources':{},'channels':{},'sum':{}}
    weight_hist=[];mass={'before':[],'after':[]}
    for ch,title,color in zip(channels,TITLES,COLORS+[ROOT.TColor.GetColor('#00AFC1')]):
        folder=BASE/'runs'/(ch+'_Roman_007');path=folder/(ch+'_Roman_007.root')
        report['sources'][str(path)]=sha(path)
        f=ROOT.TFile.Open(str(path),'READ');t=f.Get('h10');weights=[];selected=[];masses=[]
        cut=ROOT.TTreeFormula('cut_'+ch,CUT,t)
        raw_check={stage:ROOT.TH1D(ch+'_check_'+stage,'',100,0,2) for stage in ['before','after']}
        for i in range(t.GetEntries()):
            t.GetEntry(i);w=float(t.Weight);m=float(t.missmass)
            assert math.isfinite(w) and math.isfinite(m)
            weights.append(w);masses.append(m);raw_check['before'].Fill(m)
            if cut.EvalInstance():selected.append(w);raw_check['after'].Fill(m)
        assert len(weights)==1000
        assert len(selected)==t.GetEntries(CUT)
        stats=describe(weights)
        report['channels'][ch]={'before':stats,'after':describe(selected),'min_reconstructed_missing_mass_GeV':min(masses),'max_reconstructed_missing_mass_GeV':max(masses)}
        del cut;f.Close()
        tick=10.**math.floor(math.log10(max(abs(stats['min']),abs(stats['max']))))/2.
        low=math.floor(min(0.,stats['min'])/tick)*tick;high=math.ceil(max(0.,stats['max'])/tick)*tick
        if high<=stats['max']:high+=tick
        h=ROOT.TH1D(ch+'_original_Weight_all_events','',80,low,high);h.Sumw2()
        h.GetXaxis().SetNdivisions(500+round((high-low)/tick),False)
        for w in weights:h.Fill(w)
        assert h.Integral()==len(weights) and h.GetBinContent(0)==0 and h.GetBinContent(81)==0
        weight_hist.append(h)
        hist_path=folder/'plots/histograms.root';report['sources'][str(hist_path)]=sha(hist_path)
        hf=ROOT.TFile.Open(str(hist_path),'READ')
        summary=json.loads((folder/'plots/summary.json').read_text());assert summary['cut']==CUT
        for stage in ['before','after']:
            source=hf.Get('missmass_raw_'+stage)
            assert source.GetNbinsX()==100 and source.GetXaxis().GetXmin()==0 and source.GetXaxis().GetXmax()==2
            for i in range(102):close(raw_check[stage].GetBinContent(i),source.GetBinContent(i))
            hm=source.Clone(ch+'_scaled_raw_mass_'+stage);hm.Scale(scales[ch])
            for i in range(102):
                close(hm.GetBinContent(i),source.GetBinContent(i)*scales[ch])
                close(hm.GetBinError(i),source.GetBinError(i)*scales[ch])
            hm.SetLineColor(color);hm.SetLineWidth(2);hm.SetLineStyle(1);mass[stage].append(hm)
        hf.Close()
    out=ROOT.TFile(str(dest/'seven_channel_sum.root'),'RECREATE')
    labels=['#omega','#phi / 10','#rho #times 2','#eta / 22.5',"#eta' / 15",'X-phase-space #times 4',f'3#pi #times {args.three_pion_scale:g} (test)']
    for stage,hs in mass.items():
        six=hs[0].Clone('sum_six_'+stage);six.Reset('ICES')
        for h in hs[:6]:six.Add(h)
        total=six.Clone('sum_seven_'+stage);total.Add(hs[-1])
        for i in range(102):
            close(total.GetBinContent(i),math.fsum(h.GetBinContent(i) for h in hs))
            close(total.GetBinError(i)**2,math.fsum(h.GetBinError(i)**2 for h in hs))
            close(total.GetBinContent(i)-six.GetBinContent(i),hs[-1].GetBinContent(i))
        report['sum'][stage]={'six_integral':six.Integral(),'new_3pi_integral':hs[-1].Integral(),'seven_integral':total.Integral(),
                              'underflow':total.GetBinContent(0),'overflow':total.GetBinContent(101)}
        for h in [*hs,six,total]:h.Write()
        sum_canvas(dest,stage,hs,total,six,labels)
        sum_canvas(dest,stage,hs,total,six,labels,True)
    c=ROOT.TCanvas('all_seven_weights_same_setting','',2200,1150);c.Divide(4,2,.002,.002);keep=[]
    for i,(ch,h,title) in enumerate(zip(channels,weight_hist,TITLES),1):
        keep.extend(panel(c.cd(i),h,title,report['channels'][ch]['before']));h.Write()
    c.cd(8);note=ROOT.TLatex();note.SetNDC();note.SetTextFont(42);note.SetTextSize(.046)
    for y,text in [(.78,'All channels: 10.60 GeV'),(.66,'SHMS electron / HMS proton'),(.52,'Original signed Weight'),(.40,'No analysis cuts'),(.28,'No scaling or normalization')]:note.DrawLatex(.08,y,text)
    for ext in ['png','pdf']:c.SaveAs(str(dest/('all_seven_Weight_no_analysis_cuts.'+ext)))
    c.Write();c.Close()
    c=ROOT.TCanvas('new_3pi_Weight','',1050,750)
    keep=panel(c,weight_hist[-1],'X-phase-space 3#pi, 10.60 GeV',report['channels'][NEW]['before'])
    for ext in ['png','pdf']:c.SaveAs(str(dest/('three_pion_Weight_no_analysis_cuts.'+ext)))
    c.Write();c.Close()
    report['verification']='All seven inputs matched; each tree has 1000 entries; saved raw histograms match event refills bin by bin; scale factors, sums and variances checked; no source ROOT files changed.'
    assert all(sha(Path(p))==value for p,value in report['sources'].items())
    ROOT.TNamed('metadata',json.dumps(report)).Write();out.Close()
    (dest/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'new_3pi':report['channels'][NEW],'sums':report['sum'],'scales':scales,'verification':report['verification']},indent=2))


if __name__=='__main__':main()
