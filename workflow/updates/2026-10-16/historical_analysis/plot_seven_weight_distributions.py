"""Seven separate Weight diagnostics; explicitly retain the different 3pi setting."""
import json
import math
from pathlib import Path
import ROOT
from plot_weight_distributions import BASE, CHANNELS, TITLES, describe, panel, sha

OUT=BASE/'weight_distributions_seven_channels'
SMOKE=BASE/'phi007_integration/smoke_20260908_075908/xphase_threebody_smoke.root'


def main():
    OUT.mkdir(exist_ok=True)
    prior=json.loads((BASE/'weight_distributions_after_cuts/summary.json').read_text())
    smoke_prior=json.loads((SMOKE.parent/'root_validation.json').read_text())
    channels=CHANNELS+['xphsp_3pi']
    titles=TITLES[:-1]+['X-phase-space (legacy)','X-phase-space 3#pi']
    report={'expression':'Weight','histogram_fill_weight':1,'analysis_selection':'All stored h10 events; no additional analysis cuts.',
            'stage':'After SIMC. Counts are stored events, not generation attempts.',
            'scaling':'No channel scaling or normalization; signs preserved. Each panel has its own x range and power of ten.',
            'settings_warning':'Original six samples: Ebeam=10.60 GeV, SHMS/HMS. Three-pion smoke: Ebeam=4.70897 GeV, SOS/HMS. Separate diagnostics, not a combined spectrum or comparable absolute yields.',
            'six_channel_after_cuts_reference':str(BASE/'weight_distributions_after_cuts/all_six_Weight_after_cuts.png'),
            'channels':{},'source_sha256':{}}
    hs=[]
    for channel in channels:
        path=SMOKE if channel=='xphsp_3pi' else BASE/'runs'/(channel+'_Roman_007')/(channel+'_Roman_007.root')
        report['source_sha256'][str(path)]=sha(path)
        f=ROOT.TFile.Open(str(path),'READ');assert f and not f.IsZombie()
        t=f.Get('h10');assert t and t.GetBranch('Weight')
        values=[]
        for i in range(t.GetEntries()):
            t.GetEntry(i);w=float(t.Weight);assert math.isfinite(w);values.append(w)
        assert len(values)==t.GetEntries()
        counts=describe(values)
        if channel=='xphsp_3pi':
            assert counts['entries']==smoke_prior['entries']
            assert counts['negative']==smoke_prior['negative_weight']
        else:
            expected=prior['channels'][channel]['before_cuts']
            for key in ['entries','negative','positive','zero']:assert counts[key]==expected[key]
            assert math.isclose(counts['sum_weights'],expected['sum_weights'],rel_tol=1e-12)
        # Align limits/ticks to readable half-decade steps, including zero.
        tick=10.**math.floor(math.log10(max(abs(counts['min']),abs(counts['max']))))/2.
        low=math.floor(min(0.,counts['min'])/tick)*tick
        high=math.ceil(max(0.,counts['max'])/tick)*tick
        if high<=counts['max']:high+=tick
        divisions=round((high-low)/tick)
        h=ROOT.TH1D(channel+'_original_Weight_all_events','',80,low,high);h.Sumw2()
        h.GetXaxis().SetNdivisions(500+divisions,False)
        for w in values:h.Fill(w)
        assert h.Integral()==len(values) and h.GetEntries()==len(values)
        assert h.GetBinContent(0)==0 and h.GetBinContent(81)==0
        assert math.isclose(sum(h.GetBinError(i)**2 for i in range(82)),len(values),rel_tol=1e-12)
        report['channels'][channel]={'file':str(path),'tree':'h10','statistics':counts,
                                     'Ebeam_GeV':4.70897 if channel=='xphsp_3pi' else 10.60,
                                     'bins':80,'range':[low,high]}
        hs.append(h);f.Close()
    output=ROOT.TFile(str(OUT/'seven_weight_distributions.root'),'RECREATE')
    # Same-setting six-channel sheet; exclude the different-energy 3pi smoke test.
    six=ROOT.TCanvas('six_original_weight_distributions','',1800,1100);six.Divide(3,2,.002,.002)
    six_keep=[]
    for i,(channel,h,title) in enumerate(zip(channels[:6],hs[:6],titles[:6]),1):
        six_keep.extend(panel(six.cd(i),h,title,report['channels'][channel]['statistics']))
    for ext in ['png','pdf']:six.SaveAs(str(OUT/('all_six_Weight_before_analysis_cuts.'+ext)))
    six.Write();six.Close()
    canvas=ROOT.TCanvas('seven_original_weight_distributions','',2200,1150);canvas.Divide(4,2,.002,.002)
    keep=[]
    for i,(channel,h,title) in enumerate(zip(channels,hs,titles),1):
        keep.extend(panel(canvas.cd(i),h,title,report['channels'][channel]['statistics']))
    canvas.cd(8)
    note=ROOT.TLatex();note.SetNDC();note.SetTextFont(42);note.SetTextSize(.046)
    for y,text in [( .78,'Original signed Weight'),(.67,'All stored h10 events'),(.59,'No additional analysis cuts'),
                   (.45,'Six channels: 10.60 GeV'),(.37,'3#pi smoke: 4.70897 GeV'),(.23,'Separate settings; no sum')]:
        note.DrawLatex(.08,y,text)
    keep.append(note)
    for ext in ['png','pdf']:canvas.SaveAs(str(OUT/('all_seven_Weight_before_analysis_cuts.'+ext)))
    canvas.Write();canvas.Close()
    for channel,h,title in zip(channels,hs,titles):
        canvas=ROOT.TCanvas('canvas_'+channel+'_all_Weight','',1050,750)
        keep=panel(canvas,h,title,report['channels'][channel]['statistics'])
        for ext in ['png','pdf']:canvas.SaveAs(str(OUT/(channel+'_Weight_before_analysis_cuts.'+ext)))
        h.Write();canvas.Write();canvas.Close()
    assert all(sha(Path(path))==digest for path,digest in report['source_sha256'].items())
    report['verification']='Seven original ROOT trees read; counts/signs independently cross-checked; unit-count histogram integrals and Sumw2 verified; no flows or nonfinite weights; sources unchanged.'
    ROOT.TNamed('metadata',json.dumps(report)).Write();output.Close()
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({ch:obj['statistics'] for ch,obj in report['channels'].items()},indent=2))


if __name__=='__main__':main()
