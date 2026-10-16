"""Compare omega Q2, W and Bjorken x before/after event weighting."""
from pathlib import Path
import hashlib
import json
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.TH1.AddDirectory(False)
BASE=Path(__file__).resolve().parent
OUT=BASE/'omega_kinematics_weight_comparison'
SPECS=[('Q2','Q^{2} [GeV^{2}]',60,2.5,4.5),('W','W [GeV]',60,1.7,2.7),('xB','Bjorken x',60,.25,.70)]
MP=.938272


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def style(pad,h,title,ytitle):
    pad.cd();pad.SetLeftMargin(.15);pad.SetRightMargin(.055);pad.SetBottomMargin(.17);pad.SetTopMargin(.12);pad.SetTicks(1,1)
    h.SetTitle('');h.GetXaxis().SetTitle(title);h.GetYaxis().SetTitle(ytitle)
    for ax in [h.GetXaxis(),h.GetYaxis()]:ax.SetTitleSize(.047);ax.SetLabelSize(.040)
    h.GetYaxis().SetTitleOffset(1.5);h.GetXaxis().SetTitleOffset(1.2)
    h.SetMinimum(0)


def comparison(pad,raw,weighted,title,ytitle):
    a=raw.Clone(raw.GetName()+'_shape');b=weighted.Clone(weighted.GetName()+'_shape')
    a.Scale(1./a.Integral('width'));b.Scale(1./b.Integral('width'))
    for h in [a,b]:assert math.isclose(h.Integral('width'),1.,rel_tol=1e-12)
    style(pad,a,title,ytitle)
    a.SetMaximum(max(a.GetMaximum(),b.GetMaximum())*1.30)
    a.SetLineColor(ROOT.kBlue+1);a.SetLineWidth(2);b.SetLineColor(ROOT.kRed);b.SetLineWidth(2)
    a.Draw('HIST');b.Draw('HIST SAME')
    legend=ROOT.TLegend(.61,.72,.925,.87);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.036)
    legend.AddEntry(a,'No Weight','l');legend.AddEntry(b,'With Weight','l');legend.Draw();pad.RedrawAxis()
    return [a,b,legend]


def main():
    OUT.mkdir(exist_ok=True)
    path=BASE/'runs/omega_Roman_007/omega_Roman_007.root'
    input_path=path.with_suffix('.inp')
    cuts=json.loads((BASE/'pythia_vs_simc_weighted/summary.json').read_text())['cuts']
    digest=sha(path)
    f=ROOT.TFile.Open(str(path),'READ');assert f and not f.IsZombie()
    t=f.Get('h10');select=ROOT.TTreeFormula('same_acceptance',cuts,t)
    branches=[br.GetName() for br in t.GetListOfBranches()]
    assert all(br in branches for br in ['Q2','W','Weight'])
    hist={}
    for stage in ['before_cuts','after_cuts']:
        for key,title,n,low,high in SPECS:
            for mode in ['raw','weighted']:
                h=ROOT.TH1D(f'{key}_{mode}_{stage}','',n,low,high);h.Sumw2();hist[(stage,key,mode)]=h
    rows=[];after=[]
    for i in range(t.GetEntries()):
        t.GetEntry(i);q2,w,weight=float(t.Q2),float(t.W),float(t.Weight)
        xb=q2/(w*w-MP*MP+q2)
        assert all(math.isfinite(v) for v in [q2,w,weight,xb])
        row=(q2,w,xb,weight);rows.append(row)
        passed=bool(select.EvalInstance())
        if passed:after.append(row)
        for stage in ['before_cuts']+(['after_cuts'] if passed else []):
            for (key,*_),value in zip(SPECS,row):
                hist[(stage,key,'raw')].Fill(value)
                hist[(stage,key,'weighted')].Fill(value,weight)
    assert len(rows)==1000 and len(after)==595
    report={'source_file':str(path),'source_sha256':digest,'input_sha256':sha(input_path),
        'sample_source_commit':'b6025c462ef33eeb620bdf4b399cd602c92e89bf',
        'sample_note':'Existing 10.60 GeV omega run, not the new 4.70897 GeV three-pion smoke sample. Bill phi-007 does not modify the omega production source.',
        'x_definition':'x_B=Q2/(W^2-mp^2+Q2), mp=0.938272 GeV for a stationary proton; derived event by event from reconstructed Q2 and W, not a stored x branch.',
        'cuts':cuts,'comparison':'Same accepted events with fill weight 1 versus original signed Weight. No channel factors or Pythia normalization.',
        'shape_normalization':'For shape overlays only: each raw/weighted histogram separately has integral(density dX)=1 over its displayed range. The absolute comparison also preserves counts and sum Weight.',
        'stage_statistics':{}}
    for stage,rs in [('before_cuts',rows),('after_cuts',after)]:
        sw=math.fsum(row[3] for row in rs);sw2=math.fsum(row[3]**2 for row in rs)
        assert sw>0 and all(row[3]>0 for row in rs)
        stats={'entries':len(rs),'negative_weights':sum(row[3]<0 for row in rs),'sum_weight':sw,'sum_weight_squared':sw2,'variables':{}}
        for j,(key,title,n,low,high) in enumerate(SPECS):
            a=hist[(stage,key,'raw')];b=hist[(stage,key,'weighted')]
            assert a.Integral()==len(rs)
            assert math.isclose(b.Integral(),sw,rel_tol=1e-12)
            assert math.isclose(sum(b.GetBinError(i)**2 for i in range(n+2)),sw2,rel_tol=1e-12)
            assert all(h.GetBinContent(0)==0 and h.GetBinContent(n+1)==0 for h in [a,b])
            stats['variables'][key]={'unweighted_mean':math.fsum(row[j] for row in rs)/len(rs),
                'weighted_mean':math.fsum(row[j]*row[3] for row in rs)/sw,
                'min':min(row[j] for row in rs),'max':max(row[j] for row in rs),'bins':n,'range':[low,high]}
        report['stage_statistics'][stage]=stats
    del select;f.Close()
    out=ROOT.TFile(str(OUT/'omega_kinematics.root'),'RECREATE')
    for h in hist.values():h.Write()
    units=['Unit-area density [GeV^{-2}]','Unit-area density [GeV^{-1}]','Unit-area density']
    c=ROOT.TCanvas('omega_weight_shape_comparison','',1800,650);c.Divide(3,1,.002,.002);keep=[]
    for i,((key,title,*_),unit) in enumerate(zip(SPECS,units),1):
        keep.extend(comparison(c.cd(i),hist[('after_cuts',key,'raw')],hist[('after_cuts',key,'weighted')],title,unit))
    for ext in ['png','pdf']:c.SaveAs(str(OUT/('omega_Q2_W_x_weight_shapes_after_cuts.'+ext)))
    c.Write()
    for obj in keep:
        if obj.InheritsFrom('TH1'):obj.Write()
    c.Close()
    c=ROOT.TCanvas('omega_absolute_weight_comparison','',1800,1100);c.Divide(3,2,.002,.002);keep=[]
    for row,mode in enumerate(['raw','weighted']):
        for col,(key,title,*_) in enumerate(SPECS):
            h=hist[('after_cuts',key,mode)].Clone(f'absolute_display_{key}_{mode}')
            pad=c.cd(row*3+col+1)
            style(pad,h,title,'Events / bin' if mode=='raw' else 'Sum of Weight / bin')
            # Display signed sums in units of 1e-6 to avoid clipped ROOT exponents.
            if mode=='weighted':h.Scale(1e6);h.GetYaxis().SetTitle('Sum of Weight [#times10^{-6}] / bin')
            h.SetLineColor(ROOT.kBlue+1 if mode=='raw' else ROOT.kRed);h.SetLineWidth(2);h.SetMaximum(h.GetMaximum()*1.2);h.Draw('HIST')
            label=ROOT.TLatex();label.SetNDC();label.SetTextFont(42);label.SetTextSize(.046)
            label.DrawLatex(.17,.93,'No Weight' if mode=='raw' else 'With Weight')
            keep.extend([h,label])
    for ext in ['png','pdf']:c.SaveAs(str(OUT/('omega_Q2_W_x_absolute_after_cuts.'+ext)))
    c.Write();c.Close()
    for (key,title,*_),unit in zip(SPECS,units):
        c=ROOT.TCanvas('omega_'+key+'_comparison','',1000,750)
        keep=comparison(c,hist[('after_cuts',key,'raw')],hist[('after_cuts',key,'weighted')],title,unit)
        for ext in ['png','pdf']:c.SaveAs(str(OUT/(f'omega_{key}_weight_comparison.'+ext)))
        c.Write();c.Close()
    assert sha(path)==digest
    report['verification']='595 identical accepted events; positive original omega weights; zero histogram flows; exact count/sumw/Sumw2 checks; normalized areas one; source unchanged.'
    ROOT.TNamed('metadata',json.dumps(report)).Write();out.Close()
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['stage_statistics']['after_cuts'],indent=2))
    # Read-only smoke output check, kept separate from omega analysis.
    smoke_dir=BASE/'phi007_integration'
    for p in sorted(smoke_dir.glob('smoke_*/xphase_threebody_smoke.root')):
        f=ROOT.TFile.Open(str(p),'READ');t=f.Get('h10');weights=[]
        for i in range(t.GetEntries()):
            t.GetEntry(i);weights.append(float(t.Weight))
        info={'file':str(p),'entries':len(weights),'negative_weight':sum(w<0 for w in weights),
              'nonfinite_weight':sum(not math.isfinite(w) for w in weights),'sum_weight':math.fsum(weights),
              'min_weight':min(weights) if weights else None,'max_weight':max(weights) if weights else None}
        f.Close();assert info['entries']>0 and info['nonfinite_weight']==0
        (p.parent/'root_validation.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info))


if __name__=='__main__':main()
