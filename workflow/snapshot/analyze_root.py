"""ROOT-only histogram production. Preserve signed weights; never infer yields."""
from pathlib import Path
import json, math, hashlib, re
import ROOT
import numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleFontSize(0.034)
ROOT.TH1.SetDefaultSumw2(True)
BASE = Path('/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN')
EXTERNAL = Path('/path/to/EIC/SIMC_RUNS/BACKGROUND_ROOT_ANALYSIS')
DOWNLOADS = Path('/path/to/data')
CUT = ('abs(hsdelta)<8.0 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035'
       ' && ssdelta>-10.0 && ssdelta<22.0 && abs(ssxptar)<0.040 && abs(ssyptar)<0.024')
SPECS = {
    'Q2':(100,0,10,'Q^{2} [GeV^{2}]'), 'W':(100,0,5,'W [GeV]'),
    't':(100,0,3,'t branch (stored value)'),
    'missmass':(100,0,2,'Reconstructed missing mass [GeV]'),
    'mx':(100,0,2,'mx branch (definition to confirm)'),
    'hsdelta':(100,-25,30,'hsdelta [%]'), 'ssdelta':(100,-25,30,'ssdelta [%]'),
    'hsxptar':(100,-.12,.12,'hsxptar [slope]'), 'ssxptar':(100,-.12,.12,'ssxptar [slope]'),
    'hsyptar':(100,-.08,.08,'hsyptar [slope]'), 'ssyptar':(100,-.08,.08,'ssyptar [slope]'),
    'hsytar':(100,-10,10,'hsytar [cm]'), 'ssytar':(100,-10,10,'ssytar [cm]'),
}
def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''): h.update(chunk)
    return h.hexdigest()
def draw_pair(a,b,destination,title):
    canvas=ROOT.TCanvas('canvas','',1050,760)
    canvas.SetLeftMargin(.15); canvas.SetBottomMargin(.12); canvas.SetTopMargin(.12)
    a.SetTitle(title); a.SetLineColor(ROOT.kBlue+1); b.SetLineColor(ROOT.kRed+1)
    a.SetLineWidth(2); b.SetLineWidth(2)
    maximum=max(a.GetMaximum(),b.GetMaximum(),0)
    minimum=min(a.GetMinimum(),b.GetMinimum(),0)
    span=maximum-minimum
    a.SetMaximum(maximum+max(span*.20,1e-30)); a.SetMinimum(minimum-span*.08)
    a.Draw('HIST'); b.Draw('HIST SAME')
    legend=ROOT.TLegend(.60,.76,.90,.87)
    legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(a,'Before analysis cuts','l'); legend.AddEntry(b,'After requested cuts','l')
    legend.Draw()
    for extension in ('png','pdf'): canvas.SaveAs(str(destination.with_suffix('.'+extension)))
    canvas.Close()
def analyze(path,tree_name,destination,label,provisional=False):
    destination.mkdir(parents=True,exist_ok=True)
    frame=ROOT.RDataFrame(tree_name,str(path))
    branches=set(str(x) for x in frame.GetColumnNames())
    selected=frame.Filter(CUT,'Requested acceptance cuts')
    count=frame.Count(); accepted=selected.Count()
    finite_w = frame.Filter('std::isfinite(Weight)') if 'Weight' in branches else None
    stats={}
    if finite_w is not None:
        stats={'negative':finite_w.Filter('Weight<0').Count(),'zero':finite_w.Filter('Weight==0').Count(),
            'nonfinite':frame.Filter('!std::isfinite(Weight)').Count(),'sum':finite_w.Sum('Weight'),
            'min':finite_w.Min('Weight'),'max':finite_w.Max('Weight')}
    # Book all histograms before triggering the RDF event loop.
    booked=[]
    for branch,(nb,lo,hi,axis) in SPECS.items():
        if branch not in branches: continue
        if branch=='t' and provisional: axis='-t [GeV^{2}] (stored branch t)'
        for mode in (['raw','weighted'] if 'Weight' in branches else ['raw']):
            ylabel='Events / bin' if mode=='raw' else 'Sum of signed Weight / bin'
            suffix=' | PROVISIONAL' if provisional else ' | no yield normalization'
            histtitle=label+suffix+';'+axis+';'+ylabel
            handles=[]
            for state,node in [('before',frame),('after',selected)]:
                node=node.Filter('std::isfinite('+branch+')')
                name=f'{branch}_{mode}_{state}'
                model=(name,histtitle,nb,lo,hi)
                if mode=='weighted':
                    handles.append(node.Filter('std::isfinite(Weight)').Histo1D(model,branch,'Weight'))
                else: handles.append(node.Histo1D(model,branch))
            booked.append((branch,mode,histtitle,handles))
    result={'file':str(path),'tree':tree_name,'before':int(count.GetValue()),'after':int(accepted.GetValue()),
        'cut':CUT,'branches':sorted(branches),'has_weight':'Weight' in branches,
        'weight_stats':{k:v.GetValue() for k,v in stats.items()},
        'missing_requested_columns':[x for x in ['hsxtar','ssxtar','missmass','t','Weight'] if x not in branches],
        'normalization':'None. Signed weights retained verbatim; raw counts are not yields.',
        'histogram_flows':{}}
    if 'Weight' in branches:
        lo=result['weight_stats']['min']; hi=result['weight_stats']['max']
        pad=(hi-lo)*.05 if hi>lo else max(abs(lo)*.05,1e-6)
        histtitle=label+' | Weight diagnostic;Weight value;Events / bin'
        hs=[node.Filter('std::isfinite(Weight)').Histo1D((f'Weight_raw_{state}',histtitle,100,lo-pad,hi+pad),'Weight')
            for state,node in [('before',frame),('after',selected)]]
        booked.append(('Weight','raw',histtitle,hs))
    output=ROOT.TFile(str(destination/'histograms.root'),'RECREATE')
    for branch,mode,title,handles in booked:
        a,b=[h.GetValue() for h in handles]
        output.cd(); a.Write(); b.Write()
        result['histogram_flows'][branch+'_'+mode]=[
            {'underflow':h.GetBinContent(0),'overflow':h.GetBinContent(h.GetNbinsX()+1),'in_range':h.Integral()} for h in (a,b)]
        draw_pair(a,b,destination/(branch+'_'+mode),title.split(';')[0])
    ROOT.TNamed('analysis_cut',CUT).Write()
    ROOT.TNamed('normalization_warning',result['normalization']).Write()
    output.Close()
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(label,result['before'],'->',result['after'],result['weight_stats'],flush=True)
    return result
def external_shapes():
    files=[ROOT.TFile.Open(str(EXTERNAL/name/'histograms.root')) for name in ['Pythia_Background','Pythia_Signal']]
    for branch in ['mx','Q2','W','t']:
        canvas=ROOT.TCanvas('shapes','',1050,760)
        canvas.SetLeftMargin(.15); canvas.SetTopMargin(.12); canvas.SetBottomMargin(.12)
        hs=[f.Get(branch+'_raw_after').Clone('shape_'+str(i)) for i,f in enumerate(files)]
        for h,color in zip(hs,[ROOT.kBlue+1,ROOT.kRed+1]):
            integral=h.Integral()
            if integral: h.Scale(1./integral)
            h.SetLineColor(color); h.SetLineWidth(2)
        hs[0].SetTitle('Supplied ROOT samples, after cuts;'+SPECS[branch][3]+';Unit-area shape / bin')
        hs[0].SetMaximum(max(h.GetMaximum() for h in hs)*1.22); hs[0].SetMinimum(0)
        hs[0].Draw('HIST'); hs[1].Draw('HIST SAME')
        legend=ROOT.TLegend(.57,.75,.89,.87); legend.SetBorderSize(0)
        for h,label in zip(hs,['Pythia Background','Signal (duplicate lAger excluded)']): legend.AddEntry(h,label,'l')
        legend.Draw()
        for ext in ['png','pdf']: canvas.SaveAs(str(EXTERNAL/(branch+'_shapes.'+ext)))
        canvas.Close()
    for f in files:f.Close()
def main():
    EXTERNAL.mkdir(exist_ok=True)
    result={}
    for channel in ['omega','phi','rho','eta','etap','xphsp']:
        name=channel+'_Roman_007'; path=BASE/'runs'/name/(name+'.root')
        if path.exists(): result[channel]=analyze(path,'h10',BASE/'runs'/name/'plots',channel+' SIMC',True)
    (BASE/'analysis_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    external={}
    for tree in ['Background','Signal']:
        external[tree]=analyze(DOWNLOADS/'Pythia_557ab.root',tree,EXTERNAL/('Pythia_'+tree),'Pythia '+tree)
    a=ROOT.RDataFrame('Signal',str(DOWNLOADS/'Pythia_557ab.root')).AsNumpy()
    b=ROOT.RDataFrame('Signal',str(DOWNLOADS/'lAger_2zb.root')).AsNumpy()
    same=set(a)==set(b) and all(np.array_equal(a[k],b[k],equal_nan=True) for k in a)
    external['duplicate_check']={'all_branches_identical':same,'branches_compared':len(a),
        'entries':len(next(iter(a.values()))),'excluded':'lAger_2zb.root:Signal' if same else None}
    assert same, 'Previously reported duplicate is no longer identical: review before combining.'
    external['file_sha256']={n:sha256(DOWNLOADS/n) for n in ['Pythia_557ab.root','lAger_2zb.root']}
    native_name='hatta-ep-phi-10GeV-smallspec.ep-phi.composite.run00001-lumi1000.simc'
    native=DOWNLOADS/native_name
    existing=Path('/path/to/EIC/SIMC_RUNS/HENRY_SIMC_FILE_INPUT/inputs')/native_name
    external['native_simc']={'provided_sha256':sha256(native),'henry_copy_sha256':sha256(existing)}
    external['native_simc']['same_file']=external['native_simc']['provided_sha256']==external['native_simc']['henry_copy_sha256']
    (EXTERNAL/'analysis_summary.json').write_text(json.dumps(external,indent=2)+'\n')
    external_shapes()
    print(json.dumps(external['duplicate_check']),flush=True)
if __name__=='__main__': main()
