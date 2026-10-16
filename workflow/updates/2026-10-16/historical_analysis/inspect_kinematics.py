"""Read-only kinematic audit of the six saved SIMC event samples."""
import json
import math
from pathlib import Path
import ROOT

base = Path(__file__).resolve().parent
channels = ['omega', 'phi', 'rho', 'eta', 'etap', 'xphsp']
cuts = json.loads((base/'pythia_vs_simc_weighted/summary.json').read_text())['cuts']
mp = .938272  # GeV, proton target mass printed in these SIMC .hist files

def stats(rows):
    return {key: {'mean': math.fsum(row[j] for row in rows)/len(rows),
                  'min': min(row[j] for row in rows), 'max': max(row[j] for row in rows)}
            for j,key in enumerate(['Q2_GeV2','W_GeV','x_B'])}

report = {'definition': 'x_B=Q2/(W^2-mp^2+Q2), proton at rest; mp=0.938272 GeV. Means are unweighted.',
          'cuts': cuts, 'channels': {}}
for ch in channels:
    path = base/'runs'/(ch+'_Roman_007')/(ch+'_Roman_007.root')
    f = ROOT.TFile.Open(str(path),'READ'); t = f.Get('h10')
    select = ROOT.TTreeFormula('cuts_'+ch,cuts,t)
    all_rows=[]; accepted=[]
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        q2,w = float(t.Q2),float(t.W)
        xb = q2/(w*w-mp*mp+q2)
        row = (q2,w,xb); all_rows.append(row)
        if select.EvalInstance():accepted.append(row)
    report['channels'][ch] = {'entries':len(all_rows),'accepted':len(accepted),
                              'before_cuts':stats(all_rows),'after_cuts':stats(accepted)}
    del select
    f.Close()
q2 = 4*10.6*6.7*math.sin(math.radians(13)/2)**2
nu = 10.6-6.7
report['central'] = {'Ebeam_GeV':10.6,'electron_p_GeV':6.7,'electron_theta_deg':13.,
                     'proton_p_GeV':1.1,'proton_theta_deg':32.,
                     'Q2_GeV2':q2,'W_GeV':math.sqrt(mp*mp+2*mp*nu-q2),'x_B':q2/(2*mp*nu),
                     'approximation':'Neglect electron mass; nominal central electron momentum and angle.'}
print(json.dumps(report,indent=2))
