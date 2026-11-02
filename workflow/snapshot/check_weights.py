"""Read-only sign audit of existing outputs; does not alter physics models."""
from pathlib import Path
import ROOT, numpy as np, json
base=Path('/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN')
report={}
for channel in ['omega','phi','rho','eta','etap','xphsp']:
    name=channel+'_Roman_007'
    a=ROOT.RDataFrame('h10',str(base/'runs'/name/(name+'.root'))).AsNumpy(['Weight','siglab','sigcm','missmass'])
    report[channel]={key:{'negative':int(np.sum(value<0)),'nonfinite':int(np.sum(~np.isfinite(value))),
        'min':float(value.min()),'max':float(value.max())} for key,value in a.items()}
    report[channel]['weight_sign_matches_siglab']=bool(np.array_equal(np.sign(a['Weight']),np.sign(a['siglab'])))
    report[channel]['siglab_sign_matches_sigcm']=bool(np.array_equal(np.sign(a['siglab']),np.sign(a['sigcm'])))
(base/'weight_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
