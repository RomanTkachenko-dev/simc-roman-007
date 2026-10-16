from pathlib import Path
import sys, json
from analyze_root import analyze
import ROOT
dest=Path(sys.argv[1]); path=dest/'phi_baseline_10p60.root'
f=ROOT.TFile.Open(str(path))
names=[key.GetName() for key in f.GetListOfKeys() if key.ReadObj().InheritsFrom('TTree')]
assert len(names)==1, names
f.Close()
result=analyze(path,names[0],dest/'plots','Henry native phi: provisional apparatus')
result['warning']='10.60 GeV is inherited from the provisional local input; not independently encoded by the 11-column event file.'
(dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
