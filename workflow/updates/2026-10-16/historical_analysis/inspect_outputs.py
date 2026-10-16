import ROOT, json
from pathlib import Path
ROOT.gROOT.SetBatch(True)
files = list(Path('/mnt/d/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs').glob('*/*.root'))
files += [Path('/mnt/c/Users/roman/Downloads')/x for x in ['Pythia_557ab.root','lAger_2zb.root']]
for path in files:
    f = ROOT.TFile.Open(str(path))
    print(str(path),flush=True)
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom('TTree'):
            print(obj.GetName(),obj.GetEntries(),[b.GetName() for b in obj.GetListOfBranches()],flush=True)
            if obj.GetEntries():
                obj.GetEntry(0)
                print({n:getattr(obj,n) for n in ['Weight','mx','missmass','t','Q2','hsxptar'] if obj.GetBranch(n)},flush=True)
    f.Close()
