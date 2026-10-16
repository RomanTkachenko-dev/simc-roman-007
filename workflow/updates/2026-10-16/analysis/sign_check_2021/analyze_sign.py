"""Check raw ntuple signs for the isolated 2021 master build."""
from pathlib import Path
import sys, struct, json, re, tarfile
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent/'omega_cross_section_check'/'.deps'))
import numpy as np
import uproot

def read_bin(path):
    with path.open('rb') as f:
        def record():
            head=f.read(4)
            if not head:return None
            assert len(head)==4
            n,=struct.unpack('<i',head)
            assert 0<n<1000000
            data=f.read(n)
            assert len(data)==n and f.read(4)==head
            return data
        n,=struct.unpack('<i',record())
        names=[record().decode('ascii').strip(' \0') for _ in range(n)]
        assert len(set(names))==n
        rows=[]
        while (first:=record()) is not None:
            row=[first]+[record() for _ in range(n-1)]
            assert all(x is not None and len(x)==8 for x in row)
            rows.append([struct.unpack('<d',x)[0] for x in row])
    values=np.array(rows)
    return {k:values[:,i] for i,k in enumerate(names)}

def stats(x):
    return dict(negative=int((x<0).sum()),zero=int((x==0).sum()),positive=int((x>0).sum()),nonfinite=int((~np.isfinite(x)).sum()),min=float(np.min(x)),max=float(np.max(x)),mean=float(np.mean(x)))

report=dict(commit='8b43f9d0dd3be977df93b35abc528a389ca70d5c',commit_date='2021-06-24',branch='master',source='billlee77/simc_gfortran',tests={})
for name,keys,reaction in [('test_eep_h_root',['sigcc','Weight'],"H(e,e'p) elastic"),('test_pion_h_sign',['sigcm','siglab','Weight'],"H(e,e'pi+)n")]:
    a=read_bin(BASE/f'{name}.bin')
    assert len(a['Weight'])==1000
    with uproot.recreate(BASE/f'{name}.root') as f:
        t=f.mktree('h10',{k:'float64' for k in a})
        t.extend(a)
    with uproot.open(BASE/f'{name}.root') as f:
        for k in keys:
            assert np.array_equal(a[k],f['h10'][k].array(library='np'))
    log=(BASE/f'{name}.log').read_text()
    report['tests'][name]=dict(reaction=reaction,stored_events=len(a['Weight']),generation_attempts=int(re.search(r'Last Event\s+(\d+)',log)[1]),selection='All stored events; no additional analysis cuts',signs={k:stats(a[k]) for k in keys})

source=Path('//wsl.localhost/Ubuntu/opt/simc/bill-master-8b43f9d')
checked=0
with tarfile.open(BASE/'source_8b43f9d.tar') as tar:
    for m in tar:
        if m.isfile() and m.name.endswith('.f'):
            assert tar.extractfile(m).read()==(source/m.name).read_bytes(),m.name
            checked+=1
report['unchanged_fortran_source_files_verified']=checked
report['build_compatibility']='Modern compiler flags; trgDeriv compiled separately with byte-identical source; canonical trg_track.f unchanged.'
report['input_changes']='eep: existing test input with field spelling converted. pion: upstream test_pion_h; . field spelling converted to %, closing label typo corrected, ngen reduced from 10000 to 1000; physical settings retained. See .diff files.'
report['scope']='master has no exclusive-omega model. These are elastic and pion sign checks, not an omega model comparison. Positive signs do not by themselves identify or confirm the specific reported bug.'
(BASE/'sign_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
