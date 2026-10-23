"""Read this SIMC version's sequential ntuple, validate reader, plot signed model values."""
from pathlib import Path
import sys, os, struct, json, hashlib, re
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE.parent / 'omega_cross_section_check' / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_bin(path):
    with Path(path).open('rb') as f:
        def record():
            head = f.read(4)
            if not head:
                return None
            assert len(head) == 4
            n, = struct.unpack('<i', head)
            assert 0 < n < 1000000
            b, tail = f.read(n), f.read(4)
            assert len(b) == n and tail == head
            return b
        n, = struct.unpack('<i', record())
        assert 1 <= n <= 80
        names = [record().decode('ascii').strip(' \0') for _ in range(n)]
        assert len(set(names)) == n
        rows = []
        while (first := record()) is not None:
            vals = [first] + [record() for _ in range(n-1)]
            assert all(v is not None and len(v) == 8 for v in vals)
            rows.append([struct.unpack('<d', v)[0] for v in vals])
    x = np.asarray(rows)
    return {k:x[:,j] for j,k in enumerate(names)}

def cut(a):
    return ((abs(a['hsdelta']) < 8) & (abs(a['hsxptar']) < .080)
            & (abs(a['hsyptar']) < .035) & (a['ssdelta'] > -10)
            & (a['ssdelta'] < 22) & (abs(a['ssxptar']) < .040)
            & (abs(a['ssyptar']) < .024))

# Independent check against the archived official ROOT conversion, using
# our just-reproduced original-kinematics run; original ROOT used float32.
check = read_bin('//wsl.localhost/Ubuntu/opt/simc/roman007-b6025c4/worksim/setup_check_omega_007.bin')
old_path = '/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/omega_Roman_007/omega_Roman_007.root'
with uproot.open(old_path) as f:
    old = f['h10'].arrays(library='np')
verified = ['Q2','W','t','sigcm','siglab','Weight','hsdelta','hsxptar','hsyptar','ssdelta','ssxptar','ssyptar']
for k in verified:
    assert np.allclose(check[k], old[k], rtol=3e-6, atol=1e-10), k
assert cut(check).sum() == 595

a = read_bin(BASE / 'omega_W247_Q235.bin')
assert len(a['t']) == 1000
sel = cut(a)
# Keep legacy names as written; extra 'factor' header shifts trailing labels.
# No affected tail branch is used for the cross-section plot.
with uproot.recreate(BASE / 'omega_W247_Q235.root') as f:
    tree = f.mktree('h10', {k:'float64' for k in a})
    tree.extend(a)
    f['schema_note'] = 'Legacy SIMC names retained. Extra factor header shifts tail labels. W model=tprimei; epsilon model=wcmi; minus-u model=epsiloni.'

dsigma = 2*np.pi*1e6*a['sigcm']
# Source model has sigLT=sigTT=0; integrate over phi, convert MeV^2 to GeV^2.
# No abs(), event Weight, luminosity, normfac, or fitted scale is applied.
u, q2, w, eps = abs(a['epsiloni']), a['Q2i'], a['tprimei'], a['wcmi']
model = ((20.16527-4.2124*u)/q2**2 + eps*(12.2546-29.8629*u)/q2**4)/(w*w-.93827231**2)**2
assert np.allclose(model, dsigma, rtol=2e-5, atol=1e-8)
np.savez(BASE/'sample_arrays.npz', **a, pass_cut=sel, dsigma_dt_ub_GeV2=dsigma)
np.savetxt(BASE/'accepted_model_values.csv', np.column_stack([a['t'][sel],dsigma[sel],a['Q2'][sel],a['W'][sel],a['Weight'][sel]]), delimiter=',',header='minus_t_GeV2,dsigma_dt_ub_GeV2,Q2_GeV2,W_GeV,Weight',comments='')

def stats(x):
    return dict(mean=float(x.mean()), min=float(x.min()), max=float(x.max()), negative=int((x<0).sum()))
log=(BASE/'simc.log').read_text()
attempts=int(re.search(r'Last Event\s+(\d+)',log)[1])
summary = dict(source_commit='b6025c462ef33eeb620bdf4b399cd602c92e89bf',
    baseline_input='/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN/runs/omega_Roman_007/omega_Roman_007.inp',
    input_sha256=hashlib.sha256((BASE/'omega_W247_Q235.inp').read_bytes()).hexdigest(),
    nominal=dict(W_GeV=2.47,Q2_GeV2=2.35,Ebeam_GeV=10.6,electron_p_GeV=6.5656994075,electron_angle_deg=10.5432875463,proton_p_GeV=1.1,proton_angle_deg=48.6366029933),
    generation_attempts=attempts,stored_events=len(sel),after_cuts=int(sel.sum()),
    conversion='dsigma/dt = 2*pi*1e6*sigcm; signed, phi-independent model; no arbitrary scaling',
    cuts='abs(hsdelta)<8 && abs(hsxptar)<0.080 && abs(hsyptar)<0.035 && -10<ssdelta<22 && abs(ssxptar)<0.040 && abs(ssyptar)<0.024',
    additional_analysis_cuts='None on missing mass, W, Q2 or t; generator/transport constraints retained from input.',
    after_cut_stats={k:stats(a[k][sel]) for k in ['Q2','W','t','sigcm','Weight']},
    dsigma_after_cuts=stats(dsigma[sel]),epsilon_after_cuts=stats(eps[sel]),
    baseline_reader_check=dict(verified_branches=verified,after_cuts=int(cut(check).sum()),relative_tolerance=3e-6),
    model_max_abs_difference=float(abs(model-dsigma).max()))
(BASE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

t,y=a['t'][sel],dsigma[sel]
edges=np.arange(np.floor(t.min()/.02)*.02,np.ceil(t.max()/.02)*.02+.01001,.02)
centers,means=[],[]
for lo,hi in zip(edges[:-1],edges[1:]):
    mask=(t>=lo)&(t<hi)
    if mask.any():
        centers.append((lo+hi)/2); means.append(y[mask].mean())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':15,'axes.titlesize':15,'xtick.direction':'in','ytick.direction':'in'})
fig,ax=plt.subplots(figsize=(9.2,6.2))
fig.subplots_adjust(left=.15,right=.97,bottom=.14,top=.85)
ax.scatter(t,y,s=17,color='#185d96',alpha=.48,edgecolors='none',label=f'SIMC model values ({sel.sum()} accepted events)')
ax.plot(centers,means,'o-',color='black',lw=1.5,ms=4,label=r'Mean in each 0.02 GeV$^2$ bin')
ax.axhline(0,color='#777777',lw=1)
ax.set(xlabel=r'$-t$ [GeV$^2$]',ylabel=r'$d\sigma/dt$ [$\mu$b/GeV$^2$]',title=r'SIMC $\omega$ after acceptance cuts'+'\n'+r'Nominal $W=2.47$ GeV, $Q^2=2.35$ GeV$^2$')
ax.tick_params(which='both',top=True,right=True)
ax.minorticks_on(); ax.legend(loc='best',frameon=False,fontsize=10)
ax.grid(which='major',alpha=.13)
for ext in ['png','svg','pdf']:
    fig.savefig(BASE/f'simc_omega_W247_Q235_after_cuts.{ext}',dpi=180)
plt.close(fig)
print(json.dumps(summary,indent=2))
