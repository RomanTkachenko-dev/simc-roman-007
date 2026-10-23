"""Scientific diagnostic: signed model cross section and kinematic comparison.

No data points are digitized/invented from the paper. Figure-6 kinematic
settings are transcribed from its panel labels. All SIMC values are from the
original archived sample. This is not an experimental cross-section extraction.
"""
from pathlib import Path
import os, sys, json
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE / '.deps'))
os.environ['MPLCONFIGDIR'] = str(BASE / '.mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

a = np.load(BASE / 'sample_arrays.npz')
summary = json.loads((BASE / 'sample_summary.json').read_text())
c = a['pass_cut']
x = a['t'][c]
y = a['dsigma_dt_ub_GeV2'][c]
assert len(y) == 595 and (y < 0).all()
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11,
                     'axes.spines.top':False, 'axes.spines.right':False,
                     'axes.labelsize':12, 'axes.titlesize':13,
                     'savefig.facecolor':'white'})

fig, axs = plt.subplots(1,2,figsize=(12.8,6.1),gridspec_kw={'width_ratios':[1.12,1]})
fig.subplots_adjust(left=.08,right=.98,top=.80,bottom=.25,wspace=.27)
fig.suptitle('Omega cross-section check after the same acceptance cuts',
             fontsize=17,fontweight='bold',y=.965)
fig.text(.5,.89,'omega_Roman_007  |  595 accepted events  |  original signs retained',
         ha='center',fontsize=11,color='#444444')

ax=axs[0]
ax.scatter(x,y,s=15,color='#185d96',alpha=.38,edgecolors='none',label='SIMC model: accepted event values')
edges=np.linspace(.82,1.10,8)
centers=[]; means=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    use=(x>=lo)&(x<hi)
    if use.sum():
        centers.append((lo+hi)/2); means.append(y[use].mean())
ax.plot(centers,means,'o-',color='#222222',markersize=4,linewidth=1.5,label='Arithmetic mean in each t bin')
ax.axhline(0,color='#bb2b2b',linewidth=1.3)
ax.set(xlim=(.815,1.095),ylim=(-.07,.009),xlabel=r'$-t$ [GeV$^2$]',
       ylabel=r'Model $d\sigma/dt$ [$\mu$b/GeV$^2$]',title='All 595 values are negative')
ax.text(.5,.045,r'$d\sigma/dt = 2\pi\,10^6\,\mathrm{sigcm}$',transform=ax.transAxes,
        ha='center',fontsize=12,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
ax.yaxis.set_major_locator(MultipleLocator(.02))
ax.grid(alpha=.16)
ax.legend(loc='upper left',bbox_to_anchor=(0,-.20),fontsize=9,frameon=False)

ax=axs[1]
ax.scatter(a['Q2'][c],a['W'][c],s=14,color='#185d96',alpha=.4,edgecolors='none',label='SIMC after cuts')
ax.plot(1.75,2.48,'s',color='#bc4931',markersize=8)
ax.annotate('Fig. 6 left\nQ² = 1.75, W = 2.48',(1.75,2.48),xytext=(1.53,2.64),
            fontsize=10,arrowprops={'arrowstyle':'-','color':'#bc4931'})
ax.plot(2.35,2.47,'^',color='#7b4b9d',markersize=9)
ax.annotate('Fig. 6 right\nQ² = 2.35, W = 2.47',(2.35,2.47),xytext=(2.50,2.58),
            fontsize=10,arrowprops={'arrowstyle':'-','color':'#7b4b9d'})
ax.set(xlim=(1.45,4.45),ylim=(1.82,2.77),xlabel=r'$Q^2$ [GeV$^2$]',ylabel=r'$W$ [GeV]',
       title='The paper uses different kinematics')
ax.grid(alpha=.16)
ax.text(.04,.08,'SIMC arithmetic means:\nQ² = 3.75 GeV²\nW = 2.054 GeV',
        ha='left',transform=ax.transAxes,fontsize=10)
fig.text(.5,.022,'Diagnostic model values, not an acceptance-corrected measurement.  Reference: arXiv:1910.00464, Fig. 6.',
         ha='center',fontsize=9,color='#444444')
fig.savefig(BASE/'omega_cross_section_after_cuts.png',dpi=180)
fig.savefig(BASE/'omega_cross_section_after_cuts.svg')
plt.close(fig)

# A separate model-component view makes the sign failure traceable.
fig,ax=plt.subplots(figsize=(9.2,5.5))
fig.subplots_adjust(left=.13,right=.97,top=.85,bottom=.15)
ax.scatter(x,a['sigma_T_model'][c],s=12,alpha=.35,label=r'$d\sigma_T/dt$',color='#185d96')
ax.scatter(x,(a['model_epsilon']*a['sigma_L_model'])[c],s=12,alpha=.35,
           label=r'$\epsilon\,d\sigma_L/dt$',color='#b27800')
ax.scatter(x,y,s=13,alpha=.5,label=r'Total $d\sigma_T/dt+\epsilon\,d\sigma_L/dt$',color='#a82828')
ax.axhline(0,color='black',lw=1)
ax.set(xlabel=r'$-t$ [GeV$^2$]',ylabel=r'Model cross section [$\mu$b/GeV$^2$]',
       title='Active SIMC omega model: transverse and longitudinal components')
ax.legend(loc='lower left',frameon=True,fontsize=10)
ax.grid(alpha=.15)
fig.text(.5,.02,'The active formula reproduces stored sigcm to relative precision better than 10⁻⁶.\n'
         'All selected longitudinal components are negative; 231 / 595 transverse components are also negative.',
         ha='center',fontsize=9)
fig.savefig(BASE/'omega_model_components.png',dpi=180)
plt.close(fig)
print('Saved omega_cross_section_after_cuts.png / .svg and omega_model_components.png')
