"""Publication figures generated from the audited FM-AP result tables."""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='Generate FM-AP figures from result tables.')
parser.add_argument('--results',type=Path,default=Path('results_final'))
parser.add_argument('--icgc-summary',type=Path,default=Path('icgc_manifest.json'))
parser.add_argument('--out',type=Path,default=Path('R4_package/Figures'))
args=parser.parse_args()
r=args.results
s=json.loads((r/'summary.json').read_text())
ext=json.loads(args.icgc_summary.read_text())
f=args.out
f.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
colors=['#1864a1','#c44d37','#007763','#80529c']
def save(fig,n):
    for suffix in ['pdf','png']:
        fig.savefig(f/f'Figure_{n}.{suffix}',dpi=320,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def panel(ax,text):
    ax.text(-.05,1.10,text,transform=ax.transAxes,fontweight='bold',fontsize=13,va='top')

fig,ax=plt.subplots(figsize=(10.4,4.4))
ax.axis('off')
def box(x,y,t,accent=False):
    ax.text(x,y,t,ha='center',va='center',fontsize=10.3,transform=ax.transAxes,
            bbox={'boxstyle':'round,pad=.6','facecolor':'#e8f3fa' if accent else '#f6f7f8','edgecolor':'#517b96'})
def arrow(x,y,xx,yy):
    ax.annotate('',xy=(xx,yy),xytext=(x,y),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':1.4,'color':'#496372'})
box(.16,.83,'TCGA-LIHC\n371 primary HCC profiles',True)
box(.16,.52,'Rank expression within each sample\n63 gene sets',True)
box(.16,.19,f"Mutation analysis\n{s['mutation']['n']} complete cases")
box(.49,.19,f"Internal survival evaluation\n{s['n_after_os_and_clinical']} patients, {s['n_events']} deaths")
box(.81,.83,'ICGC LIRI-JP legacy expression\n445 specimens, 232 donors',True)
box(.81,.52,'243 tumor-tagged specimens\n223 donors with one tumor-tagged specimen')
box(.81,.19,f"HCCDB18 and follow-up match\n{ext['n_analysis']} patients, {ext['events']} deaths",True)
arrow(.16,.735,.16,.63);arrow(.16,.401,.16,.29);arrow(.25,.41,.423,.287)
arrow(.81,.73,.81,.632);arrow(.81,.404,.81,.295)
save(fig,1)

mut=pd.read_csv(r/'mutation_matched_purity_rank.csv')
clinical=pd.read_csv(r/'mutation_matched_clinical_rank.csv')
purity=s['mutation_matched_purity_rank']
examples=['KEGG_CELL_CYCLE','KEGG_HOMOLOGOUS_RECOMBINATION','KEGG_DNA_REPLICATION','Type_II_IFN_Reponse',
          'KEGG_PPAR_SIGNALING_PATHWAY','KEGG_FATTY_ACID_METABOLISM','CCR','LipidMetabolism_Risk_LMRG']
labels=['Cell cycle','Homologous recombination','DNA replication','Type II IFN response','PPAR signaling',
        'Fatty acid metabolism','Chemokine receptors','Lipid metabolism signature']
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.6,5.2),gridspec_kw={'width_ratios':[1.25,1.]},layout='constrained')
for j,(gene,color) in enumerate(zip(['TP53','CTNNB1'],colors)):
    for k,(table,label,face,marker) in enumerate([
        (clinical,'clinical','white','o'),(mut,'+ purity',color,'s')
    ]):
        z=table.loc[table.mutation.eq(gene)].set_index('module').loc[examples]
        yy=np.arange(len(examples))+.19*(2*j-1)+.08*(2*k-1)
        ax.errorbar(z.effect,yy,xerr=[z.effect-z.low,z.high-z.effect],fmt=marker,ms=4.2,
                    mfc=face,color=color,alpha=.60 if k==0 else 1.,label=f'{gene} {label}',
                    capsize=2,elinewidth=1.0)
ax.axvline(0,c='#7c858a',lw=1,ls='--')
ax.set(yticks=np.arange(len(examples)),yticklabels=labels,xlabel='Adjusted score difference (rank points)')
ax.invert_yaxis()
ax.legend(frameon=False,loc='lower center',bbox_to_anchor=(.53,1.005),ncol=2,fontsize=8.5)
panel(ax,'A')
pivot=mut.pivot(index='module',columns='mutation',values='effect')
families=mut.set_index('module').family.drop_duplicates()
for j,(family,z) in enumerate(mut[mut.mutation.eq('TP53')].groupby('family',sort=False)):
    k=pivot.loc[z.module]
    bx.scatter(k.TP53,k.CTNNB1,s=27,color=colors[j],alpha=.83,label=family)
bx.axhline(0,c='#a0a8ad',lw=.8);bx.axvline(0,c='#a0a8ad',lw=.8)
bx.set(xlabel='TP53 adjusted score difference',ylabel='CTNNB1 adjusted score difference')
bx.legend(frameon=False,fontsize=8,loc='lower left')
bx.set_title(f"After purity adjustment · r = {purity['r']:.2f} · P = {purity['perm_p']:.3f}",fontsize=10.5,pad=10)
panel(bx,'B')
save(fig,2)

net=pd.read_csv(r/'network_all_pairs.csv')
select=net[(net.q<.05)&(net.rho.abs()>=.3)&net.overlap.gt(0)]
fig,ax=plt.subplots(figsize=(5.7,5.05),layout='constrained')
ax.plot([-1,1],[-1,1],'--',lw=1,color='#717d84')
ax.scatter(select.rho,select.rho_without_shared_genes,c=select.overlap,cmap='viridis',s=20,alpha=.85)
ax.set(xlim=(-1,1),ylim=(-1,1),xlabel='Original Spearman ρ',ylabel='Spearman ρ after shared genes are removed')
ax.set_aspect('equal',adjustable='box')
fig.colorbar(ax.collections[0],ax=ax,shrink=.78,label='Shared genes')
save(fig,3)

fig,axs=plt.subplots(2,2,figsize=(10.8,7.6),layout='constrained')
for j,(prefix,title) in enumerate([('internal','TCGA-LIHC'),('external','ICGC LIRI-JP')]):
    a=axs[0,j]; b=axs[1,j]
    met=pd.read_csv(r/f'{prefix}_metrics.csv')
    cal=pd.read_csv(r/f'{prefix}_calibration_plot.csv')
    met=met.loc[met.model.isin(['Clinical','FM-AP','Integrated']) & met.metric.isin(['C_index','AUC_3y'])]
    for k,name in enumerate(['Clinical','FM-AP','Integrated']):
        z=met.loc[met.model.eq(name)]
        pos=np.arange(2)+.19*(1-k)
        a.errorbar(z.estimate,pos,xerr=[z.estimate-z.low,z.high-z.estimate],fmt='o',ms=5,color=colors[k],capsize=2,label=name)
        m=cal.loc[cal.model.eq(name)]
        b.plot(m.predicted,m.observed,'o-',ms=4,color=colors[k],label=name)
    a.axvline(.5,ls='--',lw=1,c='#898e91')
    a.set(xlim=(.3,.95),ylim=(-.45,1.45),yticks=[0,1],yticklabels=['C-index','3-year AUC'],xlabel='Estimate and 95% CI')
    a.set_title(title,fontsize=11.5)
    b.plot([0,1],[0,1],'--',lw=1,c='#898e91')
    b.set(xlim=(0,.75),ylim=(0,.75),xlabel='Mean predicted 3-year mortality',ylabel='Kaplan–Meier 3-year mortality')
    b.set_aspect('equal',adjustable='box')
    panel(a,'AB'[j]);panel(b,'CD'[j])
axs[0,0].legend(frameon=False,fontsize=8,loc='best')
save(fig,4)
