

import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt


CORR_FILE = "/data/share/mouls/LIHC_FMAP_R3/results/Figure4A_Correlation_Matrix_Clustered.csv"
corr = pd.read_csv(CORR_FILE)
modules = sorted(set(corr['State1']).union(set(corr['State2'])))
n = len(modules)
corr_mat = pd.DataFrame(np.eye(n), index=modules, columns=modules)


for _, row in corr.iterrows():
    s1, s2, c = row['State1'], row['State2'], row['Correlation']
    if s1 in modules and s2 in modules and not pd.isna(c):
        corr_mat.loc[s1, s2] = c
        corr_mat.loc[s2, s1] = c

print(f"Loaded correlation matrix: {corr_mat.shape}, modules={n}")


HUB_MODULES = [
    'KEGG_CELL_CYCLE',
    'Treg',
    'APC_co_stimulation',
    'KEGG_MAPK_SIGNALING_PATHWAY',
    'KEGG_FATTY_ACID_METABOLISM',
]


for h in HUB_MODULES:
    if h not in modules:
        print(f"  WARNING: {h} not in modules!")


SHORT_LABEL = {
    'KEGG_CELL_CYCLE': 'KEGG_CELL_CYCLE',
    'Treg': 'Treg',
    'APC_co_stimulation': 'APC_co_stimulation',
    'KEGG_MAPK_SIGNALING_PATHWAY': 'KEGG_MAPK',
    'KEGG_FATTY_ACID_METABOLISM': 'KEGG_FATTY_ACID',
}


fig, axes = plt.subplots(5, 1, figsize=(14*0.75, 18*4/5), sharex=False)


POS_COLOR = '#2166ac'
NEG_COLOR = '#b2182b'
NS_COLOR = '#cccccc'

for ax_idx, hub in enumerate(HUB_MODULES):
    ax = axes[ax_idx]


    corrs_with_hub = corr_mat.loc[hub].drop(hub).dropna()


    corrs_sorted = corrs_with_hub.sort_values()


    colors = []
    for r in corrs_sorted.values:
        if abs(r) > 0.3:
            colors.append(POS_COLOR if r > 0 else NEG_COLOR)
        else:
            colors.append(NS_COLOR)


    x = np.arange(len(corrs_sorted))
    ax.bar(x, corrs_sorted.values, color=colors, edgecolor='white', linewidth=0.3)


    ax.axhline(0.3, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axhline(-0.3, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axhline(0, color='black', linewidth=0.6)


    ax.set_title(f'({chr(65+ax_idx+2)}) {SHORT_LABEL.get(hub, hub)} (hub) vs all other modules',
                 fontsize=12, fontweight='bold', loc='left')
    ax.set_ylabel('Spearman correlation (ρ)', fontsize=10)
    ax.set_ylim(-1.0, 1.0)


    ax.set_xticks([])

    ax.grid(axis='y', alpha=0.3)


fig.suptitle('Figure 5C. Spearman correlation of 5 hub modules with all other functional modules\n'
             '(BH-FDR corrected; bars colored by |r|>0.3 significance: blue=positive, red=negative, gray=NS)',
             fontsize=13, y=1.005)

plt.tight_layout()


OUT_PDF = "/data/share/mouls/LIHC_FMAP_R3/results/R3_PanelC_BarPlot_Hubs.pdf"
plt.savefig(OUT_PDF, dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  → {OUT_PDF}")


OUT_PNG = OUT_PDF.replace('.pdf', '.png')
fig, axes = plt.subplots(5, 1, figsize=(14*0.75, 18*4/5), sharex=False)
for ax_idx, hub in enumerate(HUB_MODULES):
    ax = axes[ax_idx]
    corrs_with_hub = corr_mat.loc[hub].drop(hub).dropna()
    corrs_sorted = corrs_with_hub.sort_values()
    colors = [POS_COLOR if abs(r) > 0.3 and r > 0 else (NEG_COLOR if abs(r) > 0.3 and r < 0 else NS_COLOR) for r in corrs_sorted.values]
    x = np.arange(len(corrs_sorted))
    ax.bar(x, corrs_sorted.values, color=colors, edgecolor='white', linewidth=0.3)
    ax.axhline(0.3, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axhline(-0.3, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axhline(0, color='black', linewidth=0.6)
    ax.set_title(f'({chr(65+ax_idx+2)}) {SHORT_LABEL.get(hub, hub)} (hub) vs all other modules', fontsize=12, fontweight='bold', loc='left')
    ax.set_ylabel('Spearman correlation (ρ)', fontsize=10)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xticks([])
    ax.grid(axis='y', alpha=0.3)
fig.suptitle('Figure 5C. Spearman correlation of 5 hub modules with all other functional modules', fontsize=13, y=1.005)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
plt.close()
print(f"  → {OUT_PNG}")


import subprocess
MAC_PATH = "../../1FIGURES/Figure_supplementary/投稿/Illustrator_assembly/Figure5_拼合素材/r3_subpanels/R3_PanelC_BarPlot_Hubs.pdf"
r = subprocess.run(['scp', '-q', OUT_PDF, MAC_PATH], capture_output=True, text=True)
print(f"\nscp to mac: rc={r.returncode}")


import fitz
doc = fitz.open(MAC_PATH)
fonts = set()
for page in doc:
    for f in page.get_fonts():
        fonts.add((f[3], f[2]))
print(f"Fonts in PDF: {fonts}")
