import os, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
OUT_DIR = f"{PROJECT_ROOT}/results"
MAC_DIR = "../../1FIGURES/修改/Illustrator_assembly/Figure5_拼合素材/r3_subpanels"

corr = pd.read_csv(f"{OUT_DIR}/Figure4A_Correlation_Matrix_Raw.csv", index_col=0)
order = pd.read_csv(f"{OUT_DIR}/Figure4A_Clustering_Order.csv", index_col=0)['State'].tolist()
modules = [m for m in order if m in corr.index]
corr_mat = corr.loc[modules, modules].values
n = len(modules)
print(f"Loaded {n}x{n} correlation matrix")
print(f"  Range: [{corr_mat[~np.isnan(corr_mat)].min():.3f}, {corr_mat[~np.isnan(corr_mat)].max():.3f}]")

dist = 1 - corr_mat
np.fill_diagonal(dist, 0)
dist = (dist + dist.T) / 2
dist = np.clip(dist, 0, None)
dist_array = squareform(dist, checks=False)
Z = linkage(dist_array, method='average')
clustered_idx = [modules.index(m) for m in modules]
corr_clustered = corr_mat[np.ix_(clustered_idx, clustered_idx)]
modules_clustered = modules

rdbu_deep = LinearSegmentedColormap.from_list('rdbu_deep', [
    '#053061',
    '#2166ac',
    '#4393c3',
    '#92c5de',
    '#f7f7f7',
    '#f4a582',
    '#d6604d',
    '#b2182b',
    '#67001f',
], N=256)
cmap = rdbu_deep

norm = TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1.0)

fig = plt.figure(figsize=(16, 14))

title_pos = [0.05, 0.94, 0.85, 0.04]

dendro_x_pos = [0.10, 0.89, 0.78, 0.04]

heatmap_pos = [0.10, 0.10, 0.78, 0.78]

dendro_y_pos = [0.06, 0.10, 0.04, 0.78]

bottom_labels_pos = [0.10, 0.04, 0.78, 0.06]

left_labels_pos = [0.01, 0.10, 0.05, 0.78]

cbar_pos = [0.88, 0.10, 0.012, 0.78]

ax_title = fig.add_axes(title_pos)
ax_title.axis('off')
ax_title.text(0.0, 0.5, 'Spearman Correlation Matrix of Agent State Activities (TCGA)',
              ha='left', va='center', fontsize=14, fontweight='bold')

ax_dendro_x = fig.add_axes(dendro_x_pos)
dendro_x = dendrogram(Z, no_labels=True, ax=ax_dendro_x,
                      above_threshold_color='#b0b0b0', color_threshold=0)
ax_dendro_x.set_xticks([])
ax_dendro_x.set_yticks([])
for spine in ax_dendro_x.spines.values():
    spine.set_visible(False)

ax_heat = fig.add_axes(heatmap_pos)
im = ax_heat.pcolormesh(np.arange(n+1), np.arange(n+1), corr_clustered,
                          cmap=cmap, norm=norm, edgecolors='#b0b0b0', linewidth=0.3)
ax_heat.set_xlim(0, n)
ax_heat.set_ylim(n, 0)
ax_heat.set_xticks([])
ax_heat.set_yticks([])
for spine in ax_heat.spines.values():
    spine.set_visible(False)

short = [m.replace('KEGG_', '').replace('_PATHWAY', '').replace('_SIGNALING', '').replace('_CASCADE', '').replace('_BIOSYNTHESIS', '').replace('_RECOMBINATION', '').replace('_REPLICATION', '').replace('_CYCLE', '').replace('_METABOLISM', 'METAB').replace('_CYTOSKELETON', 'CYTOSK').replace('_PRESENTATION', 'PRESENT').replace('_DEGRADATION', 'DEGRAD').replace('_EXCISION_REPAIR', 'EXCIS_REP').replace('_MISMATCH_REPAIR', 'MISMATCH_REP').replace('_NUCLEOTIDE_EXCISION_REPAIR', 'NUCL_EXC_REP') for m in modules_clustered]
ax_bot = fig.add_axes(bottom_labels_pos)
ax_bot.set_xlim(0, n)
ax_bot.set_ylim(0, 1)
ax_bot.axis('off')
for i, m in enumerate(short):
    ax_bot.text(i + 0.5, 0.0, m, rotation=90, ha='right', va='top', fontsize=8)

ax_dendro_y = fig.add_axes(dendro_y_pos)
dendro_y = dendrogram(Z, no_labels=True, ax=ax_dendro_y,
                      above_threshold_color='#b0b0b0', color_threshold=0, orientation='left')
ax_dendro_y.set_xticks([])
ax_dendro_y.set_yticks([])
for spine in ax_dendro_y.spines.values():
    spine.set_visible(False)

ax_left = fig.add_axes(left_labels_pos)
ax_left.set_xlim(0, 1)
ax_left.set_ylim(n, 0)
ax_left.axis('off')
for i, m in enumerate(short):
    ax_left.text(1.0, i + 0.5, m, ha='right', va='center', fontsize=8)

ax_cbar = fig.add_axes(cbar_pos)
ax_cbar.set_xlim(0, 1)
ax_cbar.set_ylim(-0.5, 1.0)
N_rect = 100
for i in range(N_rect):
    y_bot = -0.5 + i * 1.5 / N_rect
    y_top = -0.5 + (i + 1) * 1.5 / N_rect
    y_mid = (y_bot + y_top) / 2
    color = cmap(norm(y_mid))
    ax_cbar.add_patch(plt.Rectangle((0, y_bot), 1, y_top - y_bot, facecolor=color, edgecolor='none'))
ax_cbar.set_xticks([])

ticks = [-0.5, 0, 0.5, 1.0]
ax_cbar.set_yticks(ticks)
ax_cbar.set_yticklabels([f'{t:.1f}' for t in ticks], fontsize=10)
ax_cbar.tick_params(axis='y', length=0)
ax_cbar.set_ylabel('Correlation (ρ)', fontsize=11, labelpad=4)

out_pdf = f"{OUT_DIR}/R3_Figure5A_Correlation_Heatmap.pdf"
out_png = f"{OUT_DIR}/R3_Figure5A_Correlation_Heatmap.png"
plt.savefig(out_pdf, format='pdf', pad_inches=0.2)
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"  → {out_pdf}")

if os.path.exists('/Users/'):
    os.makedirs(MAC_DIR, exist_ok=True)
    import subprocess
    subprocess.run(['scp', '-q', out_pdf, f"{MAC_DIR}/R3_PanelA_Heatmap_BH-FDR.pdf"], capture_output=True)
    subprocess.run(['scp', '-q', out_png, f"{MAC_DIR}/R3_PanelA_Heatmap_BH-FDR.png"], capture_output=True)
    print(f"  Mac: {MAC_DIR}/R3_PanelA_Heatmap_BH-FDR.pdf")
