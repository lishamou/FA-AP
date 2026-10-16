

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Patch

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
FDR = f"{PROJECT_ROOT}/results/Figure4A_PValue_FDR_Matrix.csv"
CORR_EDGES = f"{PROJECT_ROOT}/results/Figure4A_Correlation_Matrix_Clustered.csv"
OUT_DIR = f"{PROJECT_ROOT}/results"


fdr_df = pd.read_csv(FDR, index_col=0)
corr = pd.read_csv(CORR_EDGES)
modules = list(fdr_df.index)
n = len(modules)
corr_mat = np.eye(n)
for _, row in corr.iterrows():
    s1, s2, c = row['State1'], row['State2'], row['Correlation']
    if s1 in modules and s2 in modules and not pd.isna(c):
        i, j = modules.index(s1), modules.index(s2)
        corr_mat[i, j] = c
        corr_mat[j, i] = c


def get_category(name):
    if name.startswith('KEGG_'):
        return 'KEGG'
    elif name in ['CSC_Activity_PRG', 'LipidMetabolism_Risk_LMRG', 'Metabolism_General_Risk_MRG',
                  'Autophagy_Endothelial_HighRisk_AutRG', 'CytokineSignaling_HighRisk_CSIRG']:
        return 'Proprietary'
    elif name in ['Treg', 'aDCs', 'Macrophages', 'B_cells', 'Mast_cells', 'T_cell_co_stimulation',
                  'APC_co_inhibition', 'APC_co_stimulation', 'MHC_class_I', 'NK_cells',
                  'Type_I_IFN_Reponse', 'Type_II_IFN_Reponse', 'Cytolytic_activity']:
        return 'Immune'
    else:
        return None

color_map = {'KEGG': '#a6cee3', 'Proprietary': '#e31a1c', 'Immune': '#33a02c'}

def short_label(name):
    return name.replace('KEGG_', '').replace('_PATHWAY', '').replace('_SIGNALING', '')


TOP_N = 25
threshold = 0.3
sig = (fdr_df.values < 0.05) & (np.abs(corr_mat) > threshold)
np.fill_diagonal(sig, False)
degrees = sig.sum(axis=0)
top_idx = np.argsort(degrees)[::-1][:TOP_N]
top_modules = [modules[i] for i in top_idx]

G = nx.Graph()
for i, m in zip(top_idx, top_modules):
    G.add_node(m, degree=int(degrees[i]))
for i in top_idx:
    for j in top_idx:
        if i < j and sig[i, j] and abs(corr_mat[i, j]) > threshold:
            G.add_edge(modules[i], modules[j], weight=abs(corr_mat[i, j]))

print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")


LAYOUTS = {
    'spring': nx.spring_layout(G, k=1.5, iterations=200, seed=42, weight='weight'),
    'circle': nx.circular_layout(G),
    'shell': nx.shell_layout(G),
}

node_degrees = [G.degree(n) for n in G.nodes()]
max_deg = max(node_degrees) if node_degrees else 1
node_sizes = [300 + 1000 * (d / max_deg) for d in node_degrees]
node_colors = [color_map[get_category(n)] for n in G.nodes()]


for layout_name, pos in LAYOUTS.items():
    for show_kegg_label in [False, True]:
        suffix = f'v3_{layout_name}'
        suffix += '_kegg_labeled' if show_kegg_label else '_kegg_unlabeled'
        print(f"\n=== {suffix} (KEGG label: {show_kegg_label}, layout: {layout_name}) ===")

        fig, ax = plt.subplots(figsize=(13, 11))

        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color='gray', alpha=0.3,
                                width=[1 + 2 * w for w in weights])
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                               node_color=node_colors, edgecolors='black', linewidths=1.2, alpha=0.9)

        labels = {}
        for n in G.nodes():
            cat = get_category(n)
            if cat == 'Proprietary':
                labels[n] = short_label(n)
            elif cat == 'Immune':
                labels[n] = n
            elif cat == 'KEGG':
                if show_kegg_label:
                    labels[n] = short_label(n)
        if labels:
            nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7, font_weight='bold')

        legend_elements = [
            Patch(facecolor='#a6cee3', label='KEGG pathway'),
            Patch(facecolor='#e31a1c', label='Proprietary (5 R3 modules)'),
            Patch(facecolor='#33a02c', label='Immune'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=11,
                  frameon=True, facecolor='white', edgecolor='black')
        ax.set_title(f'Correlation Network (top {TOP_N} nodes, |r|>{threshold}, BH-FDR<0.05)\n'
                     f'Layout: {layout_name} | Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}', fontsize=12)
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(f"{OUT_DIR}/R3_PanelB_{suffix}.pdf", dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  → R3_PanelB_{suffix}.pdf")
