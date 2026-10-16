

import os, sys, pickle, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')

PROJECT_ROOT = '/data/share/mouls/LIHC_FMAP_R3'
DATA_PKL = f'{PROJECT_ROOT}/data/survival_prediction_data.pkl'
MODULE_LIST = f'{PROJECT_ROOT}/data/module_list.csv'
SIGNATURE = f'{PROJECT_ROOT}/data/signature.csv'
TCGA_EXPR = f'{PROJECT_ROOT}/data/tcga_expr.tsv'
OUT_DIR = f'{PROJECT_ROOT}/results'


print('=== 1. Loading GSVA scores ===')
with open(DATA_PKL, 'rb') as f:
    model_data = pickle.load(f)
if model_data['fustat'].dtype == 'O':
    model_data['fustat'] = model_data['fustat'].map({'Alive': 0, 'Dead': 1})

module_map = pd.read_csv(MODULE_LIST).set_index('Feature')['Module'].to_dict()
feature_cols = [c for c in model_data.columns if c.startswith('Feature_')]
rename_map = {c: module_map.get(c, c) for c in feature_cols}
patient_ids = model_data['PatientID_short'].astype(str).tolist()
gsva_df = model_data[feature_cols].rename(columns=rename_map)
gsva_df.index = patient_ids
gsva_df = gsva_df.T
print(f'GSVA shape: {gsva_df.shape}')


print('\n=== 2. Loading TCGA expression ===')
expr = pd.read_csv(TCGA_EXPR, sep='\t', index_col=0)
print(f'TCGA expression: {expr.shape}')


expr_short_map = {}
for c in expr.columns:
    if c.startswith('TCGA-') and len(c.split('-')) >= 4:
        short = '-'.join(c.split('-')[:3])
        expr_short_map[short] = c

common = []
for gsva_id in gsva_df.columns:
    if gsva_id in expr_short_map:
        common.append((gsva_id, expr_short_map[gsva_id]))
print(f'GSVA ↔ TCGA matched: {len(common)} samples')

gsva_aligned = gsva_df[[p[0] for p in common]]
expr_aligned = expr[[p[1] for p in common]]
print(f'Aligned: GSVA {gsva_aligned.shape}, expr {expr_aligned.shape}')


sig = pd.read_csv(SIGNATURE)
sig_modules = sig['Module_Display_Name'].tolist()
sig_genesets = []
for g in sig['Gene_Names']:
    if pd.isna(g):
        sig_genesets.append([])
    else:
        for sep in [';', ',']:
            if sep in g:
                sig_genesets.append([x.strip() for x in g.split(sep) if x.strip()])
                break
        else:
            sig_genesets.append([g.strip()])


def ssgsea_score(expr_vec, indices, alpha=0.25):
    n = len(expr_vec)
    if not indices or len(indices) < 3:
        return np.nan
    in_set = np.zeros(n, dtype=bool)
    in_set[indices] = True
    n_in = in_set.sum()
    n_out = n - n_in
    if n_in == 0 or n_out == 0:
        return np.nan

    sorted_idx = np.argsort(expr_vec)
    sorted_in = in_set[sorted_idx]
    sorted_expr = expr_vec[sorted_idx]

    hit_weights = np.abs(sorted_expr) ** alpha
    hit_sum = np.sum(hit_weights * sorted_in)
    if hit_sum == 0:
        return np.nan
    P_hit = np.cumsum(hit_weights * sorted_in) / hit_sum

    miss_cumsum = np.cumsum(~sorted_in) / n_out

    return float(np.sum(P_hit - miss_cumsum))


def zscore_score(expr_vec, indices):
    if not indices or len(indices) < 3:
        return np.nan

    all_z = (expr_vec - expr_vec.mean()) / (expr_vec.std() + 1e-10)

    return float(all_z[indices].mean())


gene_to_idx = {g: i for i, g in enumerate(expr_aligned.index.tolist())}

ssgsea_scores = np.full((len(sig_modules), len(common)), np.nan)
zscore_scores = np.full((len(sig_modules), len(common)), np.nan)

for mod_idx, (mod_name, genes) in enumerate(zip(sig_modules, sig_genesets)):
    indices = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    if len(indices) < 3:
        continue
    for sample_idx in range(len(common)):
        expr_vec = expr_aligned.iloc[:, sample_idx].values.astype(float)
        ssgsea_scores[mod_idx, sample_idx] = ssgsea_score(expr_vec, indices)
        zscore_scores[mod_idx, sample_idx] = zscore_score(expr_vec, indices)

ssgsea_df = pd.DataFrame(ssgsea_scores, index=sig_modules, columns=[p[0] for p in common])
zscore_df = pd.DataFrame(zscore_scores, index=sig_modules, columns=[p[0] for p in common])

ssgsea_df.to_csv(f'{OUT_DIR}/tcga_ssgsea_scores.csv')
zscore_df.to_csv(f'{OUT_DIR}/tcga_zscore_scores.csv')
print(f'\n  → ssgsea range: {ssgsea_scores[~np.isnan(ssgsea_scores)].min():.3f} to {ssgsea_scores[~np.isnan(ssgsea_scores)].max():.3f}')
print(f'  → zscore range: {zscore_scores[~np.isnan(zscore_scores)].min():.3f} to {zscore_scores[~np.isnan(zscore_scores)].max():.3f}')


print('\n=== 3. Computing correlation ===')
common_modules = [m for m in gsva_aligned.index if m in ssgsea_df.index and m in zscore_df.index]
common_samples = [p[0] for p in common]
print(f'Common: {len(common_modules)} modules × {len(common_samples)} samples')

gsva_sub = gsva_aligned.loc[common_modules, common_samples].values
ssgsea_sub = ssgsea_df.loc[common_modules, common_samples].values
zscore_sub = zscore_df.loc[common_modules, common_samples].values

gsva_ssgsea_corr = np.full(len(common_modules), np.nan)
gsva_zscore_corr = np.full(len(common_modules), np.nan)
ssgsea_zscore_corr = np.full(len(common_modules), np.nan)
for i in range(len(common_modules)):
    if not (np.all(np.isnan(gsva_sub[i])) or np.all(np.isnan(ssgsea_sub[i]))):
        mask = ~(np.isnan(gsva_sub[i]) | np.isnan(ssgsea_sub[i]))
        if mask.sum() > 2:
            gsva_ssgsea_corr[i] = spearmanr(gsva_sub[i, mask], ssgsea_sub[i, mask])[0]
    if not (np.all(np.isnan(gsva_sub[i])) or np.all(np.isnan(zscore_sub[i]))):
        mask = ~(np.isnan(gsva_sub[i]) | np.isnan(zscore_sub[i]))
        if mask.sum() > 2:
            gsva_zscore_corr[i] = spearmanr(gsva_sub[i, mask], zscore_sub[i, mask])[0]
    if not (np.all(np.isnan(ssgsea_sub[i])) or np.all(np.isnan(zscore_sub[i]))):
        mask = ~(np.isnan(ssgsea_sub[i]) | np.isnan(zscore_sub[i]))
        if mask.sum() > 2:
            ssgsea_zscore_corr[i] = spearmanr(ssgsea_sub[i, mask], zscore_sub[i, mask])[0]

print(f'GSVA-ssGSEA mean r: {np.nanmean(gsva_ssgsea_corr):.3f}')
print(f'GSVA-Zscore mean r: {np.nanmean(gsva_zscore_corr):.3f}')
print(f'ssGSEA-Zscore mean r: {np.nanmean(ssgsea_zscore_corr):.3f}')


fig, axes = plt.subplots(1, 3, figsize=(16, 5))
mask = ~(np.isnan(gsva_sub) | np.isnan(ssgsea_sub))
axes[0].scatter(gsva_sub[mask], ssgsea_sub[mask], alpha=0.1, s=1, color='steelblue')
axes[0].set_xlabel('GSVA score'); axes[0].set_ylabel('ssGSEA score')
axes[0].set_title(f'A. GSVA vs ssGSEA\n(mean per-module r = {np.nanmean(gsva_ssgsea_corr):.3f})')
mn, mx = gsva_sub[mask].min(), gsva_sub[mask].max()
axes[0].plot([mn, mx], [mn, mx], 'r--', alpha=0.5); axes[0].grid(alpha=0.3)

mask = ~(np.isnan(gsva_sub) | np.isnan(zscore_sub))
axes[1].scatter(gsva_sub[mask], zscore_sub[mask], alpha=0.1, s=1, color='green')
axes[1].set_xlabel('GSVA score'); axes[1].set_ylabel('Z-score mean')
axes[1].set_title(f'B. GSVA vs Z-score\n(mean per-module r = {np.nanmean(gsva_zscore_corr):.3f})')
mn, mx = gsva_sub[mask].min(), gsva_sub[mask].max()
axes[1].plot([mn, mx], [mn, mx], 'r--', alpha=0.5); axes[1].grid(alpha=0.3)

mask = ~(np.isnan(ssgsea_sub) | np.isnan(zscore_sub))
axes[2].scatter(ssgsea_sub[mask], zscore_sub[mask], alpha=0.1, s=1, color='orange')
axes[2].set_xlabel('ssGSEA score'); axes[2].set_ylabel('Z-score mean')
axes[2].set_title(f'C. ssGSEA vs Z-score\n(mean per-module r = {np.nanmean(ssgsea_zscore_corr):.3f})')
mn, mx = ssgsea_sub[mask].min(), ssgsea_sub[mask].max()
axes[2].plot([mn, mx], [mn, mx], 'r--', alpha=0.5); axes[2].grid(alpha=0.3)

plt.suptitle('Supplementary Figure 13 — FM-AP Scoring Methods Comparison (REAL TCGA expression, n=337 matched samples)', fontsize=14)
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/SuppFig13_scoring_methods_comparison.pdf', format='pdf', pad_inches=0.2)
plt.close()
print(f'\n  → SuppFig13_scoring_methods_comparison.pdf (vector + TrueType)')


pd.DataFrame({
    'Module': common_modules,
    'GSVA_ssGSEA_corr': gsva_ssgsea_corr,
    'GSVA_Zscore_corr': gsva_zscore_corr,
    'ssGSEA_Zscore_corr': ssgsea_zscore_corr,
}).to_csv(f'{OUT_DIR}/SuppFig15_scoring_methods_corr.csv', index=False)
print('Done.')
