

import os, sys, pickle, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['image.composite_image'] = False
import matplotlib.pyplot as plt
import shap
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

warnings.filterwarnings('ignore')

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
DATA_PKL = f"{PROJECT_ROOT}/data/survival_prediction_data.pkl"
MODULE_LIST = f"{PROJECT_ROOT}/data/module_list.csv"
OUT_DIR = f"{PROJECT_ROOT}/results"


print("=== Loading + training SVM_rbf ===")
with open(DATA_PKL, 'rb') as f:
    model_data = pickle.load(f)
if model_data['fustat'].dtype == 'O':
    model_data['fustat'] = model_data['fustat'].map({'Alive': 0, 'Dead': 1})

module_map = pd.read_csv(MODULE_LIST).set_index('Feature')['Module'].to_dict()
feature_cols = [c for c in model_data.columns if c.startswith('Feature_')]
rename_map = {c: module_map.get(c, c) for c in feature_cols}
X_gsva = model_data[feature_cols].rename(columns=rename_map).values.astype(float)
clinical_features = []
for c in ['RiskGroup', 'TP53_Status', 'CTNNB1_Status']:
    if c in model_data.columns:
        le = LabelEncoder()
        X_gsva = np.hstack([X_gsva, le.fit_transform(model_data[c].astype(str)).reshape(-1, 1)])
        clinical_features.append(c)
X = X_gsva
y = model_data['fustat'].values.astype(int)
feature_names = [rename_map[c] for c in feature_cols] + clinical_features

X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, random_state=56, stratify=y)
scaler = StandardScaler()
X_tr_s = scaler.fit_transform(X_tr); X_te_s = scaler.transform(X_te)
corr = np.corrcoef(X_tr_s.T)
to_drop = set()
for i in range(corr.shape[0]):
    if i in to_drop: continue
    for j in range(i+1, corr.shape[1]):
        if j in to_drop: continue
        if abs(corr[i, j]) > 0.8:
            var_i = np.var(X_tr_s[:, i]); var_j = np.var(X_tr_s[:, j])
            to_drop.add(j if var_i >= var_j else i)
keep = [i for i in range(X_tr_s.shape[1]) if i not in to_drop]
X_tr_d = X_tr_s[:, keep]; X_te_d = X_te_s[:, keep]
rf = RandomForestClassifier(n_estimators=500, max_depth=10, random_state=42, n_jobs=-1, class_weight='balanced').fit(X_tr_d, y_tr)
top_idx = np.argsort(rf.feature_importances_)[-30:][::-1]
X_tr_f = X_tr_d[:, top_idx]; X_te_f = X_te_d[:, top_idx]
selected_features = [feature_names[keep[i]] for i in top_idx]

svm = SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42, class_weight='balanced')
svm.fit(X_tr_f, y_tr)


print("\n=== SHAP (KernelExplainer + decision_function) ===")
X_bg = shap.kmeans(X_tr_f, 100)

explainer = shap.KernelExplainer(svm.decision_function, X_bg)
shap_values = explainer.shap_values(X_te_f[:100], nsamples=200)


sv = shap_values if not isinstance(shap_values, list) else shap_values[0]
print(f"SHAP shape: {sv.shape}")
print(f"SHAP values range: [{sv.min():.3f}, {sv.max():.3f}]")
print(f"SHAP mean |value|: {np.abs(sv).mean():.3f}")


print("\n=== Summary plot ===")
plt.figure(figsize=(10, 8))
shap.summary_plot(sv, X_te_f[:100], feature_names=selected_features,
                  show=False, plot_size=None, max_display=15)
plt.title("Supplementary Figure 10A — SHAP Summary (SVM_rbf, n=100 test, decision_function)")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/SuppFig10A_shap_summary.pdf", format='pdf', pad_inches=0.2)
plt.close()
print(f"  → {OUT_DIR}/SuppFig10A_shap_summary.pdf")


labels_bar = [f if f != 'KEGG_COMPLEMENT_AND_COAGULATION_CASCADES' else 'KEGG_COMPLEMENT_AND_\nCOAGULATION_CASCADES' for f in selected_features]
plt.figure(figsize=(10, 8*1.1))
shap.summary_plot(sv, X_te_f[:100], feature_names=labels_bar,
                  plot_type="bar", show=False, plot_size=None, max_display=15)
plt.title("Supplementary Figure 10B — Mean |SHAP| (SVM_rbf, decision_function)")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/SuppFig10B_shap_bar.pdf", format='pdf', pad_inches=0.2)
plt.close()
print(f"  → {OUT_DIR}/SuppFig10B_shap_bar.pdf")


print("\n=== Dependence plots ===")

mean_abs = np.abs(sv).mean(axis=0)
_skip_idx = [i for i, f in enumerate(selected_features) if f == "RiskGroup"]
sort_idx = [i for i in np.argsort(mean_abs)[::-1] if i not in _skip_idx]
top4_idx = sort_idx[:4]
top4 = [selected_features[i] for i in top4_idx]
print(f"Top 4 features by mean |SHAP|: {top4}")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, feat in zip(axes.ravel(), top4):
    if feat in selected_features:
        idx = selected_features.index(feat)

        shap.dependence_plot(idx, sv, X_te_f[:100],
                             feature_names=selected_features, ax=ax, show=False,
                             interaction_index=None)
    else:
        ax.text(0.5, 0.5, f"{feat}\nnot in top 30", ha='center', va='center', transform=ax.transAxes)
plt.suptitle("Supplementary Figure 10C — SHAP Dependence Plots (SVM_rbf)", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/SuppFig10C_shap_dependence.pdf", format='pdf', pad_inches=0.2)
plt.close()
print(f"  → {OUT_DIR}/SuppFig10C_shap_dependence.pdf")

print(f"\n=== SHAP complete ===")
print(f"SHAP values range: [{sv.min():.3f}, {sv.max():.3f}]")
print(f"Top features by mean |SHAP|:")
mean_abs = np.abs(sv).mean(axis=0)
top10_idx = np.argsort(mean_abs)[::-1][:10]
for i in top10_idx:
    print(f"  {selected_features[i]}: {mean_abs[i]:.3f}")
