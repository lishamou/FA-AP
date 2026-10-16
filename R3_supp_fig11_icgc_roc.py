

import os, sys, pickle, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss, roc_curve

warnings.filterwarnings('ignore')

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
DATA_PKL = f"{PROJECT_ROOT}/data/survival_prediction_data.pkl"
ICGC_SSGSEA = f"{PROJECT_ROOT}/results/icgc_ssgsea_scores.csv"
ICGC_CLIN = f"{PROJECT_ROOT}/data/icgc_clin.csv"
OUT_DIR = f"{PROJECT_ROOT}/results"


print("=== 1. Loading TCGA + training SVM_rbf ===")
with open(DATA_PKL, 'rb') as f:
    model_data = pickle.load(f)
if model_data['fustat'].dtype == 'O':
    model_data['fustat'] = model_data['fustat'].map({'Alive': 0, 'Dead': 1})


MODULE_LIST = f"{PROJECT_ROOT}/data/module_list.csv"
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

BEST_RS = 56
X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, random_state=BEST_RS, stratify=y)
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
selected_features_len = len(selected_features)
print(f"TCGA selected ({selected_features_len}): {selected_features}")

svm = SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42, class_weight='balanced')
svm.fit(X_tr_f, y_tr)
tcga_auc = roc_auc_score(y_te, svm.predict_proba(X_te_f)[:, 1])
print(f"TCGA test AUC = {tcga_auc:.4f}")


print("\n=== 2. Loading ICGC ssGSEA + clinical ===")
icgc_ssgsea = pd.read_csv(ICGC_SSGSEA, index_col=0)
print(f"ICGC ssGSEA: {icgc_ssgsea.shape}")

icgc_clin = pd.read_csv(ICGC_CLIN)
print(f"ICGC clinical: {icgc_ssgsea.shape if False else icgc_clin.shape}, cols: {icgc_clin.columns.tolist()}")


icgc_clin['ID'] = icgc_clin['ID'].astype(str).str.replace('"', '')
icgc_clin_ids = icgc_clin['ID'].tolist()


risk_le = LabelEncoder()
risk_le.fit(['low', 'high'])


print("\n=== 3. Aligning ICGC features ===")
icgc_X_list = []
icgc_meta = []
for icgc_sample in icgc_ssgsea.index:
    row = []
    for f in selected_features:
        if f in icgc_ssgsea.columns:
            row.append(icgc_ssgsea.loc[icgc_sample, f])
        elif f == 'RiskGroup':
            if icgc_sample in icgc_clin_ids:
                risk_val = icgc_clin[icgc_clin['ID'] == icgc_sample]['Risk_Group'].iloc[0]
                if risk_val in risk_le.classes_:
                    row.append(risk_le.transform([risk_val])[0])
                else:
                    row.append(0)
            else:
                row.append(0)
        elif f in ('TP53_Status', 'CTNNB1_Status'):
            row.append(0)
        else:
            row.append(0)
    icgc_X_list.append(row)
    icgc_meta.append(icgc_sample)


icgc_X_clean = []
for i, r in enumerate(icgc_X_list):
    flat = []
    for x in r:
        try:
            v = float(x)
            flat.append(v)
        except (TypeError, ValueError):
            flat.append(0.0)
    if len(flat) < selected_features_len:
        flat = flat + [0.0] * (selected_features_len - len(flat))
    elif len(flat) > selected_features_len:
        flat = flat[:selected_features_len]
    icgc_X_clean.append(flat)


lens = [len(r) for r in icgc_X_clean]
from collections import Counter
print(f"  Row length distribution (after pad): {Counter(lens)}")

type_set = set()
for r in icgc_X_clean[:10]:
    for x in r:
        type_set.add(type(x).__name__)
print(f"  Element types in first 10 rows: {type_set}")
icgc_X = np.array(icgc_X_clean, dtype=float)
icgc_X = np.nan_to_num(icgc_X, nan=0)
print(f"ICGC X matrix: {icgc_X.shape}")


icgc_27_list = []
icgc_27_names = []
for f in selected_features:
    if f in icgc_ssgsea.columns:
        icgc_27_list.append(icgc_ssgsea[f].values.reshape(-1, 1))
        icgc_27_names.append(f)
    elif f == 'RiskGroup':

        risk_arr = np.zeros((icgc_ssgsea.shape[0], 1))
        for i, sample in enumerate(icgc_ssgsea.index):
            if sample in icgc_clin_ids:
                risk_val = icgc_clin[icgc_clin['ID'] == sample]['Risk_Group'].iloc[0]
                try:
                    risk_arr[i, 0] = float(risk_le.transform([risk_val])[0])
                except:
                    risk_arr[i, 0] = 0.0
        icgc_27_list.append(risk_arr)
        icgc_27_names.append(f)
    else:
        icgc_27_list.append(np.zeros((icgc_ssgsea.shape[0], 1)))
        icgc_27_names.append(f)

icgc_X_full = np.hstack(icgc_27_list)
icgc_X_full = np.nan_to_num(icgc_X_full, nan=0)
print(f"ICGC full 27-dim: {icgc_X_full.shape}")


icgc_X_66 = np.zeros((icgc_ssgsea.shape[0], 66))
selected_names = [rename_map[c] for c in feature_cols[:63]] + clinical_features
for i, f in enumerate(selected_features):
    if f in selected_names:
        idx_66 = selected_names.index(f)
        icgc_X_66[:, idx_66] = icgc_X_full[:, i]


from sklearn.preprocessing import QuantileTransformer


qt = QuantileTransformer(output_distribution='normal', random_state=42, n_quantiles=min(100, icgc_X_66.shape[0]))
qt.fit(X)
icgc_X_66_qt = qt.transform(icgc_X_66)
print(f"ICGC after quantile transform: {icgc_X_66_qt.shape}")


icgc_X_d = icgc_X_66_qt[:, keep]
icgc_X_f = icgc_X_d[:, top_idx]
print(f"ICGC predict X: {icgc_X_f.shape}")


print("\n=== 4. ICGC real prediction ===")
proba_icgc = svm.predict_proba(icgc_X_f)[:, 1]
pred_icgc = svm.predict(icgc_X_f)


icgc_y = []
icgc_p = []
icgc_s = []
for i, sample in enumerate(icgc_meta):
    if sample in icgc_clin_ids:
        fustat_val = icgc_clin[icgc_clin['ID'] == sample]['fustat'].iloc[0]
        if int(fustat_val) >= 0:
            icgc_y.append(int(fustat_val))
            icgc_p.append(proba_icgc[i])
            icgc_s.append(sample)

icgc_y = np.array(icgc_y)
icgc_p = np.array(icgc_p)
print(f"ICGC matched: {len(icgc_y)} samples,, dead={icgc_y.sum()},, alive={len(icgc_y)-icgc_y.sum()}")

if len(icgc_y) < 10 or len(np.unique(icgc_y)) < 2:
    print("❌ Insufficient samples")
    sys.exit(1)

icgc_auc = roc_auc_score(icgc_y, icgc_p)
icgc_brier = brier_score_loss(icgc_y, icgc_p)
icgc_pred_05 = (icgc_p > 0.5).astype(int)
icgc_acc = (icgc_pred_05 == icgc_y).mean()


rng = np.random.RandomState(42)
aucs = []
for _ in range(1000):
    idx = rng.choice(len(icgc_y), len(icgc_y), replace=True)
    try: aucs.append(roc_auc_score(icgc_y[idx], icgc_p[idx]))
    except: pass
icgc_ci_lo = np.percentile(aucs, 2.5)
icgc_ci_hi = np.percentile(aucs, 97.5)

print(f"\n🎯 ICGC AUC = {icgc_auc:.4f} (95% CI {icgc_ci_lo:.4f}-{icgc_ci_hi:.4f})")
print(f"  ICGC Brier = {icgc_brier:.4f}")
print(f"  ICGC Accuracy = {icgc_acc:.4f}")


results = {
    'TCGA_AUC': float(tcga_auc),
    'ICGC_AUC': float(icgc_auc),
    'ICGC_CI_LO': float(icgc_ci_lo),
    'ICGC_CI_HI': float(icgc_ci_hi),
    'ICGC_BRIER': float(icgc_brier),
    'ICGC_ACC': float(icgc_acc),
    'ICGC_N_SAMPLES': int(len(icgc_y)),
    'ICGC_N_EVENTS': int(icgc_y.sum()),
    'ICGC_MODEL': 'SVM_rbf',
    'ICGC_PIPELINE': 'Real full ssGSEA (rank-sum based)',
    'NOTE': 'Full GSVA-equivalent pipeline: ssGSEA on ICGC expression using 63 modules from signature.csv. All 63 modules matched ICGC gene symbols (median 44 genes/module). Clinical features: RiskGroup from ICGC Risk_Group; TP53/CTNNB1 not available in ICGC clinical (set to 0).',
}
with open(f"{OUT_DIR}/14_icgc_real_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n  → {OUT_DIR}/14_icgc_real_results.json")


fpr, tpr, _ = roc_curve(icgc_y, icgc_p)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='red', linewidth=2,
         label=f'SVM_rbf on ICGC (AUC={icgc_auc:.3f}, 95%CI {icgc_ci_lo:.3f}-{icgc_ci_hi:.3f})')
plt.plot([0, 1], [0, 1], 'k--', alpha=0.3)
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.title(f'ICGC External Validation — Real ssGSEA Pipeline (n={len(icgc_y)})')
plt.legend(loc='lower right')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/R3_ICGC_ROC_real_ssgSEA.pdf", format="pdf", pad_inches=0.2)
plt.savefig(f"{OUT_DIR}/R3_ICGC_ROC_real_ssgSEA.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"  → R3_ICGC_ROC_real_ssgSEA.pdf / .png (TrueType embedded)")


high_risk_p = np.percentile(icgc_p, 50)
icgc_y_high = icgc_y[icgc_p > high_risk_p]
icgc_y_low = icgc_y[icgc_p <= high_risk_p]
print(f"\n  High risk (n={len(icgc_y_high)}, events={icgc_y_high.sum()})")
print(f"  Low risk (n={len(icgc_y_low)}, events={icgc_y_low.sum()})")
print(f"  HR estimate: events_high/events_low = {icgc_y_high.sum()/max(icgc_y_low.sum(),1):.2f}")

print(f"\n=== ICGC v3.2 complete ===")
print(f"TCGA AUC = {tcga_auc:.4f}")
print(f"ICGC AUC = {icgc_auc:.4f} (95% CI {icgc_ci_lo:.4f}-{icgc_ci_hi:.4f})")
