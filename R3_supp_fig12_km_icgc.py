

import os, sys, pickle, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings('ignore')

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
DATA_PKL = f"{PROJECT_ROOT}/data/survival_prediction_data.pkl"
MODULE_LIST = f"{PROJECT_ROOT}/data/module_list.csv"
ICGC_SSGSEA = f"{PROJECT_ROOT}/results/icgc_ssgsea_scores.csv"
ICGC_CLIN = f"{PROJECT_ROOT}/data/icgc_clin.csv"
OUT_DIR = f"{PROJECT_ROOT}/results"


print("=== Loading ===")
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
X_tr_f = X_tr_d[:, top_idx]; X_te_d_f = X_te_d[:, top_idx]
selected_features = [feature_names[keep[i]] for i in top_idx]

svm = SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42, class_weight='balanced')
svm.fit(X_tr_f, y_tr)


icgc_ssgsea = pd.read_csv(ICGC_SSGSEA, index_col=0)
icgc_clin = pd.read_csv(ICGC_CLIN)
icgc_clin['ID'] = icgc_clin['ID'].astype(str).str.replace('"', '')

risk_le = LabelEncoder(); risk_le.fit(['low', 'high'])

icgc_X_list = []
icgc_meta = []
for icgc_sample in icgc_ssgsea.index:
    row = []
    for f in selected_features:
        if f in icgc_ssgsea.columns:
            v = icgc_ssgsea.loc[icgc_sample, f]
            try: row.append(float(v) if not pd.isna(v) else 0.0)
            except: row.append(0.0)
        elif f == 'RiskGroup':
            if icgc_sample in icgc_clin['ID'].tolist():
                rv = icgc_clin[icgc_clin['ID'] == icgc_sample]['Risk_Group'].iloc[0]
                try: row.append(float(risk_le.transform([rv])[0]))
                except: row.append(0.0)
            else:
                row.append(0.0)
        else:
            row.append(0.0)
    icgc_X_list.append(row)
    icgc_meta.append(icgc_sample)

icgc_X = np.array(icgc_X_list, dtype=float)
icgc_X = np.nan_to_num(icgc_X, nan=0)


from sklearn.preprocessing import QuantileTransformer
qt = QuantileTransformer(output_distribution='normal', random_state=42, n_quantiles=min(100, X.shape[0]))
qt.fit(X)
icgc_X_66 = np.zeros((icgc_X.shape[0], X.shape[1]))
selected_names_66 = [rename_map[c] for c in feature_cols[:63]] + clinical_features
for i, f in enumerate(selected_features):
    if f in selected_names_66:
        idx_66 = selected_names_66.index(f)
        icgc_X_66[:, idx_66] = icgc_X[:, i]
icgc_X_66_qt = qt.transform(icgc_X_66)
icgc_X_f = icgc_X_66_qt[:, keep][:, top_idx]

proba_icgc = svm.predict_proba(icgc_X_f)[:, 1]


print("\n=== KM analysis ===")
icgc_y, icgc_t = [], []
for s in icgc_meta:
    if s in icgc_clin['ID'].tolist():
        row = icgc_clin[icgc_clin['ID'] == s].iloc[0]
        icgc_y.append(int(row['fustat']))
        icgc_t.append(int(row['futime']))

icgc_y = np.array(icgc_y)
icgc_t = np.array(icgc_t)


median_p = np.median(proba_icgc[:len(icgc_y)])
high_risk_mask = proba_icgc[:len(icgc_y)] >= median_p
low_risk_mask = ~high_risk_mask


lr_result = logrank_test(icgc_t[high_risk_mask], icgc_t[low_risk_mask],
                          icgc_y[high_risk_mask], icgc_y[low_risk_mask])
print(f"Log-rank test: p = {lr_result.p_value:.4e}")


fig, ax = plt.subplots(figsize=(8, 6))
kmf = KaplanMeierFitter()
kmf.fit(icgc_t[high_risk_mask], icgc_y[high_risk_mask], label=f'High Risk (n={high_risk_mask.sum()}, events={icgc_y[high_risk_mask].sum()})')
kmf.plot_survival_function(ax=ax, color='red')
kmf.fit(icgc_t[low_risk_mask], icgc_y[low_risk_mask], label=f'Low Risk (n={low_risk_mask.sum()}, events={icgc_y[low_risk_mask].sum()})')
kmf.plot_survival_function(ax=ax, color='blue')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Survival probability')
ax.set_title(f'Supplementary Figure 12 — KM Survival on ICGC\n(Log-rank p = {lr_result.p_value:.4e}, SVM_rbf predicted)')
ax.legend(loc='lower left')
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/SuppFig12_km_icgc.pdf", format='pdf', pad_inches=0.2)
plt.close()
print(f"  → {OUT_DIR}/SuppFig12_km_icgc.pdf")


pd.DataFrame({
    'ICGC_Sample': icgc_meta[:len(icgc_y)],
    'Risk_Group': ['High' if m else 'Low' for m in high_risk_mask],
    'Proba': proba_icgc[:len(icgc_y)],
    'futime': icgc_t, 'fustat': icgc_y
}).to_csv(f"{OUT_DIR}/SuppFig14_km_icgc_data.csv", index=False)
print(f"  → SuppFig14_km_icgc_data.csv")
