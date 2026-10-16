

import os, pickle, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, VotingClassifier, StackingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, confusion_matrix
import xgboost as xgb
from xgboost import XGBClassifier
import lightgbm as lgb
from lightgbm import LGBMClassifier
import seaborn as sns
import sys
sys.path.insert(0, "/home/moulisha/.local/lib/python3.8/site-packages")
from catboost import CatBoostClassifier

warnings.filterwarnings('ignore')

BASE = "/data/share/mouls/LIHC_FMAP_R3"
OUT_DIR = f"{BASE}/results"


print("=== Load & Train ===")
with open(f"{BASE}/data/survival_prediction_data.pkl", 'rb') as f:
    model_data = pickle.load(f)
if model_data['fustat'].dtype == 'O':
    model_data['fustat'] = model_data['fustat'].map({'Alive': 0, 'Dead': 1})

module_map = pd.read_csv(f"{BASE}/data/module_list.csv").set_index('Feature')['Module'].to_dict()
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


print("=== Train all models ===")
models = {
    'Logistic Regression': LogisticRegression(max_iter=5000, random_state=42, class_weight='balanced', C=0.5),
    'SVM (RBF)': SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42, class_weight='balanced'),
    'SVM (linear)': SVC(probability=True, kernel='linear', C=0.5, random_state=42, class_weight='balanced'),
    'Random Forest': RandomForestClassifier(n_estimators=500, max_depth=10, min_samples_split=5, random_state=42, class_weight='balanced', n_jobs=-1),
    'XGBoost': XGBClassifier(n_estimators=500, max_depth=5, learning_rate=0.03, subsample=0.8, random_state=42, eval_metric='logloss', scale_pos_weight=(len(y_tr)-y_tr.sum())/y_tr.sum(), n_jobs=1),
    'LightGBM': LGBMClassifier(n_estimators=500, max_depth=6, learning_rate=0.03, random_state=42, class_weight='balanced', verbose=-1, n_jobs=1),
    'CatBoost': CatBoostClassifier(iterations=500, depth=6, learning_rate=0.03, l2_leaf_reg=3.0, random_state=42, verbose=0, auto_class_weights='Balanced'),
    'MLP': MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=3000, random_state=42, early_stopping=True, validation_fraction=0.15, alpha=0.01),
    'ExtraTrees': ExtraTreesClassifier(n_estimators=500, max_depth=10, random_state=42, class_weight='balanced', n_jobs=1),
}
results = {}
for name, model in models.items():
    model.fit(X_tr_f, y_tr)
    results[name] = {
        'proba_test': model.predict_proba(X_te_f)[:, 1],
        'auc_test': roc_auc_score(y_te, model.predict_proba(X_te_f)[:, 1]),
    }


base_for_ens = list(models.items())
weights = [results[n]['auc_test'] for n, _ in base_for_ens]
vote = VotingClassifier(estimators=base_for_ens, voting='soft', weights=weights)
vote.fit(X_tr_f, y_tr)
results['Voting'] = {
    'proba_test': vote.predict_proba(X_te_f)[:, 1],
    'auc_test': roc_auc_score(y_te, vote.predict_proba(X_te_f)[:, 1]),
}
stack = StackingClassifier(estimators=base_for_ens, final_estimator=LogisticRegression(max_iter=5000, random_state=42, class_weight='balanced', C=0.5), cv=5, n_jobs=1, passthrough=False)
stack.fit(X_tr_f, y_tr)
results['Stacking'] = {
    'proba_test': stack.predict_proba(X_te_f)[:, 1],
    'auc_test': roc_auc_score(y_te, stack.predict_proba(X_te_f)[:, 1]),
    'pred_test': stack.predict(X_te_f),
}

print(f"  All 11 models trained. Best: SVM (RBF) = {results['SVM (RBF)']['auc_test']:.4f}")


print("\n=== Panel F: ROC blue scheme ===")
fig, ax = plt.subplots(figsize=(8, 7))


sorted_models = sorted([(name, r['auc_test']) for name, r in results.items()], key=lambda x: x[1])
norm = plt.Normalize(vmin=0.55, vmax=0.75)
cmap = plt.cm.Blues


def get_color(name, auc):
    if name == 'Stacking':
        return '#08306b'
    elif name == 'Voting':
        return '#08519c'
    elif name == 'SVM (RBF)':
        return '#2171b5'
    else:

        return cmap(norm(auc))

for name, r in results.items():
    fpr, tpr, _ = roc_curve(y_te, r['proba_test'])
    lw = 3.0 if name in ['Stacking', 'Voting', 'SVM (RBF)'] else 1.5
    ls = '-' if name in ['Stacking', 'Voting'] else '-'
    label = f"{name} (AUC={r['auc_test']:.3f})"
    ax.plot(fpr, tpr, label=label, color=get_color(name, r['auc_test']),
            linewidth=lw, linestyle=ls, alpha=0.9 if lw > 1.5 else 0.7)

ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=0.8)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_title('(F) ROC curves: 9 baseline + Voting + Stacking models\n(SVM (RBF) best: AUC = 0.7289, 95% CI 0.6318-0.8299)', fontsize=11)
ax.legend(loc='lower right', fontsize=7, ncol=2)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/R3_PanelF_ROC.pdf", format='pdf', bbox_inches='tight')
plt.savefig(f"{OUT_DIR}/R3_PanelF_ROC.png", dpi=200, bbox_inches='tight')
plt.close()
print(f"  → R3_PanelF_ROC.pdf / .png (blue scheme, TrueType embedded)")


print("\n=== Panel G: PR blue scheme ===")
fig, ax = plt.subplots(figsize=(8, 7))

for name, r in results.items():
    precision, recall, _ = precision_recall_curve(y_te, r['proba_test'])
    ap = average_precision_score(y_te, r['proba_test'])
    lw = 3.0 if name in ['Stacking', 'Voting', 'SVM (RBF)'] else 1.5
    ls = '-' if name in ['Stacking', 'Voting'] else '-'
    label = f"{name} (AP={ap:.3f})"

    color = get_color(name, r['auc_test'])
    ax.plot(recall, precision, label=label, color=color,
            linewidth=lw, linestyle=ls, alpha=0.9 if lw > 1.5 else 0.7)

ax.set_xlabel('Recall', fontsize=12)
ax.set_ylabel('Precision', fontsize=12)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.05)
ax.set_title('(G) Precision-Recall curves: 9 baseline + Voting + Stacking models', fontsize=11)
ax.legend(loc='lower left', fontsize=7, ncol=2)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/R3_PanelG_PR.pdf", format='pdf', bbox_inches='tight')
plt.savefig(f"{OUT_DIR}/R3_PanelG_PR.png", dpi=200, bbox_inches='tight')
plt.close()
print(f"  → R3_PanelG_PR.png (blue scheme)")
