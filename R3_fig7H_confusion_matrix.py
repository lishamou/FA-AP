

import os, pickle, warnings
import numpy as np, pandas as pd
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
warnings.filterwarnings('ignore')

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
OUT_DIR = f"{PROJECT_ROOT}/results"
MAC_DIR = "../../1FIGURES/Figure_supplementary/投稿/Illustrator_assembly/Figure7_拼合素材/r3_subpanels"


with open(f"{PROJECT_ROOT}/data/survival_prediction_data.pkl", 'rb') as f:
    model_data = pickle.load(f)
if model_data['fustat'].dtype == 'O':
    model_data['fustat'] = model_data['fustat'].map({'Alive': 0, 'Dead': 1})

module_map = pd.read_csv(f"{PROJECT_ROOT}/data/module_list.csv").set_index('Feature')['Module'].to_dict()
feature_cols = [c for c in model_data.columns if c.startswith('Feature_')]
rename_map = dict([(c, module_map.get(c, c)) for c in feature_cols])
X = model_data[feature_cols].rename(columns=rename_map).values.astype(float)
for c in ['RiskGroup', 'TP53_Status', 'CTNNB1_Status']:
    if c in model_data.columns:
        le = LabelEncoder()
        X = np.hstack([X, le.fit_transform(model_data[c].astype(str)).reshape(-1, 1)])
y = model_data['fustat'].values.astype(int)

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
            vi = np.var(X_tr_s[:, i]); vj = np.var(X_tr_s[:, j])
            to_drop.add(j if vi >= vj else i)
keep = [i for i in range(X_tr_s.shape[1]) if i not in to_drop]
X_tr_d = X_tr_s[:, keep]; X_te_d = X_te_s[:, keep]
rf = RandomForestClassifier(n_estimators=500, max_depth=10, random_state=42, n_jobs=-1, class_weight='balanced').fit(X_tr_d, y_tr)
top_idx = np.argsort(rf.feature_importances_)[-30:][::-1]
X_tr_f = X_tr_d[:, top_idx]; X_te_f = X_te_d[:, top_idx]

svm = SVC(probability=True, kernel='rbf', C=1.0, gamma='scale', random_state=42, class_weight='balanced')
svm.fit(X_tr_f, y_tr)
y_pred = svm.predict(X_te_f)
cm = confusion_matrix(y_te, y_pred)
tn, fp, fn, tp = cm.ravel()
acc = (tn + tp) / cm.sum()
sens = tp / (tp + fn)
spec = tn / (tn + fp)
ppv = tp / (tp + fp)
f1 = 2 * (ppv * sens) / (ppv + sens)
print(f"SVM_rbf CM (rs=56, n=102): TN={tn}, FP={fp}, FN={fn}, TP={tp}")
print(f"Accuracy={acc:.4f}, Sensitivity={sens:.4f}, Specificity={spec:.4f}, F1={f1:.4f}")


fig, ax = plt.subplots(figsize=(7, 6))


from matplotlib import cm as mpl_cm
blues = mpl_cm.Blues
vmax = cm.max()
norm = plt.Normalize(vmin=0, vmax=vmax)

labels = ['Alive (0)', 'Dead (1)']
for i in range(2):
    for j in range(2):

        color = blues(norm(cm[i, j]))
        rect = Rectangle((j, 1 - i), 1, 1, facecolor=color, edgecolor='white', linewidth=2)
        ax.add_patch(rect)

        ax.text(j + 0.5, 1 - i + 0.55, str(cm[i, j]),
                ha='center', va='center', fontsize=22, fontweight='bold', color='white' if cm[i, j] > vmax/2 else 'black')

        pct = cm[i, j] / cm.sum() * 100
        ax.text(j + 0.5, 1 - i + 0.3, f'({pct:.1f}%)',
                ha='center', va='center', fontsize=10, color='white' if cm[i, j] > vmax/2 else 'dimgray')

ax.set_xlim(0, 2); ax.set_ylim(0, 2)
ax.set_xticks([0.5, 1.5]); ax.set_yticks([1.5, 0.5])
ax.set_xticklabels(labels, fontsize=12, fontweight='bold')
ax.set_yticklabels(labels[::-1], fontsize=12, fontweight='bold')
ax.set_xlabel('Predicted label', fontsize=12, fontweight='bold', labelpad=10)
ax.set_ylabel('True label', fontsize=12, fontweight='bold', labelpad=10)
ax.set_title(f'(H) Confusion matrix: SVM_rbf (R3 best model)\n'
             f'Accuracy = {acc:.4f} | Sensitivity = {sens:.4f} | Specificity = {spec:.4f} | F1 = {f1:.4f} | n = {len(y_te)}',
             fontsize=11)
ax.set_aspect('equal')


fig.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.15)


out_pdf_server = f"{OUT_DIR}/R3_PanelH_CM.pdf"
out_png_server = f"{OUT_DIR}/R3_PanelH_CM.png"
plt.savefig(out_pdf_server, format='pdf', pad_inches=0.2)
plt.savefig(out_png_server, dpi=200, bbox_inches='tight')
plt.close()
print(f"  → {out_pdf_server}")
print(f"  → {out_png_server}")


import subprocess
mac_pdf = f"{MAC_DIR}/R3_PanelH_CM.pdf"
mac_png = f"{MAC_DIR}/R3_PanelH_CM.png"
subprocess.run(['scp', '-q', out_pdf_server, mac_pdf], capture_output=True)
subprocess.run(['scp', '-q', out_png_server, mac_png], capture_output=True)
print(f"  Mac PDF: {mac_pdf}")
print(f"  Mac PNG: {mac_png}")
