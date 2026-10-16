

import os, pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

BASE = "../.."
OUT_DIR = f"{BASE}/results"
DST = f"{BASE}/1FIGURES/Figure_supplementary/投稿/Illustrator_assembly/Figure7_拼合素材/r3_subpanels"

perf = pd.read_csv(f"{OUT_DIR}/3_model_performance_test_v31.csv")


paper_legend = ['LogisticRegression', 'SVM_rbf', 'RandomForest', 'XGBoost', 'MLP', 'LightGBM', 'CatBoost']
display_names = {
    'LogisticRegression': 'Logistic Regression',
    'SVM_rbf': 'SVM (RBF)',
    'RandomForest': 'Random Forest',
    'XGBoost': 'XGBoost',
    'MLP': 'MLP',
    'LightGBM': 'LightGBM',
    'CatBoost': 'CatBoost',
}
sub = perf[perf['Model'].isin(paper_legend)].copy()
sub['Display'] = sub['Model'].map(display_names)
sub = sub.set_index('Display').loc[list(display_names.values())].reset_index()
print(sub[['Display', 'AUC_train', 'AUC_test']].to_string(index=False))


fig, ax = plt.subplots(figsize=(12, 6))
models = sub['Display'].tolist()
x = np.arange(len(models))
width = 0.4

ax.bar(x - width/2, sub['AUC_train'], width, label='Train',
       color='#08519c', edgecolor='white', linewidth=0.5)
ax.bar(x + width/2, sub['AUC_test'], width, label='Test',
       color='#6baed6', edgecolor='white', linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels(models, rotation=0, ha='center', fontsize=11)
ax.set_ylabel('AUC', fontsize=12)
ax.set_ylim(0, 1.05)
ax.axhline(0.7, color='red', linestyle='--', alpha=0.5, label='AUC=0.7')
ax.set_title('(E) Training vs Test AUC of 7 baseline models (Logistic Regression, SVM, Random Forest, XGBoost, MLP, LightGBM, CatBoost)\n(SVM best: test AUC = 0.7289, 95% CI 0.6318-0.8299)', fontsize=11)
ax.legend(loc='lower right', fontsize=10)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/R3_PanelE_TrainTest.png", dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  → R3_PanelE_TrainTest.png (7 baseline, blue scheme, x-axis full names)")
