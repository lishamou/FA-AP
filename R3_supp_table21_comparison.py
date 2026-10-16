import os, pandas as pd

PROJECT_ROOT = "../../../../../.."
OUT_DIR = "../../results"

models = [
    {
        'Model': 'Zhang et al. (2026, Discov Oncol) [31]',
        'Year': 2026,
        'Approach': 'Palmitoylation-related genes',
        'TCGA_AUC': 'N/A',
        'Validation_AUC': '0.68-0.72',
        'Cohort': 'TCGA-LIHC + ICGC',
        'Method': 'Cox + LASSO',
        'N_features': 14,
        'R3_FM_AP_AUC': '0.7289 (TCGA) / 0.6525 (ICGC)',
    },
    {
        'Model': 'Yao et al. (2026, Faseb j) [32]',
        'Year': 2026,
        'Approach': 'Sialylation-related genes',
        'TCGA_AUC': 'N/A',
        'Validation_AUC': '0.70',
        'Cohort': 'TCGA-LIHC',
        'Method': 'Single-cell + Bulk integrated',
        'N_features': 12,
        'R3_FM_AP_AUC': '0.7289 / 0.6525',
    },
    {
        'Model': 'Shen et al. (2025, Front Immunol) [33]',
        'Year': 2025,
        'Approach': 'Pyroptosis-related genes',
        'TCGA_AUC': 'N/A',
        'Validation_AUC': '0.68',
        'Cohort': 'TCGA-LIHC',
        'Method': 'Single-cell RNA-seq + Cox',
        'N_features': 18,
        'R3_FM_AP_AUC': '0.7289 / 0.6525',
    },
    {
        'Model': 'Liu et al. (2025, Recent Pat) [34]',
        'Year': 2025,
        'Approach': 'Autophagy-related gene signature',
        'TCGA_AUC': 'N/A',
        'Validation_AUC': '0.69',
        'Cohort': 'TCGA-LIHC',
        'Method': 'Bioinformatics + Cox',
        'N_features': 15,
        'R3_FM_AP_AUC': '0.7289 / 0.6525',
    },
    {
        'Model': 'Lin et al. (2025, Discov Oncol) [35]',
        'Year': 2025,
        'Approach': 'Immune-metabolic prognostic model',
        'TCGA_AUC': 'N/A',
        'Validation_AUC': '0.70',
        'Cohort': 'TCGA-LIHC',
        'Method': 'Machine learning ensemble',
        'N_features': 20,
        'R3_FM_AP_AUC': '0.7289 / 0.6525',
    },
    {
        'Model': 'R3 Ours (FM-AP Stacking / SVM_rbf)',
        'Year': 2026,
        'Approach': 'Knowledge-based FM-AP (63 modules)',
        'TCGA_AUC': '0.7289 (test set, n=102)',
        'Validation_AUC': '0.6525 (ICGC, n=445, ssGSEA)',
        'Cohort': 'TCGA-LIHC + ICGC LIRI-JP',
        'Method': 'SVM_rbf + 8 baseline models',
        'N_features': '27 (selected) from 66 (63 GSVA + 3 clinical)',
        'R3_FM_AP_AUC': '0.7289 / 0.6525',
    },
]

df = pd.DataFrame(models)
df.to_csv(f"{OUT_DIR}/SuppTable21_recent_models_comparison.csv", index=False)
print(f"  → SuppTable21_recent_models_comparison.csv ({len(df)} models)")

print("\n=== Supp Table 21: Recent HCC prognostic models ===")
print(df.to_string(index=False))

try:
    df.to_excel(f"{OUT_DIR}/SuppTable21_recent_models_comparison.xlsx", index=False, engine='openpyxl')
    print(f"  → SuppTable21_recent_models_comparison.xlsx")
except Exception as e:
    print(f"  Excel export failed: {e}")
