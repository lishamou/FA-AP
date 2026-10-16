import os
import pickle
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from lifelines.statistics import proportional_hazard_test
import scipy.stats as stats
import openpyxl
from openpyxl import Workbook

PROJECT_ROOT = "/data/share/mouls/LIHC_FMAP_R3"
SUPP_DIR = f"{PROJECT_ROOT}/results/Supplementary_Tables"
os.makedirs(SUPP_DIR, exist_ok=True)

with open(f"{PROJECT_ROOT}/data/survival_prediction_data.pkl", "rb") as f:
    data = pickle.load(f)
modules = pd.read_csv(f"{PROJECT_ROOT}/data/module_list.csv")
feat_to_mod = dict(zip(modules["Feature"], modules["Module"]))
print(f"Loaded data: {data.shape[0]} samples × {data.shape[1]} cols")
print(f"Module mapping: {len(modules)} rows ({len(set(feat_to_mod.values()))} unique modules)")

df = pd.DataFrame()
df["time"] = data["futime"].astype(float)
df["event"] = (data["fustat"] == "Dead").astype(int)

selected_features = [f"Feature_{i}" for i in range(6, 64)]
selected_modules = [feat_to_mod[f] for f in selected_features]
for f, m in zip(selected_features, selected_modules):
    df[m] = data[f].astype(float)
df = df.dropna()
print(f"After dropna: {df.shape[0]} samples; event rate = {df['event'].mean():.3f}")
print(f"Modules (n={len(selected_modules)}): {selected_modules[:5]} ... {selected_modules[-3:]}")

cph = CoxPHFitter()
cph.fit(df, duration_col="time", event_col="event")

summary_df = cph.summary
s17 = summary_df.reset_index().rename(columns={"index": "covariate"})
s17["covariate"] = s17["covariate"].astype(str)

expected_cols = ["covariate", "coef", "exp(coef)", "se(coef)",
                 "coef lower 95%", "coef upper 95%",
                 "exp(coef) lower 95%", "exp(coef) upper 95%",
                 "cmp to", "z", "p", "-log2(p)"]
s17 = s17[expected_cols]
print(f"\nS17 rows: {len(s17)}")

assert len(s17) == 58, f"Expected 58 covariates, got {len(s17)}"
assert "KEGG_CELL_CYCLE" in s17["covariate"].values
assert "Type_II_IFN_Reponse" in s17["covariate"].values

N = len(df)
events = int(df["event"].sum())
censored = N - events
c_index = concordance_index(df["time"], -cph.predict_partial_hazard(df).values, df["event"])

ph = proportional_hazard_test(cph, df, time_transform="rank")
chi2_global = float(ph.summary["test_statistic"].sum())
dof_global = len(ph.summary)
p_global = 1 - stats.chi2.cdf(chi2_global, dof_global)

s18_rows = [
    {"Section": "Dataset & Model Overview", "Metric": "Total samples (N)", "Value": N},
    {"Section": "Dataset & Model Overview", "Metric": "Death events", "Value": events},
    {"Section": "Dataset & Model Overview", "Metric": "Censored", "Value": censored},
    {"Section": "Dataset & Model Overview", "Metric": "Censoring rate", "Value": round(censored / N, 3)},
    {"Section": "Dataset & Model Overview", "Metric": "Features (modules)", "Value": len(selected_modules)},
    {"Section": "Dataset & Model Overview", "Metric": "C-index", "Value": round(c_index, 4)},
    {"Section": "Dataset & Model Overview", "Metric": "Significant covariates (p < 0.05)", "Value": int((s17["p"] < 0.05).sum())},
    {"Section": "PH Assumption Test", "Metric": "Test", "Value": "Schoenfeld Residuals Global Test"},
    {"Section": "PH Assumption Test", "Metric": "Statistic (χ²)", "Value": round(chi2_global, 4)},
    {"Section": "PH Assumption Test", "Metric": "p-value", "Value": round(p_global, 4)},
    {"Section": "PH Assumption Test", "Metric": "Conclusion", "Value": "PH assumption satisfied" if p_global > 0.05 else "PH assumption violated"},
    {"Section": "PH Assumption Test", "Metric": "Interpretation", "Value": "Cox PH model is appropriate" if p_global > 0.05 else "Stratified or time-varying model recommended"},
]
s18 = pd.DataFrame(s18_rows)
print(f"\nS18:")
print(s18.to_string(index=False))

s17.to_csv(f"{SUPP_DIR}/Table S19_continuous_cox_R3_rerun.csv", index=False)
s18.to_csv(f"{SUPP_DIR}/Table S20_cox_summary_R3_rerun.csv", index=False)

def write_xlsx(csv_path, xlsx_path, sheet_name):
    df_in = pd.read_csv(csv_path)
    wb = Workbook()
    ws = wb.active
    ws.title = sheet_name
    ws.append(list(df_in.columns))
    for _, row in df_in.iterrows():
        ws.append([row[c] if not pd.isna(row[c]) else None for c in df_in.columns])
    wb.save(xlsx_path)

write_xlsx(f"{SUPP_DIR}/Table S19_continuous_cox_R3_rerun.csv",
           f"{SUPP_DIR}/Table S19_continuous_cox_R3_rerun.xlsx",
           "Table S19_continuous_cox")
write_xlsx(f"{SUPP_DIR}/Table S20_cox_summary_R3_rerun.csv",
           f"{SUPP_DIR}/Table S20_cox_summary_R3_rerun.xlsx",
           "Table S20_cox_summary")

print(f"\n✓ Wrote:")
for f in sorted(os.listdir(SUPP_DIR)):
    if "R3_rerun" in f and ("S17" in f or "S18" in f):
        print(f"  {SUPP_DIR}/{f}")
