

import os, sys, pandas as pd
import openpyxl

PROJECT_ROOT = "../../../../../.."
SIGNATURE = f"{PROJECT_ROOT}/data/signature.csv"
OUT_DIR = "../../results"

sig = pd.read_csv(SIGNATURE)
print(f"signature.csv: {sig.shape}")
print(f"columns: {sig.columns.tolist()}")


records = []
for i, row in sig.iterrows():
    module_name = row['Module_Display_Name']
    std_name = row['Module_Standard_Name']
    cat1 = row['Category_Level_1']
    cat2 = row['Category_Level_2']
    genes_raw = row['Gene_Names']
    if pd.isna(genes_raw):
        genes_list = []
    else:
        for sep in [';', ',']:
            if sep in genes_raw:
                genes_list = [x.strip() for x in genes_raw.split(sep) if x.strip()]
                break
        else:
            genes_list = [genes_raw.strip()]
    records.append({
        'No': int(row['no']),
        'Module_Name': module_name,
        'Standard_Name': std_name,
        'Category_1': cat1,
        'Category_2': cat2,
        'N_Genes': len(genes_list),
        'Genes': ';'.join(genes_list),
        'Source_Ref': row.get('Source_Reference', ''),
        'Notes': row.get('Consistency_Notes', ''),
    })

df = pd.DataFrame(records)
df.to_csv(f"{OUT_DIR}/SuppTable2_63_modules.csv", index=False)
print(f"  → SuppTable2_63_modules.csv ({len(df)} modules)")


try:
    df.to_excel(f"{OUT_DIR}/SuppTable2_63_modules.xlsx", index=False, engine='openpyxl')
    print(f"  → SuppTable2_63_modules.xlsx")
except Exception as e:
    print(f"  Excel export failed: {e}")


print(f"\n=== Category summary ===")
print(df['Category_1'].value_counts())
print(f"\n=== 5 proprietary modules ===")
proprietary = ['CSC_Activity_PRG', 'LipidMetabolism_Risk_LMRG', 'Metabolism_General_Risk_MRG',
               'Autophagy_Endothelial_HighRisk_AutRG', 'CytokineSignaling_HighRisk_CSIRG']
prop_df = df[df['Module_Name'].isin(proprietary)]
print(prop_df[['Module_Name', 'Category_2', 'N_Genes', 'Genes']].to_string(index=False))
