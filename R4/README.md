# Reproduce the R4 analysis

This directory contains the analysis and figure code for the fourth revision
of manuscript BSPC-D-25-06229, *Functional module profiles of hepatocellular
carcinoma transcriptomes and their associations with mutation status and
survival*. The code is identical to version 4.0.0 in the journal supplement.

## Run the numerical verification

Use Python 3.9 and extract the journal supplement,
`BSPC_R4_Supplementary_Data_and_Code.zip`, to a local directory. That
directory should contain `Analysis_Inputs`, `Reproducibility`, and
`Supplementary_Data`. The analysis-input snapshot and specimen crosswalk
are distributed with the journal supplement.

Run these commands from the repository root. Set `SUPPLEMENT` to the
directory you extracted and choose a new path for `replay_results`.

```bash
SUPPLEMENT=/path/to/extracted-supplement
python -m venv ../fa-ap-r4-venv
source ../fa-ap-r4-venv/bin/activate
python -m pip install -r R4/Reproducibility/requirements.txt
python R4/Reproducibility/reproduce.py --package "$SUPPLEMENT" --out replay_results
python R4/Reproducibility/draw.py --results replay_results \
  --icgc-summary "$SUPPLEMENT/Supplementary_Data/S1_Module_definitions_and_cohort/icgc_manifest.json" \
  --out replay_figures
```

A successful run verifies all 44 numerical tables in the journal supplement
and the analysis summary. The figure command regenerates Figures 1–4 in
PDF and PNG formats. Keep the generated outputs local because the workflow
also produces intermediate predictions by record.

The [verified analysis summary](Reproducibility/summary.json) and
[verification report](Reproducibility/replay_verification.json) record the R4
reference outputs. See [source inputs](Reproducibility/SOURCES.md) and
[input checksums](Reproducibility/input_sha256.txt) to audit the cohort links,
rank scoring, and historical source files. The extracted supplement also
includes the commands for running the workflow from those source files.
