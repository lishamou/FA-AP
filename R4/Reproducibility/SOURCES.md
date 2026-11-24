# Source inputs for FM-AP R4

The numerical verification route uses the included compact snapshot. The
source-input route uses the original historical files and read-only SHA-256
checks. A current export from the same portal can contain different data.
`input_sha256.txt` identifies each exact local source file.

| File | Upstream resource and preparation |
| --- | --- |
| `TCGA-LIHC_expression_data.txt` | [NCI GDC TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC), archived project expression export in gene-by-aliquot format. The R4 script uses the primary-tumor columns in the archived GSVA matrix to identify the original 371 eligible RNA aliquots, averages duplicate uppercased symbols and ranks genes within each sample. |
| `hcc_tcga_gdc_clinical_data.tsv` | Project clinical export from the [cBioPortal GDC TCGA-LIHC study](https://www.cbioportal.org/study/summary?id=hcc_tcga_gdc). Use rows labeled primary tumor and hepatocellular carcinoma. The analysis reads patient code, overall survival time and status, and AJCC pathologic T. |
| TCGA clinical matrix workbook | Archived annotated workbook based on UCSC Xena `TCGA.LIHC.sampleMap/LIHC_clinicalMatrix` ([Xena data hub](https://tcga.xenahubs.net)). The code reads the third-row header and extracts primary-tumor sample identifiers, age, sex, stage, grade and recorded hepatic risk factors. |
| `alterations_across_samples_TP53.csv`, `alterations_across_samples_CTNNB1.csv` | Archived [cBioPortal](https://www.cbioportal.org/) gene-query exports, with `Patient ID` and the recorded binary `Altered` field. Inspect the actual export before replacing a historical file with a new download. |
| `signature.csv` | R3 archive of the module name and gene lists. The packaged S1 module definitions include the original lists and the lists after measured-gene coverage checks. The analysis snapshot contains the identical 63 sets for computational replay. |
| `gsva_scores.csv` | Archived TCGA GSVA sensitivity input from the earlier analysis. The exact analysis values are preserved in `Analysis_Inputs/tcga_gsva.npz`. The full matrix also supplies primary-tumor aliquot identifiers for the source-input route. |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | [NCI GDC PanCanAtlas](https://gdc.cancer.gov/about-data/publications/pancanatlas) [ABSOLUTE source file](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5). Pass the original TSV. The analysis matches DNA and RNA aliquots at the tissue-portion level. |
| `23.ICGC-LIRI-JP-sampleExp.txt`, `24.ICGC-LIRI-JP-time.txt` | Historical project exports of LIRI-JP data. The upstream resource is [Legacy ICGC 25K](https://docs.icgc-argo.org/docs/data-access/icgc-25k-data). These are frozen project exports. The current ICGC download endpoint does not provide a verified byte-identical replacement under these filenames. |
| `bulk_ICGC_samples.csv`, `bulk_ICGC_survival_input.csv` | Derived from [HCCDB v2.0](https://doi.org/10.1093/gpbjnl/qzae011), dataset HCCDB18. The extraction is described below. R4 uses these files to validate tumor assignment at donor level and to obtain clinical covariates. |
| `TCGA_time.txt` | Archived TCGA overall survival data in days. Used to audit the month-to-day conversion on 341 overlapping records and is separate from the modeling inputs. |

## HCCDB18 extraction and linkage

Read the HCCDB18 sample and patient annotation tables as transposed metadata
matrices. The first column of each source table lists field names, and the
other columns identify records. Locate the single row for NCEH1 in the gene
matrix, accepting the aliases AADACL1 and KIAA1363 when needed. Match the sample
annotation ID against an expression column, keep finite NCEH1 values, and
retain PATIENT, TYPE, SAMPLE_NAME and NCEH1 for the specimen export. For the
patient export, retain tumor-linked PATIENT_ID records with positive SUR
(months) and known STATUS. Join the patient annotation by PATIENT_ID. Keep
PATIENT, AGE, GENDER, SUR, TNM_STAGE_T and event (Dead = 1).

R4 `prepare_icgc.py` extracts donor IDs from the legacy expression headers,
keeps donors with a unique tumor-suffix specimen and a unique annotated HCC
specimen in HCCDB18, and cross-checks paired NCEH1 and outcomes. It imports
follow-up days from the legacy ICGC file. Donor-level linkage and expression
concordance support the tissue selection. The legacy and HCCDB identifiers
are maintained in separate specimen naming systems.

## Compact snapshot

The snapshot retains gene ranks only for genes used by the 63 module sets,
after ranking against the full TCGA background. It includes the archived
GSVA values and the already calculated common-background transfer scores.
It preserves the matched purity and the minimum clinical and event fields
needed for the analysis. Its row keys replace the original patient barcodes.
Raw expression, clinical source files and specimen-level joins are upstream
of the snapshot. For auditing those transformations, use the source route.
