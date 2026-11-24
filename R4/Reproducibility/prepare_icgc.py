"""Cross-check legacy ICGC tissue labels against a separately annotated HCCDB extract."""
import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

parser = argparse.ArgumentParser()
parser.add_argument('data',type=Path,help='Project input directory')
parser.add_argument('hccdb',type=Path,help='Directory with annotated HCCDB18 per-specimen and clinical extracts')
parser.add_argument('out',type=Path)
a=parser.parse_args()
times=pd.read_csv(a.data/'24.ICGC-LIRI-JP-time.txt',sep='\t').set_index('id',verify_integrity=True)
expr_path=a.data/'23.ICGC-LIRI-JP-sampleExp.txt'
with expr_path.open() as stream:
    header=stream.readline().rstrip().split('\t')[1:]
    expr_row=next(line for line in stream if line.startswith('NCEH1\t')).rstrip().split('\t')[1:]
spec=pd.read_csv(a.hccdb/'bulk_ICGC_samples.csv')
clinical=pd.read_csv(a.hccdb/'bulk_ICGC_survival_input.csv').set_index('PATIENT',verify_integrity=True)
dat=pd.DataFrame({'sample_id':header,'NCEH1_legacy':pd.to_numeric(expr_row)})
dat['patient_id']=dat.sample_id.str.extract(r'^(?:SP\d+)-(DO\d+)(?:-T)?$')[0]
assert dat.patient_id.notna().all() and dat.sample_id.is_unique
dat['t_suffix']=dat.sample_id.str.endswith('-T')
dat['t_count']=dat['t_suffix'].groupby(dat.patient_id).transform('sum')
spec['id_count']=spec.groupby(['PATIENT','TYPE']).PATIENT.transform('size')
l=dat.pivot_table(index='patient_id',columns='t_suffix',values='NCEH1_legacy',aggfunc='mean').rename(columns={False:'legacy_unmarked',True:'legacy_T'})
h=spec.pivot_table(index='PATIENT',columns='TYPE',values='NCEH1',aggfunc='mean')
pair=l.join(h).dropna()
agreement={}
for n, col in [('legacy_T','HCC'),('legacy_unmarked','Adjacent')]:
    valid=l.join(h[[col]]).dropna(subset=[n,col])
    agreement[n]={'n':len(valid),'rho':float(spearmanr(valid[n],valid[col]).statistic)}
agreement['paired_contrast']={'n':len(pair),'rho':float(spearmanr(pair.legacy_T-pair.legacy_unmarked,pair.HCC-pair.Adjacent).statistic)}
print(agreement)
tum=spec.loc[spec.TYPE.eq('HCC') & spec.id_count.eq(1),['PATIENT','SAMPLE_NAME','NCEH1']].set_index('PATIENT',verify_integrity=True)
one_tumor=dat.loc[dat.t_suffix & dat.t_count.eq(1)].set_index('patient_id',verify_integrity=True)
by_annotation=one_tumor.join(tum,how='inner')
by_clinical=by_annotation.join(clinical[['AGE','GENDER','SUR','event','TNM_STAGE_T']],how='inner')
cohort=by_clinical.join(times,how='inner')
assert (cohort.fustat==cohort.event).all(), 'ICGC outcome discrepancy between sources'
months_delta=(cohort.futime/30-cohort.SUR).abs()
assert months_delta.le(1).all()
cohort['time']=cohort.futime
cohort['age']=pd.to_numeric(cohort.AGE,errors='coerce')
cohort['sex']=cohort.GENDER.str.lower().map({'male':0,'female':1})
# HCCDB18 supplies primary tumor T category for a separate sensitivity analysis.
cohort['t_category']=pd.to_numeric(cohort.TNM_STAGE_T,errors='coerce').where(lambda x:x.isin([1,2,3,4]))
cohort['tissue']='HCC'
cohort=cohort.dropna(subset=['time','event','age','sex'])
assert (cohort.time>0).all()
cohort[['sample_id','tissue','time','event','age','sex','SAMPLE_NAME','t_category']].to_csv(a.out,index_label='patient_id')
summary={'n_matrix_samples':len(dat),'n_matrix_donors':dat.patient_id.nunique(),
         'n_t_tagged':int(dat.t_suffix.sum()),'n_t_donors':dat.loc[dat.t_suffix,'patient_id'].nunique(),
         'n_after_single_t':int((dat.t_suffix & dat.t_count.eq(1)).sum()),'n_analysis':len(cohort),
         'events':int(cohort.event.sum()),'score_label_corroboration':agreement,
         'attrition':{'one_legacy_tumor':len(one_tumor),
                      'one_hccdb_hcc_for_same_donor':len(by_annotation),
                      'linked_hccdb_clinical':len(by_clinical),
                      'linked_legacy_followup':len(months_delta),
                      'evaluable_with_clinical':len(cohort)},
         'endpoint_agreement':{'event_agreement_n':len(months_delta),
                               'exact_under_30_day_conversion':int(months_delta.eq(0).sum()),
                               'within_one_month':int(months_delta.le(1).sum())-int(months_delta.eq(0).sum())}}

print(json.dumps(summary,indent=2))
a.out.with_suffix('.json').write_text(json.dumps(summary,indent=2))
