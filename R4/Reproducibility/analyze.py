"""FM-AP reanalysis on patient-level HCC samples."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
from scipy.linalg import qr
from scipy.special import logit
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

SEED = 20261123
HORIZON = 1095
PENALTIES = (0.01, 0.1, 1., 10.)
CLIN = ['age', 'sex', 'stage', 'grade']
ROOT = Path(__file__).parent

parser = argparse.ArgumentParser()
parser.add_argument('project', nargs='?', type=Path, help='Project root with original source inputs')
parser.add_argument('--clinical', type=Path, help='Original TCGA LIHC clinical matrix workbook')
parser.add_argument('--purity', type=Path, help='Original GDC PanCanAtlas ABSOLUTE purity/ploidy TSV')
parser.add_argument('--mutation-only', action='store_true', help='Exit after producing mutation-association outputs')
parser.add_argument('--external-expression', type=Path, help='TSV: first column gene symbol, other columns verified specimen IDs')
parser.add_argument('--external-manifest', type=Path, help='CSV: sample_id,patient_id,tissue,time,event,age,sex and optional t_category; days, 0/1 event, female=1')
parser.add_argument('--snapshot', type=Path, help='Recompute from supplied analysis-input snapshot')
parser.add_argument('--write-snapshot', type=Path, help='Write minimized analysis-input snapshot from source inputs')
parser.add_argument('--out', type=Path, default=ROOT / 'results')
args = parser.parse_args()
out = args.out
out.mkdir(parents=True, exist_ok=True)
if args.snapshot:
    if any((args.project,args.clinical,args.purity,args.external_expression,args.external_manifest,args.write_snapshot)):
        parser.error('Snapshot replay takes its inputs exclusively from --snapshot')
    definition = pd.read_csv(args.snapshot/'module_sets.csv')
    mods = {name:set(genes.split(';')) for name,genes in zip(definition.module,definition.genes)}
else:
    if not args.project or not args.clinical:
        parser.error('Supply the project root and --clinical, or use --snapshot')
    if bool(args.external_expression) != bool(args.external_manifest):
        parser.error('External expression and manifest must be supplied together')
    data = args.project/'data'
    signature = pd.read_csv(data/'signature.csv')
    mods = {name: set(x.strip().upper() for x in genes.split(',') if x.strip())
            for name, genes in zip(signature['基因集作图名称'], signature['基因名称'])}
assert len(mods) == 63
category = {m: 'Prognostic signature' if i < 5 else 'KEGG pathway' if i < 45 else 'Immune signature'
            for i, m in enumerate(mods)}
public = list(mods)[5:]
def collapse_expression(expr):
    expr.index = expr.index.astype(str).str.upper().str.strip()
    return expr.groupby(level=0).mean().dropna()

def rank_expression(expr):
    u = rankdata(expr.to_numpy(), axis=0) / (len(expr) + 1) - .5
    return pd.DataFrame(100 * u, index=expr.index, columns=expr.columns)

def scoring(ranks, signatures):
    return pd.DataFrame({m: ranks.loc[sorted(gs.intersection(ranks.index))].mean(axis=0)
                         for m, gs in signatures.items()})

def absolute_purity(path, manifest):
    calls = pd.read_csv(path,sep='\t',usecols=['array','sample','call status','purity'])
    calls = calls.rename(columns={'sample':'dna_aliquot','call status':'call_status'})
    calls['dna_aliquot'] = calls.dna_aliquot.str.strip().str.upper()
    # Positions 0–18 identify patient, sample type, vial and tissue portion.
    # The DNA and RNA aliquots have different analyte identifiers after this key.
    calls = calls.loc[calls.dna_aliquot.str.fullmatch(
        r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z]-[0-9]{2}D-[A-Z0-9]{4}-[A-Z0-9]{2}'
    )].copy()
    calls['portion'] = calls.dna_aliquot.str[:19]
    calls['sample_code_matches'] = calls.array.str.upper().eq(calls.dna_aliquot.str[:15])
    rna = manifest.assign(portion=manifest.sample_id.str[:19])
    assert rna.portion.is_unique and rna.sample_id.str.fullmatch(
        r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[A-Z0-9]{2}'
    ).all()
    matched = rna.merge(calls,on='portion',how='left',validate='one_to_many')
    assert matched.loc[matched.dna_aliquot.notna(),'sample_code_matches'].all()
    matched['valid_purity'] = (matched.call_status.eq('called') &
                              matched.purity.between(0,1,inclusive='both'))
    valid = matched.loc[matched.valid_purity,['patient_id','purity']].drop_duplicates()
    conflicts = valid.patient_id.duplicated(keep=False)
    clean = valid.loc[~conflicts].set_index('patient_id',verify_integrity=True).purity
    dna_matched = matched.loc[matched.dna_aliquot.notna(),'patient_id'].nunique()
    audit = {
        'rna_primary_tumors':len(rna),
        'dna_portion_matches':int(dna_matched),
        'unmatched_portions':int(len(rna)-dna_matched),
        'usable_purity':len(clean),
        'conflicting_purity_patients':int(valid.loc[conflicts,'patient_id'].nunique()),
    }
    return clean,audit


def read_matrix(path):
    with np.load(path,allow_pickle=False) as f:
        return pd.DataFrame(f['values'],index=f['index'],columns=f['columns'])

def write_matrix(path, frame):
    np.savez_compressed(path,values=frame.to_numpy(),index=frame.index.to_numpy(dtype=str),
                        columns=frame.columns.to_numpy(dtype=str))

if args.snapshot:
    cache = args.snapshot
    metadata = json.loads((cache/'metadata.json').read_text())
    baseline = pd.read_csv(cache/'tcga_baseline.csv',index_col=0,float_precision='round_trip')
    mut = pd.read_csv(cache/'tcga_mutation.csv',index_col=0,float_precision='round_trip')
    ranks = read_matrix(cache/'tcga_gene_ranks.npz')
    old = read_matrix(cache/'tcga_gsva.npz')
    assert baseline.index.equals(mut.index) and mut.index.equals(old.index) and old.index.equals(ranks.columns)
    purity_audit = metadata['purity_audit']
    tcga_t = pd.read_csv(cache/'tcga_t_category.csv',index_col=0,float_precision='round_trip')['t_category']
    if metadata['has_external']:
        common = (cache/'transfer_gene_background.txt').read_text().splitlines()
        ts = read_matrix(cache/'tcga_transfer_scores.npz')
        es = read_matrix(cache/'icgc_transfer_scores.npz')
        ext_meta = pd.read_csv(cache/'icgc_baseline.csv',index_col=0,float_precision='round_trip')
        assert ts.index.equals(baseline.index) and es.index.equals(ext_meta.index)
else:
    clin = pd.read_excel(args.clinical, header=2)
    clin = clin.loc[clin['sample_type'].eq('Primary Tumor')].copy()
    clin['patient_id'] = clin['sampleID'].str[:12]
    clin = clin.set_index('patient_id', verify_integrity=True)
    clin['stage'] = clin.pathologic_stage.str.extract(r'^Stage (IV|III|II|I)(?:[ABC]?)(?:$|\s)')[0].map({'I':1,'II':2,'III':3,'IV':4})
    clin['grade'] = clin.neoplasm_histologic_grade.str.extract(r'^G([1-4])$')[0].astype(float)
    clin['age'] = pd.to_numeric(clin.age_at_initial_pathologic_diagnosis, errors='coerce')
    clin['sex'] = clin.gender.map({'MALE':0,'FEMALE':1})
    risk = clin.hist_hepato_carc_fact
    for col, token in [('HBV','Hepatitis B'),('HCV','Hepatitis C'),('alcohol','Alcohol consumption')]:
        clin[col] = risk.str.contains(token,regex=False).astype(float)
    clin['NAFLD'] = risk.str.contains('Non-Alcoholic Fatty Liver Disease',regex=False).astype(float)
    clin['other_risk_history'] = risk.str.contains(
        r'(?:^|[|])(?:Other|Hemochromatosis|Alpha-1 Antitrypsin Deficiency)(?:[|]|$)',regex=True
    ).astype(float)
    for gene in ['TP53','CTNNB1']:
        maf = pd.read_csv(data / f'alterations_across_samples_{gene}.csv')
        clin[gene] = maf.set_index('Patient ID',verify_integrity=True)['Altered']
    source_os = pd.read_csv(data/'hcc_tcga_gdc_clinical_data.tsv',sep='\t')
    source_os = source_os.loc[source_os['Sample Type'].eq('Primary Tumor') &
                              source_os['Cancer Type Detailed'].eq('Hepatocellular Carcinoma')]
    surv = source_os.set_index('Patient ID',verify_integrity=True)
    surv['time'] = np.rint(surv['Overall Survival (Months)'] * 30.44)
    surv['event'] = surv['Overall Survival Status'].map({'0:LIVING':0,'1:DECEASED':1})
    surv = surv.loc[surv.time.gt(0), ['time','event']]
    old = pd.read_csv(args.project / 'results/gsva_scores.csv',index_col=0).T
    expr = pd.read_csv(data / 'TCGA-LIHC_expression_data.txt',sep='\t',index_col=0)
    # Legacy list provides an auditable one-primary-tumor-per-patient sample manifest.
    assert all(x[13:15]=='01' for x in old.index) and not old.index.str[:12].duplicated().any()
    assert old.index.isin(expr.columns).all()
    manifest = pd.DataFrame({'sample_id':old.index,'patient_id':old.index.str[:12]})
    if args.purity:
        clin['purity'], purity_audit = absolute_purity(args.purity,manifest)
    expr = collapse_expression(expr.loc[:,old.index])
    expr.columns = old.index.str[:12]
    old.index = old.index.str[:12]
    ranks = rank_expression(expr)
    baseline = clin[CLIN].join(surv).reindex(ranks.columns)
    mut = clin[['stage','grade','age','sex','HBV','HCV','alcohol','NAFLD','other_risk_history','TP53','CTNNB1'] + (['purity'] if args.purity else [])].reindex(ranks.columns)
    manifest = manifest.join(baseline[['time','event']],on='patient_id')
    manifest['in_prediction'] = manifest.patient_id.isin(baseline.dropna().index)
    manifest['in_mutation'] = manifest.patient_id.isin(mut.drop(columns=['purity'],errors='ignore').dropna().index)
    manifest.to_csv(out/'tcga_manifest.csv',index=False)
    tcga_t = source_os.set_index('Patient ID')['AJCC Pathologic T-Stage'].str.extract(
        r'^T([1-4])(?:[a-z])?$',expand=False).astype(float)
    if args.external_manifest:
        m = pd.read_csv(args.external_manifest)
        assert m.tissue.eq('HCC').all() and m.patient_id.is_unique and m.sample_id.is_unique
        assert m.event.isin([0,1]).all() and (m.time>0).all()
        assert m[['age','sex']].notna().all().all()
        ext = pd.read_csv(args.external_expression,sep='\t',index_col=0, low_memory=False).apply(pd.to_numeric,errors='coerce')
        assert m.sample_id.isin(ext.columns).all()
        ext = collapse_expression(ext[m.sample_id])
        ext.columns = m.patient_id
        common = sorted(set(expr.index)&set(ext.index))
        es = scoring(rank_expression(ext.loc[common]),mods)
        ts = scoring(rank_expression(expr.loc[common]),mods)
        ext_meta = m.set_index('patient_id')
    if args.write_snapshot:
        cache = args.write_snapshot
        cache.mkdir(parents=True,exist_ok=False)
        # Preserve row order for CV and lexical order for patient-row permutations.
        tcga_ids = {bc:f'TCGA_{i:04d}' for i,bc in enumerate(sorted(baseline.index),1)}
        pd.DataFrame([(name,';'.join(sorted(genes))) for name,genes in mods.items()],
                     columns=['module','genes']).to_csv(cache/'module_sets.csv',index=False)
        for name, tab in [('tcga_baseline',baseline),('tcga_mutation',mut),
                          ('tcga_t_category',tcga_t.reindex(baseline.index).rename('t_category').to_frame())]:
            tab.rename(index=tcga_ids).to_csv(cache/f'{name}.csv',index_label='record')
        keep = sorted(set.union(*mods.values()) & set(ranks.index))
        write_matrix(cache/'tcga_gene_ranks.npz',ranks.loc[keep].rename(columns=tcga_ids))
        write_matrix(cache/'tcga_gsva.npz',old.rename(index=tcga_ids))
        if args.external_manifest:
            icgc_ids = {bc:f'ICGC_{i:04d}' for i,bc in enumerate(ext_meta.index,1)}
            variables = [v for v in ['time','event','age','sex','t_category'] if v in ext_meta]
            ext_meta[variables].rename(index=icgc_ids).to_csv(cache/'icgc_baseline.csv',index_label='record')
            write_matrix(cache/'tcga_transfer_scores.npz',ts.rename(index=tcga_ids))
            write_matrix(cache/'icgc_transfer_scores.npz',es.rename(index=icgc_ids))
            (cache/'transfer_gene_background.txt').write_text('\n'.join(common)+'\n')
        (cache/'metadata.json').write_text(json.dumps(
            {'has_external':bool(args.external_manifest),
             'purity_audit':purity_audit if args.purity else None},indent=2))
scores = scoring(ranks,mods)
predict = baseline.dropna().join(scores)
summary = {'n_scored':len(scores),'n_os':int(baseline.time.notna().sum()),
           'n_after_os_and_clinical':len(predict),'n_events':int(predict.event.sum()),
           'os_events':int(baseline.event.sum()),'missing_baseline':baseline.isna().sum().to_dict(),
           'module_size_min':min(len(gs) for gs in mods.values())}

def adjusted(scores_in, cov, permutations, eligible=None):
    data_rows = mut if eligible is None else mut.loc[eligible]
    frame = data_rows[cov+['TP53','CTNNB1']].dropna().sort_index()
    Y = scores_in.loc[frame.index,list(mods)].to_numpy()
    x = pd.get_dummies(frame[['stage','grade']].astype(int).astype(str),drop_first=True,dtype=float)
    x = pd.concat([x,frame[[c for c in cov if c not in ('stage','grade')]]],axis=1)
    W = sm.add_constant(x,has_constant='add').to_numpy(dtype=float)
    D = frame[['TP53','CTNNB1']].to_numpy(dtype=float)
    X = np.column_stack([W,D])
    assert np.linalg.matrix_rank(X) == X.shape[1]
    rows = []
    for i,m in enumerate(mods):
        fit = sm.OLS(Y[:,i],X).fit(cov_type='HC3',use_t=True)
        for j,g in enumerate(['TP53','CTNNB1']):
            k = X.shape[1]-2+j
            rows.append((g,m,category[m],fit.params[k],*fit.conf_int()[k],fit.pvalues[k]))
    tab = pd.DataFrame(rows,columns=['mutation','module','family','effect','low','high','p'])
    tab['q'] = multipletests(tab.p,method='fdr_bh')[1]
    pivot = tab.pivot(index='module',columns='mutation',values='effect')
    observed = float(pivot.corr().iloc[0,1])
    Q = qr(W,mode='economic')[0]
    # einsum avoids spurious floating-point kernel warnings from Accelerate matmul.
    res = Y-np.einsum("ij,jk->ik",Q,np.einsum("ji,jk->ik",Q,Y))
    Dr = D-np.einsum("ij,jk->ik",Q,np.einsum("ji,jk->ik",Q,D))
    pinv = np.linalg.pinv(Dr)
    generated = np.random.default_rng(SEED)
    large = 0
    for _ in range(permutations):
        coef = np.einsum("ij,jk->ik",pinv,res[generated.permutation(len(res))])
        assert np.isfinite(coef).all()
        r = np.corrcoef(coef)[0,1]
        large += abs(r) >= abs(observed)
    meta = {'n':len(frame),'tp53':int(frame.TP53.sum()),'ctnnb1':int(frame.CTNNB1.sum()),
            'joint':int(((frame.TP53==1)&(frame.CTNNB1==1)).sum()),'r':observed,
            'perm_p':(1+large)/(1+permutations),'permutations':permutations,
            'discoveries':tab[tab.q<.05].groupby('mutation').size().to_dict()}
    return tab,meta

cov = ['age','sex','stage','grade','HBV','HCV','alcohol']
a, summary['mutation'] = adjusted(scores,cov,9999)
a.to_csv(out/'mutation_adjusted.csv',index=False)
g, summary['mutation_gsva'] = adjusted(old,cov,9999)
g.to_csv(out/'mutation_gsva_sensitivity.csv',index=False)
summary['covariates'] = cov
extended_etiology = cov + ['NAFLD','other_risk_history']
for name, score_matrix in [('rank',scores),('gsva',old)]:
    estimates, summary[f'mutation_etiology_{name}'] = adjusted(score_matrix,extended_etiology,9999)
    estimates.to_csv(out/f'mutation_etiology_{name}.csv',index=False)
if 'purity' in mut:
    eligible = mut[extended_etiology+['purity','TP53','CTNNB1']].dropna().index
    if not len(eligible):
        raise ValueError('No complete TCGA mutation cases matched tumor purity by primary-tumor tissue portion')
    complete = mut[extended_etiology+['TP53','CTNNB1']].dropna()
    matched_purity = complete.join(mut.purity).loc[eligible,'purity']
    summary['purity_audit'] = dict(purity_audit, clinical_mutation_complete=len(complete),
                                   matched_mutation_cases=len(eligible),
                                   purity_median=float(matched_purity.median()),
                                   purity_q1=float(matched_purity.quantile(.25)),
                                   purity_q3=float(matched_purity.quantile(.75)),
                                   purity_min=float(matched_purity.min()),
                                   purity_max=float(matched_purity.max()))
    paired_discoveries = []
    for name, score_matrix in [('rank',scores),('gsva',old)]:
        pair = {}
        for label,covars in [('matched_clinical',extended_etiology),
                             ('matched_purity',extended_etiology+['purity'])]:
            result_key = f'mutation_{label}_{name}'
            estimates, summary[result_key] = adjusted(score_matrix,covars,9999,eligible=eligible)
            estimates.to_csv(out/f'{result_key}.csv',index=False)
            pair[label] = estimates.set_index(['mutation','module'])
        for gene in ('TP53','CTNNB1'):
            clinical = pair['matched_clinical'].loc[gene,'q'].lt(.05)
            with_purity = pair['matched_purity'].loc[gene,'q'].lt(.05)
            paired_discoveries.append((name,gene,len(eligible),int(clinical.sum()),int(with_purity.sum()),
                                       int((clinical & with_purity).sum())))
    pd.DataFrame(paired_discoveries,columns=['scoring','mutation','matched_n',
                                            'clinical_q_lt_0_05','purity_q_lt_0_05','both_q_lt_0_05']).to_csv(
        out/'mutation_purity_discovery_comparison.csv',index=False)
    pd.DataFrame([summary['purity_audit']]).to_csv(out/'tcga_purity_audit.csv',index=False)
if args.mutation_only:
    (out/'mutation_summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item()))
    raise SystemExit

rho,pvalues = spearmanr(scores,axis=0)
i,j = np.triu_indices(len(mods),1)
qvals = multipletests(pvalues[i,j],method='fdr_bh')[1]
network = pd.DataFrame({'module_1':np.array(list(mods))[i],'module_2':np.array(list(mods))[j],
                        'rho':rho[i,j],'p':pvalues[i,j],'q':qvals})
adjrho = []
overlap = []
for a,b,r in zip(network.module_1,network.module_2,network.rho):
    ga,gb = mods[a]&set(ranks.index),mods[b]&set(ranks.index)
    shared = ga&gb
    overlap.append(len(shared))
    if shared and ga-shared and gb-shared:
        xa = ranks.loc[sorted(ga-shared)].mean(axis=0)
        xb = ranks.loc[sorted(gb-shared)].mean(axis=0)
        adjrho.append(float(spearmanr(xa,xb).statistic))
    else:
        adjrho.append(np.nan if shared else r)
network['overlap'] = overlap
network['rho_without_shared_genes'] = adjrho
edges = network.query('q < .05 and abs(rho) >= .3')
network.to_csv(out/'network_all_pairs.csv',index=False)
summary['network'] = {'total_pairs':len(network),'display_edges':len(edges),
                      'positive':int((edges.rho>0).sum()),'negative':int((edges.rho<0).sum()),
                      'shared_edges':int((edges.overlap>0).sum()),
                      'preserved':int((edges.rho_without_shared_genes.abs().ge(.3)&
                                      np.sign(edges.rho).eq(np.sign(edges.rho_without_shared_genes))).sum())}
pd.DataFrame([{'module':m,'family':category[m],'total_genes':len(gs),'retained_genes':len(gs&set(ranks.index)),
               'genes_used':';'.join(sorted(gs&set(ranks.index)))} for m,gs in mods.items()]).to_csv(out/'modules.csv',index=False)
coverage = pd.read_csv(out/'modules.csv')
assert (coverage.retained_genes>=3).all()

def fit_ridge(df, cols, penalty):
    model = CoxPHFitter(penalizer=penalty)
    return model.fit(df[['time','event']+cols],duration_col='time',event_col='event')

def risk_at(model, df, cols):
    hazard = np.exp(model.predict_log_partial_hazard(df[cols]).to_numpy())
    baseline_hazard = model.baseline_cumulative_hazard_.iloc[:,0]
    eligible = baseline_hazard.loc[baseline_hazard.index<=HORIZON]
    return -np.expm1(-(eligible.iloc[-1] if len(eligible) else 0.)*hazard)

def tune(df, cols):
    folds = list(StratifiedKFold(4,shuffle=True,random_state=SEED).split(df,df.event))
    measures = []
    for penalty in PENALTIES:
        values = []
        for ii,jj in folds:
            tr,te = df.iloc[ii],df.iloc[jj]
            model = fit_ridge(tr,cols,penalty)
            values.append(concordance_index(te.time,-model.predict_log_partial_hazard(te[cols]),te.event))
        measures.append(np.mean(values))
    return PENALTIES[int(np.argmax(measures))]

def censoring_weights(df):
    km = KaplanMeierFitter().fit(df.time,event_observed=1-df.event)
    target = (df.time<=HORIZON) & df.event.eq(1)
    control = df.time>HORIZON
    weights = np.zeros(len(df))
    weights[target] = 1/np.array(km.predict(np.nextafter(df.loc[target,'time'].to_numpy(dtype=float),-np.inf)))
    weights[control] = 1/km.predict(HORIZON)
    return target.to_numpy(dtype=int),weights

def metrics(df, risk, target, weights, groups):
    concordance, concordance_pairs, auc, auc_pairs = 0., 0, 0., 0.
    for fold in np.unique(groups):
        i = np.flatnonzero(groups == fold)
        times, event = df.time.to_numpy()[i], df.event.to_numpy()[i]
        earlier = (times[:, None] < times[None, :]) | (
            (times[:, None] == times[None, :]) & (event[None, :] == 0))
        valid = earlier & (event[:, None] == 1)
        n = int(valid.sum())
        if n:
            r = risk[i]
            concordance += n * concordance_index(times, -r, event)
            concordance_pairs += n
        weight_cases = np.sum(weights[i] * target[i])
        weight_controls = np.sum(weights[i] * (1-target[i]))
        if weight_cases and weight_controls:
            n_weight = weight_cases * weight_controls
            auc += n_weight * roc_auc_score(target[i],risk[i],sample_weight=weights[i])
            auc_pairs += n_weight
    return [concordance/concordance_pairs, auc/auc_pairs, np.mean(weights*(target-risk)**2)]

def score_table(df, profiles, groups, *, seed):
    yy,ww = censoring_weights(df)
    models = list(profiles)
    pp = np.asarray([profiles[name] for name in models])
    estimates = np.array([metrics(df,x,yy,ww,groups) for x in pp])
    rng = np.random.default_rng(seed)
    parts = [np.flatnonzero(groups==f) for f in np.unique(groups)]
    all_boot = []
    for _ in range(1500):
        idx = np.concatenate([rng.choice(part,size=len(part),replace=True) for part in parts])
        d = df.iloc[idx]
        y,w = censoring_weights(d)
        all_boot.append([metrics(d,x[idx],y,w,groups[idx]) for x in pp])
    draws = np.array(all_boot)
    rows = []
    for k, name in enumerate(models):
        for j, metric in enumerate(['C_index','AUC_3y','Brier_3y']):
            lo,hi = np.percentile(draws[:,k,j],[2.5,97.5])
            if name=='Null' and j!=2:
                continue
            rows.append((name,metric,estimates[k,j],lo,hi))
    for after,before in [('Integrated','Clinical'),('FM-AP (all 63)','FM-AP')]:
        if after not in models or before not in models:
            continue
        a,b = models.index(after),models.index(before)
        for j, metric in enumerate(['C_index','AUC_3y','Brier_3y']):
            lo,hi = np.percentile(draws[:,a,j]-draws[:,b,j],[2.5,97.5])
            rows.append((f'{after} − {before}',metric,estimates[a,j]-estimates[b,j],lo,hi))
    return pd.DataFrame(rows,columns=['model','metric','estimate','low','high'])

def calibration_coefficients(target, weights, logits):
    mod = sm.GLM(target,sm.add_constant(logits),family=sm.families.Binomial(),freq_weights=weights).fit()
    citl = sm.GLM(target,np.ones(len(target)),family=sm.families.Binomial(),
                  freq_weights=weights,offset=logits).fit()
    if not (mod.converged and citl.converged):
        raise RuntimeError('Calibration model did not converge')
    return float(citl.params[0]),float(mod.params[1])

def calibration(df, profiles, prefix, groups, *, seed):
    target,weights = censoring_weights(df)
    rows, points = [], []
    models = [name for name in profiles if name!='Null']
    logit_risks = {name:logit(np.clip(profiles[name],1e-6,1-1e-6)) for name in models}
    draws = {name:[] for name in models}
    rng = np.random.default_rng(seed)
    parts = [np.flatnonzero(groups==f) for f in np.unique(groups)]
    for _ in range(1500):
        idx = np.concatenate([rng.choice(part,size=len(part),replace=True) for part in parts])
        y,w = censoring_weights(df.iloc[idx])
        for name in models:
            draws[name].append(calibration_coefficients(y,w,logit_risks[name][idx]))
    for name in models:
        intercept,slope = calibration_coefficients(target,weights,logit_risks[name])
        lo,hi = np.percentile(draws[name],[2.5,97.5],axis=0)
        rows.append((name,intercept,slope,lo[0],hi[0],lo[1],hi[1],len(draws[name])))
        risk = profiles[name]
        bins = pd.qcut(risk,4,labels=False,duplicates='drop')
        for k in np.unique(bins):
            i = np.flatnonzero(bins==k)
            q = df.iloc[i]
            km = KaplanMeierFitter().fit(q.time,event_observed=q.event)
            points.append((name,int(k)+1,len(q),risk[i].mean(),1-km.predict(HORIZON)))
    pd.DataFrame(rows,columns=['model','intercept_fixed_slope','slope',
                              'intercept_low','intercept_high','slope_low','slope_high',
                              'bootstrap_n']).to_csv(out/f'{prefix}_calibration.csv',index=False)
    pd.DataFrame(points,columns=['model','quartile','n','predicted','observed']).to_csv(out/f'{prefix}_calibration_plot.csv',index=False)
    followup = pd.DataFrame([{'n':len(df),'horizon_days':HORIZON,
                              'deaths_by_horizon':int(sum(target)),
                              'censored_by_horizon':int(((df.time<=HORIZON)&df.event.eq(0)).sum()),
                              'observed_beyond_horizon':int(df.time.gt(HORIZON).sum())}])
    followup.to_csv(out/f'{prefix}_followup.csv',index=False)

def ph_diagnostics(df, feature_sets, prefix):
    frames = []
    for name,cols in feature_sets.items():
        train = df[['time','event']+cols]
        alpha = tune(df,cols)
        model = fit_ridge(df,cols,alpha)
        tests = proportional_hazard_test(model,train,time_transform='rank').summary
        tests = tests.rename_axis('term').reset_index().rename(columns={'test_statistic':'chi_square_1df'})
        tests['q'] = multipletests(tests.p,method='fdr_bh')[1]
        tests.insert(0,'model',name)
        tests['penalty'] = alpha
        tests['n'] = len(train)
        tests['events'] = int(train.event.sum())
        frames.append(tests)
    pd.concat(frames,ignore_index=True).to_csv(out/f'{prefix}_ph_diagnostics.csv',index=False)

def null_risk(df):
    return 1-KaplanMeierFitter().fit(df.time,event_observed=df.event).predict(HORIZON)

def predictor(df, feature_sets, prefix, external=None, external_prefix='external'):
    ids = df.index
    profiles = {k:np.zeros(len(df)) for k in [*feature_sets,'Null']}
    outer = StratifiedKFold(5,shuffle=True,random_state=SEED)
    folds = []
    groups = np.zeros(len(df),dtype=int)
    for fold,(i,j) in enumerate(outer.split(df,df.event),1):
        tr,te = df.iloc[i],df.iloc[j]
        groups[j] = fold
        profiles['Null'][j] = null_risk(tr)
        for name,cols in feature_sets.items():
            alpha = tune(tr,cols)
            model = fit_ridge(tr,cols,alpha)
            profiles[name][j] = risk_at(model,te,cols)
            c = concordance_index(te.time,-model.predict_log_partial_hazard(te[cols]),te.event)
            folds.append((fold,name,alpha,len(tr),len(te),int(te.event.sum()),float(c)))
    pd.DataFrame(folds,columns=['fold','model','penalty','n_train','n_assess','events','C_index']).to_csv(out/f'{prefix}_folds.csv',index=False)
    pd.DataFrame(profiles,index=ids).assign(fold=groups).join(df[['time','event']]).to_csv(out/f'{prefix}_predictions.csv')
    score_table(df,profiles,groups,seed=SEED).to_csv(out/f'{prefix}_metrics.csv',index=False)
    calibration(df,profiles,prefix,groups,seed=SEED)
    ph_diagnostics(df,feature_sets,prefix)
    if external is not None:
        exprofiles = {}
        final_fits = []
        for name,cols in feature_sets.items():
            alpha = tune(df,cols)
            final_fits.append((name,alpha,len(df),int(df.event.sum())))
            exprofiles[name] = risk_at(fit_ridge(df,cols,alpha),external,cols)
        exprofiles['Null']=np.repeat(null_risk(df),len(external))
        pd.DataFrame(exprofiles,index=external.index).join(external[['time','event']]).to_csv(out/f'{external_prefix}_predictions.csv')
        score_table(external,exprofiles,np.zeros(len(external),dtype=int),seed=SEED+1).to_csv(out/f'{external_prefix}_metrics.csv',index=False)
        calibration(external,exprofiles,external_prefix,np.zeros(len(external),dtype=int),seed=SEED+1)
        pd.DataFrame(final_fits,columns=['model','penalty','n_train','training_events']).to_csv(
            out/f'{external_prefix}_fit.csv',index=False)

feature_sets = {'Clinical':CLIN,'FM-AP':public,'Integrated':CLIN+public,'FM-AP (all 63)':list(mods)}
predictor(predict,feature_sets,'internal')
if (args.snapshot and metadata['has_external']) or args.external_manifest:
    assert min(len(gs&set(common)) for gs in mods.values())>=3
    pd.DataFrame([{'module':m,'retained_genes_common':len(gs&set(common)),'genes_used_common':';'.join(sorted(gs&set(common)))} for m,gs in mods.items()]).to_csv(out/'transfer_modules.csv',index=False)
    shared_base = [x for x in CLIN if x in ext_meta and ext_meta[x].notna().all()]
    assert shared_base
    preds = {'Clinical':shared_base,'FM-AP':public,'Integrated':shared_base+public}
    predictor(predict[['time','event']+CLIN].join(ts),preds,'transfer_internal',
              ext_meta[['time','event']+shared_base].join(es))
    if 't_category' in ext_meta:
        # Use primary tumor T category (not overall stage) and collapse T3/T4
        # to reduce dependence on staging edition and sparse category counts.
        dev = predict[['time','event','age','sex']].join(tcga_t.rename('t_category')).dropna().join(ts)
        val = ext_meta[['time','event','age','sex','t_category']].dropna().join(es)
        for cohort in (dev,val):
            cohort['t2'] = cohort.t_category.eq(2).astype(int)
            cohort['t3_4'] = cohort.t_category.ge(3).astype(int)
        t_base = ['age','sex','t2','t3_4']
        t_models = {'Clinical':t_base,'FM-AP':public,'Integrated':t_base+public}
        predictor(dev,t_models,'tcategory_internal',val,external_prefix='tcategory_external')
        cohort_rows = []
        for name, cohort in [('TCGA',dev),('ICGC',val)]:
            for category_name, group in cohort.groupby(cohort.t_category.clip(upper=3)):
                cohort_rows.append((name,{1:'T1',2:'T2',3:'T3-T4'}[category_name],len(group),int(group.event.sum())))
        pd.DataFrame(cohort_rows,columns=['cohort','T_category','n','deaths']).to_csv(
            out/'tcategory_cohorts.csv',index=False)
        summary['t_category'] = {'n_tcga':len(dev),'events_tcga':int(dev.event.sum()),
                                  'n_icgc':len(val),'events_icgc':int(val.event.sum())}
    targets=['CSC_Activity_PRG','KEGG_CELL_CYCLE','LipidMetabolism_Risk_LMRG',
             'KEGG_FATTY_ACID_METABOLISM','Macrophages','Autophagy_Endothelial_HighRisk_AutRG']
    assoc_rows=[]
    assoc_ph=[]
    for cohort, matrix in [('TCGA',predict[['time','event']+shared_base].join(ts)),
                            ('ICGC',ext_meta[['time','event']+shared_base].join(es))]:
        for module in targets:
            one=matrix[['time','event',*shared_base,module]].copy()
            one[module]=(one[module]-ts.loc[predict.index,module].mean())/ts.loc[predict.index,module].std()
            m=CoxPHFitter().fit(one,'time','event')
            stats=m.summary.loc[module]
            diagnostic=proportional_hazard_test(m,one,time_transform='rank').summary.loc[module]
            assoc_ph.append((cohort,module,len(one),int(one.event.sum()),diagnostic.test_statistic,diagnostic.p))
            assoc_rows.append((cohort,module,len(one),int(one.event.sum()),stats['exp(coef)'],
                               stats['exp(coef) lower 95%'],stats['exp(coef) upper 95%'],stats['p']))
    assoc=pd.DataFrame(assoc_rows,columns=['cohort','module','n','events','HR','low','high','p'])
    assoc['q']=assoc.groupby('cohort').p.transform(lambda v:multipletests(v,method='fdr_bh')[1])
    assoc.to_csv(out/'module_cox_associations.csv',index=False)
    ph_assoc=pd.DataFrame(assoc_ph,columns=['cohort','module','n','events','chi_square_1df','p'])
    ph_assoc['q']=ph_assoc.groupby('cohort').p.transform(lambda v:multipletests(v,method='fdr_bh')[1])
    ph_assoc.to_csv(out/'module_cox_ph_diagnostics.csv',index=False)
    summary['external'] = {'n':len(ext_meta),'events':int(ext_meta.event.sum()),
                           'common_genes':len(common),'clinical':shared_base}
(out/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item()))
print(json.dumps(summary,indent=2,default=lambda x:x.item()))
