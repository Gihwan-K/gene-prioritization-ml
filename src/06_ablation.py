"""Feature-group ablation under both CV schemes."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, warnings
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
import scipy.sparse as sp, scipy.sparse.csgraph as csg
warnings.filterwarnings('ignore')
exec(open(os.path.join(PROJ,"src","_common.py")).read())

r,c,v=[],[],[]
for q,d in hom.items():
    if q not in Bi: continue
    for s,b in d.items():
        if s in Bi and b>=80: r.append(Bi[q]); c.append(Bi[s]); v.append(1)
ncc,cc=csg.connected_components(sp.coo_matrix((v,(r,c)),shape=(len(B),len(B))),directed=False)
groups=pd.Series(cc,index=B)

lab=y.dropna(); yl=lab.values.astype(int); idx=lab.index
AA=[c for c in X.columns if c.startswith('aa_')]
PHYS=['hydropathy','frac_pos','frac_neg','net_charge','frac_aromatic','frac_disorder',
      'low_complexity','prot_len','log_prot_len','has_protein']
STRUCT=['gene_span','exon_len','n_transcripts','n_exon_rec','exons_per_tx','log_span','log_exon_len']
EXPR=[c for c in X.columns if c not in AA+PHYS+STRUCT]
GROUPS={'expression':EXPR,'structure':STRUCT,'physicochem':PHYS,'aa_composition':AA}
FOLD=['seed_coexpr','hom_max_pos','hom_sum_pos','hom_n_pos']

def evaluate(cols, use_fold, cv, gr):
    ap=[];au=[]
    for tr,te in cv.split(idx,yl,gr):
        ff=fold_feats(idx[tr][yl[tr]==1])
        Xf=X[cols].join(ff[FOLD]) if use_fold else X[cols]
        Xt=Xf.loc[idx[tr]]; med=Xt.median(); Xt=Xt.fillna(med)
        Xe=Xf.loc[idx[te]].fillna(med)
        m=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.06,max_leaf_nodes=15,
            l2_regularization=1.0,random_state=0,class_weight='balanced').fit(Xt,yl[tr])
        p=m.predict_proba(Xe)[:,1]
        ap.append(average_precision_score(yl[te],p)); au.append(roc_auc_score(yl[te],p))
    return np.mean(ap),np.std(ap),np.mean(au)

configs=[('expression only',EXPR,False),
         ('+ gene structure',EXPR+STRUCT,False),
         ('+ physicochemical',EXPR+STRUCT+PHYS,False),
         ('+ aa composition',EXPR+STRUCT+PHYS+AA,False),
         ('+ co-expr & homology (full)',EXPR+STRUCT+PHYS+AA,True),
         ('homology+coexpr only',[],True)]
rows=[]
for tag,cols,uf in configs:
    a1=evaluate(cols,uf,StratifiedKFold(5,shuffle=True,random_state=0),None)
    a2=evaluate(cols,uf,StratifiedGroupKFold(5,shuffle=True,random_state=0),groups.loc[idx].values)
    rows.append(dict(config=tag,n_feat=len(cols)+(len(FOLD) if uf else 0),
                     pr_std=a1[0],pr_std_sd=a1[1],roc_std=a1[2],
                     pr_grp=a2[0],pr_grp_sd=a2[1],roc_grp=a2[2]))
    print(f"{tag:30s} n={rows[-1]['n_feat']:3d}  standardPR={a1[0]:.3f}±{a1[1]:.3f}  groupedPR={a2[0]:.3f}±{a2[1]:.3f}")
ab=pd.DataFrame(rows); ab.to_csv(O+'ablation.csv',index=False)
print(f"\nrandom baseline PR = {yl.mean():.3f}")
