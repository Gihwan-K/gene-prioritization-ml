"""Two operating modes: retrieval (with homology) vs discovery (sequence+expression)."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, warnings
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
import scipy.sparse as sp, scipy.sparse.csgraph as csg
warnings.filterwarnings('ignore')
exec(open(os.path.join(PROJ,"src","_common.py")).read())

r,c,v=[],[],[]
for q,d in hom.items():
    if q not in Bi: continue
    for s,b in d.items():
        if s in Bi and b>=80: r.append(Bi[q]);c.append(Bi[s]);v.append(1)
ncc,cc=csg.connected_components(sp.coo_matrix((v,(r,c)),shape=(len(B),len(B))),directed=False)
groups=pd.Series(cc,index=B)
lab=y.dropna(); yl=lab.values.astype(int); idx=lab.index

def score(use_hom, cv, gr, seed_list=(0,1,2)):
    """average model score over repeated CV; returns full-universe score + metrics"""
    sc=np.zeros(len(B)); ap=[];au=[]; n=0
    for sd in seed_list:
        C = cv(5,shuffle=True,random_state=sd)
        for tr,te in C.split(idx,yl,gr):
            ff=fold_feats(idx[tr][yl[tr]==1])
            Xf=X.join(ff) if use_hom else X.join(ff[['seed_coexpr']])
            Xt=Xf.loc[idx[tr]]; med=Xt.median(); Xt=Xt.fillna(med)
            m=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.06,max_leaf_nodes=15,
                l2_regularization=1.0,random_state=0,class_weight='balanced').fit(Xt,yl[tr])
            p=m.predict_proba(Xf.loc[idx[te]].fillna(med))[:,1]
            ap.append(average_precision_score(yl[te],p)); au.append(roc_auc_score(yl[te],p))
            sc+=m.predict_proba(Xf.fillna(med))[:,1]; n+=1
    return pd.Series(sc/n,index=B), np.mean(ap),np.std(ap),np.mean(au)

print("random baseline PR = %.3f\n"%yl.mean())
s_ret,a1,s1,r1 = score(True , StratifiedKFold      , None)
print(f"RETRIEVAL mode (homology on, standard CV) : PR {a1:.3f}±{s1:.3f}  ROC {r1:.3f}  lift {a1/yl.mean():.1f}x")
s_dis,a2,s2,r2 = score(False, StratifiedGroupKFold , groups.loc[idx].values)
print(f"DISCOVERY mode (homology off, grouped CV) : PR {a2:.3f}±{s2:.3f}  ROC {r2:.3f}  lift {a2/yl.mean():.1f}x")

F_=pd.read_csv(O+'feature_matrix_v1.csv.gz',index_col=0)
e=F_.loc[B,['elim_vs_noelim','AR_vs_Root','anther_vs_leaf','tau']]
elim=(e.rank(pct=True)*[0.4,0.4,0.1,0.1]).sum(1)
res=pd.DataFrame({'func_retrieval':s_ret,'func_discovery':s_dis,'elim_score':elim})
res['func_disc_pct']=res.func_discovery.rank(pct=True)
res['func_ret_pct'] =res.func_retrieval.rank(pct=True)
res['combined']=res.func_disc_pct*res.elim_score
res['rank']=res.combined.rank(ascending=False).astype(int)
res.sort_values('rank').to_csv(O+'B_ranking_final.csv')

def rep(col,label):
    rk={g:i+1 for i,g in enumerate(res.sort_values(col,ascending=False).index)}
    out=[f"{label:32s}"]
    for nm,g in [('SYN2-B',SYN2),('CENH3-B',CENH3)]:
        out.append(f"{nm} top{100*rk[g]/len(B):5.2f}%")
    for nm,s in [('T10',t10-pos),('T6',t6-pos)]:
        s=(set(s)&set(B))-set(held); rr=np.array([rk[g] for g in s])
        hit=(rr<=200).sum(); exp=len(s)*200/len(B)
        out.append(f"{nm}@200 {hit}/{len(s)} ({hit/max(exp,1e-9):.1f}x)")
    print('  '+'   '.join(out))
print("\n=== HELD-OUT RECOVERY (SYN2-B / CENH3-B never trained on) ===")
for c,l in [('func_retrieval','functional (retrieval)'),('func_discovery','functional (discovery)'),
            ('elim_score','elimination expression'),('combined','COMBINED')]: rep(c,l)
