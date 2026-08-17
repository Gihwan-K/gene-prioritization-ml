"""v2: + sequence composition, guilt-by-homology, dual CV schemes, two-axis score."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, pickle, warnings
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
warnings.filterwarnings('ignore')


F=pd.read_csv(O+'feature_matrix_v1.csv.gz',index_col=0)
S=pd.read_csv(O+'seq_features.csv.gz',index_col=0)
tpm=pd.read_csv(D+'tpm.csv.gz',index_col=0)
go=pd.read_csv(O+'go_B.csv').set_index('gene_id')['GO']
t5=set(pd.read_csv(O+'t5.csv').gene_id); t6=set(pd.read_csv(O+'t6.csv').gene_id)
t10=set(pd.read_csv(O+'t10.csv').gene_id)
hom=pickle.load(open(O+'homology.pkl','rb'))
SYN2,CENH3='Aesp_chrB_ipkv1.63465','Aesp_chrB_ipkv1.64685'

B=list(F.index[F.is_B==1]); Bi={g:i for i,g in enumerate(B)}
gob=go.reindex(B).fillna('')
SEG=['GO:0007059','GO:0000775','GO:0008278','GO:0000776','GO:0000777','GO:0007062',
     'GO:0051301','GO:0000280','GO:0140014','GO:0045132','GO:0007076','GO:0034085',
     'GO:0007080','GO:0051276','GO:0000278','GO:0005819','GO:0072686']
pos=set(gob[gob.str.contains('|'.join(SEG))].index)|(t5&set(B))
neg=(set(gob[gob!=''].index)|t5)-pos
held=[SYN2,CENH3]
y=pd.Series(np.nan,index=B); y.loc[list(pos)]=1; y.loc[list(neg)]=0; y.loc[held]=np.nan

num=[c for c in F.select_dtypes(include=[np.number]).columns if c not in('is_B','is_scaffold')]
X=F.loc[B,num].join(S,how='left').replace([np.inf,-np.inf],np.nan)
X['has_protein']=S.reindex(B).notna().any(axis=1).astype(int)

# co-expression basis
L=np.log2(tpm.loc[B]+1); Lc=L.sub(L.mean(1),axis=0)
Ln=(Lc.div(np.sqrt((Lc**2).sum(1)).replace(0,np.nan),axis=0).fillna(0)).values

# homology clusters (for conservative CV)
import scipy.sparse as sp, scipy.sparse.csgraph as csg
r,c,v=[],[],[]
for q,d in hom.items():
    if q not in Bi: continue
    for s,b in d.items():
        if s in Bi and b>=80: r.append(Bi[q]); c.append(Bi[s]); v.append(1)
A=sp.coo_matrix((v,(r,c)),shape=(len(B),len(B)))
ncc,lab_cc=csg.connected_components(A,directed=False)
groups=pd.Series(lab_cc,index=B)
print(f"universe={len(B)} pos={int((y==1).sum())} neg={int((y==0).sum())} homology-clusters={ncc}")
print(f"positives spread over {groups[list(pos)].nunique()} clusters")

def fold_feats(train_pos):
    """features that must be recomputed per fold to stay leak-free"""
    idx=[Bi[g] for g in train_pos]
    s=Ln[idx].mean(0); s/= (np.linalg.norm(s)+1e-9)
    co=pd.Series(Ln@s,index=B,name='seed_coexpr')
    tp=set(train_pos)
    mx=[];sm=[];nh=[]
    for g in B:
        d=hom.get(g,{})
        bs=[b for s2,b in d.items() if s2 in tp]
        mx.append(max(bs) if bs else 0.0); sm.append(sum(bs)); nh.append(len(bs))
    return pd.concat([co,
        pd.Series(mx,index=B,name='hom_max_pos'),
        pd.Series(np.log1p(sm),index=B,name='hom_sum_pos'),
        pd.Series(nh,index=B,name='hom_n_pos')],axis=1)

lab=y.dropna(); yl=lab.values.astype(int); idx=lab.index

def run(cv, gr, tag):
    sc=np.zeros(len(B)); ap=[];au=[]; imps=[]
    for tr,te in cv.split(idx,yl,gr):
        ff=fold_feats(idx[tr][yl[tr]==1]); Xf=X.join(ff)
        Xtr=Xf.loc[idx[tr]]; med=Xtr.median()
        Xtr=Xtr.fillna(med); Xte=Xf.loc[idx[te]].fillna(med); Xall=Xf.fillna(med)
        m=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.06,
              max_leaf_nodes=15,l2_regularization=1.0,random_state=0,
              class_weight='balanced')
        m.fit(Xtr,yl[tr])
        p=m.predict_proba(Xte)[:,1]
        ap.append(average_precision_score(yl[te],p)); au.append(roc_auc_score(yl[te],p))
        sc+=m.predict_proba(Xall)[:,1]/cv.get_n_splits()
    ap=np.array(ap);au=np.array(au)
    print(f"\n[{tag}]  PR-AUC {ap.mean():.3f}±{ap.std():.3f}   ROC-AUC {au.mean():.3f}±{au.std():.3f}   (random PR={yl.mean():.3f}, lift {ap.mean()/yl.mean():.1f}x)")
    return pd.Series(sc,index=B)

s_std = run(StratifiedKFold(5,shuffle=True,random_state=0), None, "standard 5-fold  (optimistic: paralogs may co-occur)")
s_grp = run(StratifiedGroupKFold(5,shuffle=True,random_state=0), groups.loc[idx].values,
            "homology-grouped 5-fold  (conservative: novel-family discovery)")

res=pd.DataFrame({'func_std':s_std,'func_grp':s_grp})
# axis 2: elimination-associated expression (no labels involved)
e=F.loc[B,['elim_vs_noelim','AR_vs_Root','anther_vs_leaf','tau']]
res['elim_score']=(e.rank(pct=True)*[0.4,0.4,0.1,0.1]).sum(1)
for k in ['func_std','func_grp']:
    res[k+'_pct']=res[k].rank(pct=True)
res['combined']=(res.func_grp_pct*res.elim_score)
res=res.sort_values('combined',ascending=False)
res['rank']=np.arange(1,len(res)+1)
res.to_csv(O+'B_ranking_v2.csv')

def report(col):
    r=res.sort_values(col,ascending=False); rk={g:i+1 for i,g in enumerate(r.index)}
    print(f"\n--- ranking by {col} ---")
    for nm,g in [('SYN2-B',SYN2),('CENH3-B',CENH3)]:
        print(f"   {nm:8s} rank {rk[g]:5d}/{len(B)}  top {100*rk[g]/len(B):.2f}%")
    for nm,s in [('sorghum-conserved T10',t10-pos),('final candidates T6',t6-pos)]:
        s=(set(s)&set(B))-set(held); rr=np.array([rk[g] for g in s])
        line=f"   {nm:22s}"
        for k in (100,300,600):
            hit=(rr<=k).sum(); exp=len(s)*k/len(B)
            line+=f"  top{k}:{hit}/{len(s)}({hit/max(exp,1e-9):.1f}x)"
        print(line)
for c in ['func_std','func_grp','elim_score','combined']: report(c)
