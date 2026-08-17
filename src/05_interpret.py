"""Final model fit, SHAP interpretation, ranked candidate table."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, warnings, shap
from sklearn.ensemble import HistGradientBoostingClassifier
warnings.filterwarnings('ignore')
exec(open(os.path.join(PROJ,"src","_common.py")).read())

lab=y.dropna(); yl=lab.values.astype(int); idx=lab.index
ff=fold_feats(idx[yl==1])
Xf=X.join(ff); med=Xf.loc[idx].median()
Xtr=Xf.loc[idx].fillna(med); Xall=Xf.fillna(med)
m=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.06,max_leaf_nodes=15,
    l2_regularization=1.0,random_state=0,class_weight='balanced').fit(Xtr,yl)

ex=shap.TreeExplainer(m)
sv=ex.shap_values(Xtr)
sv=np.array(sv); 
if sv.ndim==3: sv=sv[:,:,1] if sv.shape[2]==2 else sv[1]
imp=pd.Series(np.abs(sv).mean(0),index=Xtr.columns).sort_values(ascending=False)
imp.to_csv(O+'shap_importance.csv')
print("=== TOP 15 FEATURES (mean |SHAP|) ===")
for k,v in imp.head(15).items(): print(f"  {k:22s} {v:.4f}")

svall=np.array(ex.shap_values(Xall))
if svall.ndim==3: svall=svall[:,:,1] if svall.shape[2]==2 else svall[1]
SV=pd.DataFrame(svall,index=Xall.index,columns=Xall.columns)
SV.to_csv(O+'shap_all.csv.gz'); Xtr.to_csv(O+'shap_X.csv.gz')
np.save(O+'shap_values.npy',sv)
print("\n=== WHY THESE GENES? ===")
for nm,g in [('SYN2-B',SYN2),('CENH3-B',CENH3)]:
    s=SV.loc[g].sort_values(key=abs,ascending=False).head(6)
    print(f"  {nm}:")
    for k,v in s.items(): print(f"     {k:22s} {v:+.3f}  (value={Xall.loc[g,k]:.2f})")

res=pd.read_csv(O+'B_ranking_v2.csv',index_col=0)
t3=pd.read_excel(D+'suppl.xlsx','Supplementary Table 3')
t3.columns=['tx','loc','txlen','orflen','homolog','desc']
t3['gene_id']=t3.tx.astype(str).str.replace(r'\.\d+$','',regex=True)
ann=t3.assign(short=t3.desc.astype(str).str.split(';').str[0]).groupby('gene_id')['short'].first()
top=res.sort_values('combined',ascending=False).head(50).copy()
top['annotation']=ann.reindex(top.index).fillna('not annotated')
top['is_validated']=[g in (SYN2,CENH3) for g in top.index]
top['in_T6_final']=[g in t6 for g in top.index]
top['in_T10_sorghum']=[g in t10 for g in top.index]
top[['rank','combined','func_grp_pct','elim_score','annotation',
     'is_validated','in_T6_final','in_T10_sorghum']].to_csv(O+'top50_candidates.csv')
print("\n=== TOP 15 COMBINED-SCORE B GENES ===")
for g,r in top.head(15).iterrows():
    fl=('★' if r.is_validated else '')+('⁶' if r.in_T6_final else '')+('¹⁰' if r.in_T10_sorghum else '')
    print(f"  {int(r['rank']):3d}. {g:26s} {str(r.annotation)[:50]:50s} {fl}")
