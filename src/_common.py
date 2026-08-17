
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, pickle

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
L=np.log2(tpm.loc[B]+1); Lc=L.sub(L.mean(1),axis=0)
Ln=(Lc.div(np.sqrt((Lc**2).sum(1)).replace(0,np.nan),axis=0).fillna(0)).values
def fold_feats(train_pos):
    i=[Bi[g] for g in train_pos]; s=Ln[i].mean(0); s=s/(np.linalg.norm(s)+1e-9)
    co=pd.Series(Ln@s,index=B,name='seed_coexpr')
    tp=set(train_pos); mx=[];sm=[];nh=[]
    for g in B:
        bs=[b for s2,b in hom.get(g,{}).items() if s2 in tp]
        mx.append(max(bs) if bs else 0.0); sm.append(sum(bs)); nh.append(len(bs))
    return pd.concat([co,pd.Series(mx,index=B,name='hom_max_pos'),
        pd.Series(np.log1p(sm),index=B,name='hom_sum_pos'),
        pd.Series(nh,index=B,name='hom_n_pos')],axis=1)
