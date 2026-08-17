"""Protein composition features + homology graph (external-DB-free)."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, itertools, pickle

AA='ACDEFGHIKLMNPQRSTVWY'

seqs={}
g=None
for line in open(D+'prot_longest.faa'):
    if line[0]=='>': g=line[1:].strip(); seqs[g]=[]
    else: seqs[g].append(line.strip())
seqs={k:''.join(v) for k,v in seqs.items()}

# Kyte-Doolittle hydropathy / charge / etc.
KD=dict(zip('AVLIPFMWGSTCYNQDEKRH',[1.8,4.2,3.8,4.5,-1.6,2.8,1.9,-0.9,-0.4,-0.8,-0.7,2.5,-1.3,-3.5,-3.5,-3.5,-3.5,-3.9,-4.5,-3.2]))
pos_aa,neg_aa='KRH','DE'
rows={}
for g,s in seqs.items():
    n=len(s)
    if n==0: continue
    comp={f'aa_{a}': s.count(a)/n for a in AA}
    comp['prot_len']=n
    comp['log_prot_len']=np.log10(n+1)
    comp['hydropathy']=np.mean([KD.get(c,0) for c in s])
    comp['frac_pos']=sum(s.count(a) for a in pos_aa)/n
    comp['frac_neg']=sum(s.count(a) for a in neg_aa)/n
    comp['net_charge']=(comp['frac_pos']-comp['frac_neg'])
    comp['frac_aromatic']=sum(s.count(a) for a in 'FWY')/n
    comp['frac_disorder']=sum(s.count(a) for a in 'PESQKA')/n   # disorder-promoting
    comp['low_complexity']=1-len(set(s))/20
    rows[g]=comp
S=pd.DataFrame.from_dict(rows,orient='index')
S.index.name='gene_id'
S.to_csv(O+'seq_features.csv.gz')
print('seq features:',S.shape)

# homology graph: gene -> {hit: bitscore}
h=pd.read_csv(D+'self_hits.tsv',sep='\t',header=None,
              names=['q','s','pid','len','ev','bits'])
h=h[h.q!=h.s]
hom={}
for q,s,b in zip(h.q,h.s,h.bits):
    hom.setdefault(q,{})[s]=max(hom.get(q,{}).get(s,0),b)
pickle.dump(hom,open(O+'homology.pkl','wb'))
print('genes with >=1 homolog:',len(hom),
      '| median n_homologs:',int(np.median([len(v) for v in hom.values()])))
