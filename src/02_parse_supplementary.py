"""Parse Supplementary Tables -> GO map for B genes + curated candidate sets."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
os.makedirs(O, exist_ok=True)
import pandas as pd

xl = pd.ExcelFile(D+'suppl.xlsx')

# --- Table 3: B transcripts + AEG-9674-1 homolog description carrying GO terms
t3 = pd.read_excel(xl, 'Supplementary Table 3')
t3.columns = ['tx','loc','tx_len','orf_len','homolog','desc']
t3['gene_id'] = t3.tx.astype(str).str.replace(r'\.\d+$', '', regex=True)
t3['GO'] = t3.desc.astype(str).str.extract(r'Ontology_term=([^;]+)')
go = (t3.dropna(subset=['GO']).groupby('gene_id')['GO']
        .apply(lambda s: ','.join(sorted({x for v in s for x in v.split(',')}))))
go.reset_index().to_csv(O+'go_B.csv', index=False)
print(f"Table 3 -> GO for {len(go)} B genes "
      f"(mean {go.str.count(',').add(1).mean():.1f} terms/gene)")

# --- Table 5 / 6 / 10: curated sets
def sheet(name, cols):
    d = pd.read_excel(xl, name, skiprows=1)
    d.columns = cols
    return d.dropna(subset=[cols[0]])

t5 = sheet('Supplementary Table 5', ['gene_id','annot'])
t6 = sheet('Supplementary Table 6', ['gene_id','annot'])
t10 = sheet('Supplementary Table 10', ['sorghum','aesp_tx','annot'])
t10['gene_id'] = t10.aesp_tx.astype(str).str.replace(r'\.\d+\.p\d+$', '', regex=True)

t5.to_csv(O+'t5.csv', index=False)
t6.to_csv(O+'t6.csv', index=False)
t10.to_csv(O+'t10.csv', index=False)
print(f"Table 5  chromosome-segregation B genes : {len(t5)}")
print(f"Table 6  final candidates              : {len(t6)}")
print(f"Table 10 sorghum-conserved candidates   : {len(t10)}")
