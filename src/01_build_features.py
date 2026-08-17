"""Build the gene-level feature matrix from RSEM + GTF."""
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np



cnt  = pd.read_csv(D+'counts.csv.gz', index_col=0)
tpm  = pd.read_csv(D+'tpm.csv.gz',    index_col=0)
meta = pd.read_csv(D+'meta.csv')
st   = pd.read_csv(D+'gene_struct.tsv', sep='\t', index_col=0)
ntx  = pd.read_csv(D+'gene_ntx.tsv',    sep='\t', index_col=0)
deg  = pd.read_csv(D+'deg_full.csv.gz')

genes = cnt.index
F = pd.DataFrame(index=genes)

# ---- 1. location / structure -------------------------------------------
st = st.reindex(genes); ntx = ntx.reindex(genes)
F['chrom']         = st['chrom']
F['is_B']          = (st['chrom'] == 'chrB').astype(int)
F['is_scaffold']   = st['chrom'].str.startswith('scaffold').astype(int)
F['gene_span']     = st['span']
F['exon_len']      = st['total_exon_len']
F['n_transcripts'] = ntx['n_transcripts']
F['n_exon_rec']    = st['n_exon_records']
F['exons_per_tx']  = st['n_exon_records'] / ntx['n_transcripts']
F['log_span']      = np.log10(st['span'] + 1)
F['log_exon_len']  = np.log10(st['total_exon_len'] + 1)

# ---- 2. expression per tissue x condition ------------------------------
ltpm = np.log2(tpm + 1)
grp  = meta.set_index('sample')
tissues = sorted(meta.tissue.unique())

for t in tissues:
    for cond in ['B', 'B0']:
        cols = grp[(grp.tissue == t) & (grp.condition == cond)].index
        F[f'ltpm_{t}_{cond}'] = ltpm[cols].mean(axis=1)
        F[f'cv_{t}_{cond}']   = tpm[cols].std(axis=1) / (tpm[cols].mean(axis=1) + 1e-6)

Bcols  = grp[grp.condition == 'B'].index
ltpmB  = ltpm[Bcols]
tis_means = pd.DataFrame({t: ltpm[grp[(grp.tissue==t)&(grp.condition=='B')].index].mean(axis=1)
                          for t in tissues})

# ---- 3. tissue-specificity (tau) ---------------------------------------
x   = tis_means.clip(lower=0)
mx  = x.max(axis=1)
tau = ((1 - x.div(mx.replace(0, np.nan), axis=0)).sum(axis=1)) / (x.shape[1] - 1)
F['tau']            = tau.fillna(0)
F['max_ltpm']       = mx
F['mean_ltpm']      = x.mean(axis=1)
F['n_tissue_expr']  = (tis_means > np.log2(2)).sum(axis=1)
F['top_tissue']     = tis_means.idxmax(axis=1)

# ---- 4. elimination-tissue contrasts (biology-driven) ------------------
# elimination-active: youngEm, MatEm, AR, LCM_Root ; non-elimination: Leaf, Root ; drive: Anther
elim   = ['youngEm', 'MatEm', 'AR', 'LCM_Root']
noelim = ['Leaf', 'Root']
F['ltpm_elim_mean']   = tis_means[elim].mean(axis=1)
F['ltpm_noelim_mean'] = tis_means[noelim].mean(axis=1)
F['elim_vs_noelim']   = F['ltpm_elim_mean'] - F['ltpm_noelim_mean']
F['ltpm_anther']      = tis_means['Anther']
F['anther_vs_leaf']   = tis_means['Anther'] - tis_means['Leaf']
F['root_vs_leaf']     = tis_means['Root']  - tis_means['Leaf']
F['AR_vs_Root']       = tis_means['AR']    - tis_means['Root']   # AR = B present, Root = B eliminated

# ---- 5. edgeR statistics per tissue ------------------------------------
for t in tissues:
    d = deg[deg.tissue == t].set_index('gene_id')
    F[f'logFC_{t}'] = d['logFC'].reindex(genes)
    F[f'FDR_{t}']   = d['FDR'].reindex(genes)
    F[f'sig_{t}']   = ((d['FDR'].reindex(genes) < 0.01) &
                       (d['logFC'].reindex(genes).abs() > 1)).astype(float)
    F[f'tested_{t}'] = d['logFC'].reindex(genes).notna().astype(int)

F['n_tissue_sig']  = F[[f'sig_{t}' for t in tissues]].sum(axis=1)
F['n_tissue_up']   = sum(((F[f'logFC_{t}'] > 1) & (F[f'FDR_{t}'] < 0.01)).astype(int) for t in tissues)
F['mean_logFC']    = F[[f'logFC_{t}' for t in tissues]].mean(axis=1)

F.to_csv(O+'feature_matrix_v1.csv.gz')
print('feature matrix:', F.shape)
print('numeric features:', F.select_dtypes(include=[np.number]).shape[1])
print('\nchrB vs A summary:')
print(F.groupby('is_B')[['tau','max_ltpm','elim_vs_noelim','AR_vs_Root','n_tissue_sig']].mean().round(3))
