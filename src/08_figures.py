
import os
PROJ = os.environ.get("PROJ_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(PROJ, "data") + os.sep
O = os.path.join(PROJ, "out")  + os.sep
FG = os.path.join(PROJ, "figs") + os.sep
os.makedirs(O, exist_ok=True); os.makedirs(FG, exist_ok=True)
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SYN2,CENH3='Aesp_chrB_ipkv1.63465','Aesp_chrB_ipkv1.64685'

TH={'light':dict(surf='#fcfcfb',tp='#0b0b0b',ts='#52514e',grid='#e3e2de',
                 s1='#2a78d6',s2='#eb6834',s3='#1baf7a',mute='#c9c8c2'),
    'dark' :dict(surf='#1a1a19',tp='#ffffff',ts='#c3c2b7',grid='#3a3a37',
                 s1='#3987e5',s2='#d95926',s3='#199e70',mute='#55544f')}

def base(T,figsize):
    f,ax=plt.subplots(figsize=figsize,dpi=200)
    f.patch.set_facecolor(T['surf']); ax.set_facecolor(T['surf'])
    for s in ('top','right'): ax.spines[s].set_visible(False)
    for s in ('left','bottom'): ax.spines[s].set_color(T['grid']); ax.spines[s].set_linewidth(1)
    ax.tick_params(colors=T['ts'],labelsize=8,length=3,width=1)
    return f,ax
def title(f,ax,T,t,sub=None):
    """figure-level title block with reserved space; never collides with axes"""
    f.subplots_adjust(top=0.795 if sub else 0.865,left=0.155,right=0.975,bottom=0.165)
    f.text(0.012,0.975,t,color=T['tp'],fontsize=10.5,fontweight='600',va='top',ha='left',wrap=True)
    if sub: f.text(0.012,0.905,sub,color=T['ts'],fontsize=8.2,va='top',ha='left',wrap=True)

# ---------------- Fig 1 : biology QC ----------------
deg=pd.read_csv(D+'deg_full.csv.gz')
F=pd.read_csv(O+'feature_matrix_v1.csv.gz',index_col=0)
Bset=set(F.index[F.is_B==1])
cnts=deg.groupby('tissue')['gene_id'].apply(lambda s:len(set(s)&Bset)).sort_values(ascending=False)
for mode,T in TH.items():
    f,ax=base(T,(6.8,3.7))
    cols=[T['s2'] if t=='Root' else T['s1'] for t in cnts.index]
    b=ax.barh(range(len(cnts))[::-1],cnts.values,color=cols,height=.62)
    for i,(t,v) in enumerate(cnts.items()):
        ax.text(v+30,len(cnts)-1-i,f'{v:,}',va='center',color=T['tp'],fontsize=8.5,fontweight='600')
    ax.set_yticks(range(len(cnts))[::-1]); ax.set_yticklabels(cnts.index,color=T['ts'],fontsize=9)
    ax.set_xlim(0,2100); ax.set_xlabel('B-chromosome genes passing expression filter',color=T['ts'],fontsize=8.5)
    ax.xaxis.grid(True,color=T['grid'],lw=.8); ax.set_axisbelow(True)
    title(f,ax,T,'B-gene detection collapses in the primary root',
          'The one tissue where B chromosomes are already eliminated — recovered from the data alone')
    ax.text(0.98,0.06,'primary root',transform=ax.transAxes,ha='right',color=T['s2'],
            fontsize=8.5,fontweight='600')
    f.savefig(f'{FG}fig1_bgene_detection_{mode}.png',facecolor=T['surf']); plt.close(f)

# ---------------- Fig 2 : ablation ----------------
ab=pd.read_csv(O+'ablation.csv')
ab=ab[ab.config!='homology+coexpr only']
for mode,T in TH.items():
    f,ax=base(T,(7.4,4.1))
    x=np.arange(len(ab)); w=.36
    ax.bar(x-w/2,ab.pr_std,w,yerr=ab.pr_std_sd,color=T['s1'],label='standard 5-fold',
           error_kw=dict(ecolor=T['ts'],lw=1,capsize=2.5))
    ax.bar(x+w/2,ab.pr_grp,w,yerr=ab.pr_grp_sd,color=T['s2'],label='homology-grouped 5-fold',
           error_kw=dict(ecolor=T['ts'],lw=1,capsize=2.5))
    ax.axhline(0.034,color=T['ts'],ls=(0,(4,3)),lw=1.2)
    ax.text(-0.42,0.042,'random  0.034',color=T['ts'],fontsize=8,ha='left')
    for i,(a,sa,b_,sb) in enumerate(zip(ab.pr_std,ab.pr_std_sd,ab.pr_grp,ab.pr_grp_sd)):
        ax.text(i-w/2,a+sa+.016,f'{a:.2f}',ha='center',color=T['tp'],fontsize=8,fontweight='600')
        ax.text(i+w/2,b_+sb+.016,f'{b_:.2f}',ha='center',color=T['tp'],fontsize=8,fontweight='600')
    ax.set_xticks(x); ax.set_xticklabels(['expression\nonly','+ gene\nstructure','+ physico-\nchemical','+ aa\ncomposition','+ co-expr &\nhomology'],fontsize=8,color=T['ts'])
    ax.set_ylabel('PR-AUC',color=T['ts'],fontsize=8.5); ax.set_ylim(0,.60)
    ax.yaxis.grid(True,color=T['grid'],lw=.8); ax.set_axisbelow(True)
    lg=ax.legend(frameon=False,fontsize=8.5,loc='upper left',labelcolor=T['ts'])
    title(f,ax,T,'Protein sequence carries the signal',
          'Grouping folds by homology cluster removes paralog leakage — the gap between bars is that leakage')
    f.savefig(f'{FG}fig2_ablation_{mode}.png',facecolor=T['surf']); plt.close(f)

# ---------------- Fig 3 : two-axis map ----------------
res=pd.read_csv(O+'B_ranking_final.csv',index_col=0)
t10=set(pd.read_csv(O+'t10.csv').gene_id)&set(res.index)
for mode,T in TH.items():
    f,ax=base(T,(6.6,5.8))
    ax.scatter(res.func_disc_pct,res.elim_score,s=5,c=T['mute'],alpha=.45,lw=0)
    s=res.loc[list(t10)]
    ax.scatter(s.func_disc_pct,s.elim_score,s=34,c=T['s1'],lw=1.2,edgecolor=T['surf'],zorder=3)
    for nm,g,col in [('SYN2-B',SYN2,T['s2']),('CENH3-B',CENH3,T['s2'])]:
        p=res.loc[g]
        ax.scatter(p.func_disc_pct,p.elim_score,s=150,marker='*',c=col,lw=1.2,
                   edgecolor=T['surf'],zorder=5)
        ax.annotate(nm,(p.func_disc_pct,p.elim_score),textcoords='offset points',
                    xytext=(10,7),color=col,fontsize=9.5,fontweight='700',zorder=6)
    ax.set_xlabel('functional score  (percentile, discovery model)',color=T['ts'],fontsize=8.5)
    ax.set_ylabel('elimination-associated expression  (percentile)',color=T['ts'],fontsize=8.5)
    ax.grid(True,color=T['grid'],lw=.8); ax.set_axisbelow(True)
    ax.text(.985,.975,'top-right = high on both axes',transform=ax.transAxes,ha='right',va='top',
            color=T['ts'],fontsize=8,style='italic')
    h=[Line2D([],[],marker='*',ls='',ms=12,mfc=T['s2'],mec=T['surf'],label='validated candidates (held out)'),
       Line2D([],[],marker='o',ls='',ms=6,mfc=T['s1'],mec=T['surf'],label='conserved in sorghum (never trained on)'),
       Line2D([],[],marker='o',ls='',ms=5,mfc=T['mute'],mec=T['mute'],label='other B genes (n=3,176)')]
    ax.legend(handles=h,frameon=False,fontsize=8,loc='lower left',labelcolor=T['ts'])
    title(f,ax,T,'Neither axis alone recovers both candidates',
          'Functional score and expression pattern are complementary; their product ranks both in the top 4%')
    f.savefig(f'{FG}fig3_two_axis_{mode}.png',facecolor=T['surf']); plt.close(f)

# ---------------- Fig 4 : held-out recovery ----------------
schemes=[('functional\n(retrieval)','func_retrieval'),('functional\n(discovery)','func_discovery'),
         ('elimination\nexpression','elim_score'),('COMBINED','combined')]
rows=[]
for lab_,c in schemes:
    rk={g:i+1 for i,g in enumerate(res.sort_values(c,ascending=False).index)}
    rows.append((lab_,100*rk[SYN2]/len(res),100*rk[CENH3]/len(res)))
R=pd.DataFrame(rows,columns=['scheme','SYN2-B','CENH3-B'])
for mode,T in TH.items():
    f,ax=base(T,(7.0,4.1))
    x=np.arange(len(R)); w=.36
    ax.bar(x-w/2,R['SYN2-B'],w,color=T['s1'],label='SYN2-B')
    ax.bar(x+w/2,R['CENH3-B'],w,color=T['s2'],label='CENH3-B')
    for i in range(len(R)):
        ax.text(i-w/2,R['SYN2-B'][i]+.35,f"{R['SYN2-B'][i]:.1f}%",ha='center',color=T['tp'],fontsize=8,fontweight='600')
        ax.text(i+w/2,R['CENH3-B'][i]+.35,f"{R['CENH3-B'][i]:.1f}%",ha='center',color=T['tp'],fontsize=8,fontweight='600')
    ax.axhline(5,color=T['s3'],ls=(0,(4,3)),lw=1.4)
    ax.text(3.45,5.4,'top 5%',color=T['s3'],fontsize=8.5,ha='right',fontweight='600')
    ax.set_xticks(x); ax.set_xticklabels(R.scheme,fontsize=8.5,color=T['ts'])
    ax.set_ylabel('rank percentile  (lower = better)',color=T['ts'],fontsize=8.5)
    ax.set_ylim(0,20); ax.yaxis.grid(True,color=T['grid'],lw=.8); ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=8.5,labelcolor=T['ts'],loc='upper right')
    title(f,ax,T,'Both validated candidates land in the top 4% of 3,176 B genes',
          'Neither gene was seen during training at any stage')
    f.savefig(f'{FG}fig4_recovery_{mode}.png',facecolor=T['surf']); plt.close(f)
print('figures written')
