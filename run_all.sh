#!/usr/bin/env bash
# Full pipeline (local / login node). See slurm/ for cluster submission.
# Expects input files in $PROJ_DIR/data — see docs/SLURM_ko.md § 2단계.
set -euo pipefail
export PROJ_DIR="${PROJ_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$PROJ_DIR"

# longest ORF per gene, then all-vs-all homology
python - <<'PY'
import os, re
D = os.path.join(os.environ["PROJ_DIR"], "data")
best, name, seq = {}, None, []
def flush():
    if name is None: return
    g = re.sub(r'\.\d+\.p\d+$','',name.split()[0]); s=''.join(seq)
    if g not in best or len(s) > len(best[g]): best[g] = s
for line in open(os.path.join(D,"prot.pep")):
    if line[0]=='>': flush(); name=line[1:].strip(); seq=[]
    else: seq.append(line.strip())
flush()
with open(os.path.join(D,"prot_longest.faa"),"w") as o:
    for g,s in best.items(): o.write(f">{g}\n{s.replace('*','')}\n")
print("proteins:", len(best))
PY
diamond makedb --in data/prot_longest.faa -d data/self --quiet
diamond blastp -q data/prot_longest.faa -d data/self -o data/self_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 50 --evalue 1e-5 --threads 8 --quiet

python src/01_build_features.py
python src/02_parse_supplementary.py
python src/03_seq_features.py
python src/04_model_v2.py
python src/06_ablation.py
python src/07_final.py
python src/05_interpret.py
python src/08_figures.py
echo "done — results in out/, figures in figs/"
