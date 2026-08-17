#!/bin/bash
#SBATCH --job-name=gp_diamond
#SBATCH --output=logs/diamond_%j.out
#SBATCH --error=logs/diamond_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
# --partition=<your_partition>          # <-- ask IPK IT / check `sinfo`

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geneprio

: "${PROJ_DIR:?set PROJ_DIR to the repository root}"
cd "$PROJ_DIR"

# 1) longest ORF per gene from the TransDecoder peptides
python - <<'PY'
import os, re
D = os.path.join(os.environ["PROJ_DIR"], "data")
best, name, seq = {}, None, []
def flush():
    if name is None: return
    g = re.sub(r'\.\d+\.p\d+$', '', name.split()[0]); s = ''.join(seq)
    if g not in best or len(s) > len(best[g]): best[g] = s
src = os.path.join(D, "prot.pep")
for line in open(src):
    if line[0] == '>': flush(); name = line[1:].strip(); seq = []
    else: seq.append(line.strip())
flush()
with open(os.path.join(D, "prot_longest.faa"), "w") as o:
    for g, s in best.items(): o.write(f">{g}\n{s.replace('*','')}\n")
print("proteins:", len(best))
PY

# 2) all-vs-all homology graph
diamond makedb --in data/prot_longest.faa -d data/self --threads "$SLURM_CPUS_PER_TASK"
diamond blastp \
  -q data/prot_longest.faa -d data/self -o data/self_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 50 --evalue 1e-5 \
  --threads "$SLURM_CPUS_PER_TASK"

wc -l data/self_hits.tsv
