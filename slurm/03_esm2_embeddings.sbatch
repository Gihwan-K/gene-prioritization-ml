#!/bin/bash
#SBATCH --job-name=gp_esm2
#SBATCH --output=logs/esm2_%j.out
#SBATCH --error=logs/esm2_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
# --partition=<gpu_partition>
#
# OPTIONAL EXTENSION — could not be run in the offline build environment.
# Produces protein-language-model embeddings to replace amino-acid composition.
# Requires internet on the compute or login node to fetch model weights once.

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geneprio
: "${PROJ_DIR:?set PROJ_DIR}"
cd "$PROJ_DIR"

# one-time: pip install torch fair-esm
python - <<'PY'
import os, torch, esm, numpy as np
D = os.path.join(os.environ["PROJ_DIR"], "data")
O = os.path.join(os.environ["PROJ_DIR"], "out")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval().cuda()
bc = alphabet.get_batch_converter()

seqs, g = [], None
for line in open(os.path.join(D, "prot_longest.faa")):
    if line[0] == '>': g = line[1:].strip(); seqs.append([g, ""])
    else: seqs[-1][1] += line.strip()
seqs = [(g, s[:1022]) for g, s in seqs if s]

ids, embs, B = [], [], 8
with torch.no_grad():
    for i in range(0, len(seqs), B):
        lab, _, toks = bc(seqs[i:i+B])
        out = model(toks.cuda(), repr_layers=[33])["representations"][33]
        for j, (name, s) in enumerate(seqs[i:i+B]):
            ids.append(name)
            embs.append(out[j, 1:len(s)+1].mean(0).float().cpu().numpy())
        if i % 800 == 0: print(i, "/", len(seqs), flush=True)
np.savez_compressed(os.path.join(O, "esm2_embeddings.npz"),
                    ids=np.array(ids), emb=np.vstack(embs))
print("saved", len(ids), "embeddings")
PY
