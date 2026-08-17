#!/bin/bash
#SBATCH --job-name=gp_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# --partition=<your_partition>

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geneprio

: "${PROJ_DIR:?set PROJ_DIR to the repository root}"
cd "$PROJ_DIR"

# keep BLAS single-threaded; sklearn parallelism is handled per-model
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python src/01_build_features.py
python src/02_parse_supplementary.py
python src/03_seq_features.py
python src/04_model_v2.py
python src/06_ablation.py
python src/07_final.py
python src/05_interpret.py
python src/08_figures.py

echo "done. results in out/ , figures in figs/"
