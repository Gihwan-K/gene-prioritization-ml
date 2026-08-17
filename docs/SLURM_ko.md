# IPK 클러스터(SLURM)에서 돌리는 법

## 먼저 — 솔직한 이야기

**이 파이프라인은 클러스터가 필요 없습니다.** 전체가 노트북에서 20분이면 끝납니다. 가장 무거운 DIAMOND all-vs-all이 16코어로 27초입니다.

그래도 클러스터에서 돌리는 게 의미 있는 이유는 두 가지입니다.

1. **면접에서 물어봅니다.** "HPC 경험 있나요"에 "네, 이 프로젝트를 IPK SLURM에서 돌렸습니다"라고 답할 수 있는 것과 없는 것의 차이가 큽니다.
2. **확장할 때는 진짜로 필요합니다.** ESM-2 임베딩(GPU 필요)과 A 유전자 전체 확장은 노트북에서 못 합니다. 그리고 그 두 가지가 이 프로젝트의 가장 큰 개선 여지입니다.

그러니 **먼저 로그인 노드에서 그냥 돌려보시고**, 되는 걸 확인한 다음 sbatch로 옮기시길 권합니다.

---

## SLURM 기본 개념 (5분)

클러스터는 여러 대의 컴퓨터(노드) 묶음이고, SLURM은 그 순번을 관리하는 프로그램입니다.

| 개념 | 뜻 |
|---|---|
| **로그인 노드** | ssh로 접속했을 때 처음 닿는 곳. 여기서 무거운 계산을 돌리면 안 됩니다 (다른 사람들도 씁니다) |
| **계산 노드** | 실제 작업이 도는 곳. 직접 접속하지 않고 SLURM에 부탁합니다 |
| **작업(job)** | "이 명령을 이만큼의 자원으로 돌려줘"라는 요청서 |
| **파티션(partition)** | 계산 노드 그룹. 용도별로 나뉩니다 (일반/대용량메모리/GPU 등) |

명령어는 사실상 이 네 개면 됩니다.

```bash
sbatch  script.sbatch    # 작업 제출
squeue -u $USER          # 내 작업 상태 확인
scancel <jobid>          # 작업 취소
sinfo                    # 어떤 파티션이 있는지 확인
```

`.sbatch` 파일은 그냥 셸 스크립트인데, 맨 위에 `#SBATCH`로 시작하는 주석이 붙습니다. 그게 SLURM에게 보내는 요청서입니다.

```bash
#SBATCH --cpus-per-task=16    # CPU 16개 주세요
#SBATCH --mem=16G             # 메모리 16GB 주세요
#SBATCH --time=00:30:00       # 30분 안에 끝납니다
```

> ⚠️ `--time`을 실제보다 짧게 쓰면 **작업이 도중에 강제 종료됩니다.** 넉넉하게 잡으세요. 다만 너무 길게 잡으면 순번이 늦게 옵니다.

---

## 0단계 — 파티션 이름 확인

이 저장소의 sbatch 파일에는 `--partition`이 주석 처리되어 있습니다. IPK의 실제 파티션 이름을 모르기 때문입니다. 먼저 확인하세요.

```bash
sinfo -o "%P %a %l %D %c %m %G"
```

출력의 첫 열이 파티션 이름입니다 (`*`가 붙은 게 기본값). GPU 파티션은 마지막 `GRES` 열에 `gpu:...`가 표시됩니다.

확인했으면 각 sbatch 파일의 주석을 풀고 이름을 넣으세요.

```bash
#SBATCH --partition=fat        # 예시 — 실제 이름으로
```

---

## 1단계 — 환경 만들기 (한 번만)

로그인 노드에서 실행합니다.

```bash
# conda가 없다면 (있으면 건너뛰기)
module avail conda          # IPK에 모듈로 있는지 확인
module load anaconda3       # 혹은 miniconda3 — 이름은 환경마다 다릅니다

cd ~/gene-prioritization-ml
conda env create -f env/environment.yml
conda activate geneprio
```

`conda env create`는 5~10분 걸립니다. 끝나면 확인:

```bash
python -c "import sklearn, shap, pandas; print('ok')"
diamond --version
```

> conda가 아예 없고 설치도 번거로우면 `module load python/3.11` 후 `pip install --user -r requirements.txt`로도 됩니다. 다만 DIAMOND는 별도로 필요합니다 (`module avail diamond`로 확인).

---

## 2단계 — 데이터 배치

```
~/gene-prioritization-ml/
├── data/                     ← 여기에 입력 파일을 넣습니다
│   ├── counts.csv.gz             R 스크립트 출력
│   ├── tpm.csv.gz                R 스크립트 출력
│   ├── meta.csv                  R 스크립트 출력
│   ├── deg_full.csv.gz           R 스크립트 출력
│   ├── gene_struct.tsv           ← 아래 명령으로 생성
│   ├── gene_ntx.tsv              ← 아래 명령으로 생성
│   ├── prot.pep                  TransDecoder 결과 (압축 해제 상태)
│   └── suppl.xlsx                Supplementary Tables
├── src/  slurm/  env/  docs/
└── logs/                     ← 만들어두세요: mkdir logs
```

GTF에서 구조 테이블 두 개를 만드는 명령입니다. 파일러의 GTF를 직접 참조하면 됩니다.

```bash
cd ~/gene-prioritization-ml
mkdir -p data logs
GTF=/path/to/02_Aespelt3B_trans_v1.strand_filtered.gtf.gz

zcat $GTF | awk -F'\t' '$3=="exon"{
  match($9,/gene_id "[^"]*"/); g=substr($9,RSTART+9,RLENGTH-10);
  len=$5-$4+1; exlen[g]+=len; exn[g]++; chrom[g]=$1; strand[g]=$7;
  if(!(g in mn) || $4<mn[g]) mn[g]=$4;
  if(!(g in mx) || $5>mx[g]) mx[g]=$5 }
END{ print "gene_id\tchrom\tstrand\tgene_start\tgene_end\tspan\ttotal_exon_len\tn_exon_records";
  for(k in exlen) print k"\t"chrom[k]"\t"strand[k]"\t"mn[k]"\t"mx[k]"\t"mx[k]-mn[k]+1"\t"exlen[k]"\t"exn[k] }' \
  > data/gene_struct.tsv

zcat $GTF | awk -F'\t' '$3=="transcript"{
  match($9,/gene_id "[^"]*"/); g=substr($9,RSTART+9,RLENGTH-10); n[g]++ }
END{ print "gene_id\tn_transcripts"; for(k in n) print k"\t"n[k] }' \
  > data/gene_ntx.tsv

wc -l data/gene_struct.tsv data/gene_ntx.tsv   # 각각 40,147줄이어야 정상
```

---

## 3단계 — 제출

```bash
export PROJ_DIR=$HOME/gene-prioritization-ml
cd $PROJ_DIR

# ① 상동성 그래프 (~5분)
sbatch --export=ALL,PROJ_DIR=$PROJ_DIR slurm/01_diamond.sbatch
```

`Submitted batch job 1234567` 같은 메시지가 나옵니다. 상태 확인:

```bash
squeue -u $USER
```

`ST` 열이 `PD`면 대기 중, `R`이면 실행 중입니다. 사라지면 끝난 겁니다.

```bash
# 로그 확인 — 에러가 없는지 반드시 보세요
cat logs/diamond_1234567.out
cat logs/diamond_1234567.err
```

`data/self_hits.tsv`가 생겼고 30만 줄 이상이면 성공입니다.

```bash
# ② 본 파이프라인 (~30분)
sbatch --export=ALL,PROJ_DIR=$PROJ_DIR slurm/02_pipeline.sbatch
```

### 의존성 걸어서 한 번에 제출하기

①이 끝나야 ②가 의미 있으므로, 한 번에 넣고 싶으면:

```bash
JID=$(sbatch --parsable --export=ALL,PROJ_DIR=$PROJ_DIR slurm/01_diamond.sbatch)
sbatch --dependency=afterok:$JID --export=ALL,PROJ_DIR=$PROJ_DIR slurm/02_pipeline.sbatch
```

`afterok`는 "①이 **성공적으로** 끝나면 ② 시작"이라는 뜻입니다. ①이 실패하면 ②는 자동 취소됩니다.

---

## 4단계 — 결과 확인

```bash
ls -la out/ figs/
head -20 out/top50_candidates.csv
```

그림을 로컬 PC로 가져오려면 (본인 컴퓨터 터미널에서):

```bash
scp -r <사용자명>@<클러스터주소>:~/gene-prioritization-ml/figs ./
```

---

## 자주 만나는 문제

| 증상 | 원인과 해결 |
|---|---|
| `Invalid partition name specified` | `--partition` 이름이 틀림. `sinfo`로 확인 |
| 작업이 `PD`에서 안 움직임 | 순번 대기 중. `squeue -u $USER --start`로 예상 시작 시각 확인. 자원 요청을 줄이면 빨라집니다 |
| `CANCELLED ... DUE TO TIME LIMIT` | `--time`이 부족. 늘려서 재제출 |
| `slurmstepd: Exceeded job memory limit` | `--mem` 부족. 32G → 64G로 |
| `conda: command not found` | `module load anaconda3` 를 sbatch 안에 추가 |
| `ModuleNotFoundError` | `conda activate geneprio`가 sbatch 안에 있는지 확인. 로그인 노드의 환경은 계산 노드로 자동 전달되지 않습니다 |
| `PROJ_DIR: unbound variable` | `--export=ALL,PROJ_DIR=...` 를 빠뜨림 |

### 디버깅 요령 — 대화형 세션

에러가 나면 sbatch로 반복 제출하지 말고 계산 노드를 직접 잡아서 손으로 돌려보세요. 훨씬 빠릅니다.

```bash
srun --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash
conda activate geneprio
export PROJ_DIR=$HOME/gene-prioritization-ml
cd $PROJ_DIR
python src/01_build_features.py      # 한 단계씩 확인
exit
```

---

## 확장 — 여기서부터가 클러스터가 진짜 필요한 부분

### ① ESM-2 단백질 언어모델 (GPU)

`slurm/03_esm2_embeddings.sbatch`에 준비해뒀습니다. 제가 작업한 환경은 외부 네트워크가 막혀 모델 가중치를 받을 수 없었습니다. IPK에서는 될 겁니다.

```bash
conda activate geneprio
pip install torch fair-esm          # 로그인 노드에서 한 번
sbatch --export=ALL,PROJ_DIR=$PROJ_DIR slurm/03_esm2_embeddings.sbatch
```

650M 모델로 3만 개 단백질이 GPU에서 1~2시간입니다. 결과 `out/esm2_embeddings.npz`를 아미노산 조성 대신(또는 함께) 피처로 넣으면 됩니다. **성능이 가장 크게 오를 지점입니다.**

### ② InterProScan (Pfam 도메인)

IPK에 설치되어 있을 가능성이 높습니다 (`module avail interproscan`). 도메인 정보는 기능 예측에 강력한 피처입니다.

```bash
interproscan.sh -i data/prot_longest.faa -f TSV -o out/ipr.tsv \
  -appl Pfam -cpu 16 -goterms
```

`-goterms` 옵션을 주면 **GO 항목도 같이 나옵니다.** 즉 이걸로 A 유전자 GO 문제도 동료분 도움 없이 해결됩니다. 3만 개 단백질에 몇 시간 걸리니 `--time=24:00:00`으로 잡으세요.

### ③ A 유전자로 확장

②가 되면 라벨을 A 유전자까지 넓힐 수 있습니다. `src/_common.py`에서 `B = list(F.index[F.is_B==1])`를 전체 유전자로 바꾸고, "A로 학습 → B에 적용" 구조로 평가를 재구성하면 됩니다. 학습 세트가 2,371개에서 약 25,000개로 늘어납니다.

이 세 가지가 README의 "What this does not show"에 적어둔 한계를 그대로 해소합니다. 시간이 되면 순서대로 하시고, 안 되면 지금 상태로도 지원에 충분합니다.
