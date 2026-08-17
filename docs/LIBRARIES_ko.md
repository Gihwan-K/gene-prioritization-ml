# 어떤 라이브러리를 왜 썼나 — 그리고 선생님 이해가 맞는지 채점

## 채점부터

> **"RSEM까지는 무조건 내가 한 것처럼 직접 소프트웨어 돌려야 한다, 그치?"**

✅ **정확합니다.** `fastp → STAR → RSEM`은 하나도 안 바뀝니다. 실제 실행파일이고, 클러스터에서 돌려야 하고, 제가 한 일은 **선생님이 이미 만들어 놓은 RSEM 결과를 그대로 받아 쓴 것**입니다. 이 부분은 머신러닝과 아무 상관이 없습니다.

> **"근데 R에서 edgeR 대신에 다른 라이브러리 이용해서 결과 도출한 거다, 그치?"**

❌ **여기가 틀렸습니다. 두 군데가요.**

### 오해 ① — edgeR을 대체하지 않았습니다. **그대로 썼습니다.**

edgeR 결과는 이 프로젝트의 **입력**입니다. 111개 피처 중 **31개가 선생님의 edgeR 출력에서 직접 나왔습니다.**

| 출처 | 피처 수 |
|---|---:|
| RSEM TPM (제가 계산) | 39 |
| **edgeR 결과 (선생님이 만든 것)** | **31** |
| TransDecoder 단백질 서열 | 30 |
| GTF 어노테이션 | 7 |
| DIAMOND 상동성 | 3 |
| 공발현 | 1 |
| **합계** | **111** |

`logFC_AR`, `FDR_Anther`, `sig_LCM_Root`, `n_tissue_sig`… 전부 선생님이 R에서 돌린 `glmLRT` 결과입니다. **edgeR이 빠진 게 아니라, 최종 단계에서 입력 재료로 자리를 옮긴 겁니다.**

### 오해 ② — "대신"이 아니라 "그 다음"입니다

```
선생님 파이프라인
fastq → fastp → STAR → RSEM → edgeR → DEG 표 → [사람이 손으로 필터링] → 후보 8개
                                                 ↑
                                          여기서 끝났습니다

이 프로젝트
fastq → fastp → STAR → RSEM → edgeR → DEG 표 ─┐
                                              │
                                  GTF ────────┤
                      TransDecoder 단백질 ────┼──→ 피처 표 ──→ ML ──→ 3,176개 순위
                     GO (Supplementary) ──────┤    (111열)
                              DIAMOND ────────┘
                                                 ↑
                          [사람이 손으로 필터링] 을 여기로 교체
```

**edgeR까지는 완전히 동일합니다.** 바뀐 건 그 다음 — 선생님이 손으로 하시던 필터링 단계를 학습된 랭커로 바꾸고, 거기에 edgeR이 애초에 볼 수 없는 정보 4가지(유전자 구조, 단백질 서열, GO, 서열 상동성)를 추가로 넣은 겁니다.

edgeR은 발현량만 봅니다. 단백질에 라이신이 몇 %인지, 알려진 cohesin과 서열이 비슷한지는 edgeR이 알 방법이 없습니다. 그 정보를 넣는 게 이 프로젝트의 절반입니다.

### 오해 ③ (덤) — R이 아니라 Python입니다

R로도 **할 수 있습니다** (뒤에 §5에서 설명). 다만 저는 Python으로 했습니다.

---

## 정리하면 선생님 문장은 이렇게 고치면 맞습니다

> ~~"R에서 edgeR 대신 다른 라이브러리를 써서 결과를 도출했다"~~
>
> **"RSEM과 edgeR까지는 그대로 쓰고, 그 결과에 단백질 서열·유전자 구조·GO·상동성 정보를 추가로 붙여서, 손으로 하던 후보 선별을 Python 머신러닝 라이브러리로 대체했다"**

---

## 라이브러리 하나하나

### 요약표

| 라이브러리 | R로 치면 | 어디에 | 무엇을 |
|---|---|---|---|
| **pandas** | `dplyr` + `data.frame` | 전부 | 표 읽고 합치고 가공 |
| **numpy** | R 기본 벡터·행렬 | 전부 | 수치 계산 |
| **scipy.sparse.csgraph** | `igraph` | 06, 07 | 상동성 그래프에서 유전자 패밀리 찾기 |
| **scikit-learn** | `tidymodels` / `caret` | 04, 06, 07 | **모델 학습·교차검증·평가** |
| **shap** | `fastshap` / `shapr` | 05 | "왜 이 점수인가" 분해 |
| **matplotlib** | `ggplot2` | 08 | 그림 |
| **openpyxl** | `readxl` | 02 | Supplementary xlsx 읽기 |
| **DIAMOND** | `BLAST` | 03 앞단 | 단백질 서열 비교 (유일한 실행파일) |

### ① pandas — 표를 다루는 라이브러리

R의 `data.frame`과 `dplyr`을 합친 것입니다. **머신러닝과는 무관하고, 코드의 대부분이 이겁니다.**

```r
# R
counts <- read_csv("counts.csv.gz")
mean_tpm <- tpm %>% group_by(gene_id) %>% summarise(m = mean(value))
```
```python
# Python
counts = pd.read_csv("counts.csv.gz", index_col=0)
mean_tpm = tpm.mean(axis=1)
```

RSEM 매트릭스 읽기, GTF에서 만든 구조 표 붙이기, edgeR 결과를 조직별로 옆에 붙이기 — 전부 pandas입니다.

### ② numpy — 수치 계산

R은 벡터 연산이 언어에 내장되어 있지만, Python은 라이브러리로 합니다. `log2(tpm + 1)`, 상관계수 계산, 평균/표준편차 같은 것들입니다.

τ(조직 특이성) 계산이 numpy로 한 줄입니다:

```python
tau = ((1 - x.div(x.max(axis=1), axis=0)).sum(axis=1)) / (x.shape[1] - 1)
```

### ③ scipy.sparse.csgraph — 유전자 패밀리 찾기

**이게 왜 필요했냐면** — 파라로그 누출을 막으려면 "어떤 유전자들이 한 가족인가"를 먼저 정의해야 합니다.

DIAMOND가 "A와 B가 비슷하다, B와 C가 비슷하다"를 알려주면, 이걸 그래프로 보고 **연결 요소(connected component)** 를 찾으면 A-B-C가 한 가족이 됩니다. R의 `igraph::components()`와 같은 일입니다.

```python
ncc, labels = connected_components(adjacency_matrix)
# 3,176개 B 유전자 → 2,386개 패밀리
```

이 패밀리 번호를 `StratifiedGroupKFold`에 넘겨서 "한 가족은 통째로 같은 fold에" 넣었습니다.

### ④ scikit-learn — 여기가 진짜 머신러닝

R의 `tidymodels`나 `caret`에 해당합니다. 실제로 쓴 부품은 네 종류뿐입니다.

**모델**
```python
from sklearn.ensemble import HistGradientBoostingClassifier
```

**교차검증 분할기**
```python
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
```
`Stratified`는 각 fold에 양성 비율을 똑같이 유지한다는 뜻입니다. 양성이 3.4%뿐이라 무작위로 나누면 어떤 fold엔 양성이 거의 안 들어갈 수 있거든요.

**평가 지표**
```python
from sklearn.metrics import average_precision_score, roc_auc_score
```

**전처리 (v1 baseline에서만)**
```python
from sklearn.impute import SimpleImputer      # 결측값 채우기
from sklearn.preprocessing import StandardScaler  # 정규화
from sklearn.linear_model import LogisticRegression
```
최종 모델은 결측값과 스케일을 자체 처리해서 이건 안 씁니다.

### ⑤ shap — 해석

모델이 "이 유전자 점수 0.87"이라고만 하면 생물학적으로 쓸모가 없습니다. SHAP은 그 0.87을 **피처별 기여도로 분해**합니다.

```
CENH3-B의 점수 구성:
  aa_K            +0.804   ← 라이신 비율이 높다
  aa_T            +0.792
  tau             +0.707   ← 조직 특이적이다
  ltpm_LCM_Root_B +0.481   ← LCM 뿌리에서 발현이 높다
```

게임이론의 Shapley value를 적용한 방법입니다. R에도 `fastshap`, `shapr`이 있지만 Python 쪽이 더 성숙합니다.

### ⑥ DIAMOND — 유일하게 선생님한테 익숙한 도구

BLAST와 같은 일을 100~1000배 빠르게 합니다. 3만 개 단백질을 자기 자신과 전부 비교하는 데 **27초** 걸렸습니다. BLAST였으면 몇 시간이었을 겁니다.

```bash
diamond makedb --in prot_longest.faa -d self
diamond blastp -q prot_longest.faa -d self -o self_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 50 --evalue 1e-5 --threads 16
```

이건 conda로 설치하고 명령줄에서 돌립니다. 선생님이 체코 서버에서 하시던 것과 완전히 같은 방식입니다.

---

## 왜 이 모델을 골랐나 — 실제로 비교해봤습니다

세 가지를 다 돌려보고 골랐습니다 (v1 단계, 표준 CV 기준):

| 모델 | PR-AUC | 왜 이런 결과인가 |
|---|---|---|
| Logistic Regression | 0.078 | 선형이라 "발현이 높으면서 **동시에** 라이신이 많은 경우" 같은 조합을 못 잡음 |
| Random Forest | 0.091 | 트리를 독립적으로 많이 만들어 평균. 나쁘지 않음 |
| **HistGradientBoosting** | **0.104** | 트리를 **순차적으로** 쌓아 앞 트리의 실수를 뒤 트리가 보정 |

**HistGradientBoosting을 최종 선택한 이유:**

1. **성능이 제일 좋았습니다** (근소하지만)
2. **결측값을 자체 처리합니다.** 단백질 예측이 없는 유전자가 9,239개인데, 별도 처리 없이 그냥 넣으면 됩니다
3. **피처 스케일을 안 맞춰도 됩니다.** TPM은 0~수천, 아미노산 비율은 0~1인데 그냥 넣으면 됩니다
4. **scikit-learn 내장입니다.** XGBoost나 LightGBM도 성능은 비슷한데, 별도 설치가 필요합니다. 재현성 측면에서 의존성이 적은 게 낫습니다

**안 쓴 것:**

- **딥러닝(신경망)** — 2,371행 × 111열짜리 표 데이터에서는 트리 기반 모델이 거의 항상 이깁니다. 신경망은 이미지·텍스트·서열처럼 구조가 있는 대용량 데이터용입니다. 여기서 썼으면 성능도 나쁘고 "유행 따라간다"는 인상만 줬을 겁니다

### 왜 PR-AUC로 평가했나

양성이 3.4%뿐이라 **정확도는 무의미합니다** — 전부 "음성"이라 찍어도 96.6% 정확합니다. ROC-AUC도 불균형에서는 후하게 나옵니다.

PR-AUC는 "상위 N개를 뽑았을 때 그 안에 진짜가 얼마나 있나"를 종합한 값이라, 순위 매기기가 목적일 때 정직한 지표입니다.

---

## R로도 할 수 있나? — 네, 그리고 선생님한테는 그게 나을 수도 있습니다

Python이 필수는 아닙니다. R에도 다 있습니다.

| Python | R 대응 |
|---|---|
| scikit-learn | `tidymodels` (또는 `mlr3`, `caret`) |
| HistGradientBoosting | `xgboost`, `lightgbm`, `ranger` |
| StratifiedGroupKFold | `rsample::group_vfold_cv()` |
| average_precision_score | `yardstick::pr_auc()` |
| shap | `fastshap`, `shapr` |
| pandas | `dplyr` |

같은 분석을 R로 쓰면 이 정도입니다:

```r
library(tidymodels)

# df : gene_id, 111개 피처, label(factor), family(상동성 패밀리 번호)
folds <- group_vfold_cv(df, group = family, v = 5)

spec <- boost_tree(trees = 400, learn_rate = 0.06) |>
  set_engine("xgboost") |>
  set_mode("classification")

wf <- workflow() |>
  add_formula(label ~ .) |>
  add_model(spec)

res <- fit_resamples(wf, folds,
                     metrics = metric_set(pr_auc, roc_auc))
collect_metrics(res)
```

**진지하게 권합니다** — 시간이 되시면 이 프로젝트를 **R로 다시 짜보세요.** 이유가 세 가지입니다.

1. R은 이미 편하시니 코드가 눈에 들어옵니다. 개념 이해가 훨씬 빨라집니다
2. 결과가 제 Python 버전과 비슷하게 나오면, **본인이 이해했다는 증거**가 됩니다
3. 면접에서 "R로도 구현해봤습니다"는 강한 문장입니다

다만 지원까지 2주라면 지금 상태로 내시고, R 버전은 그 다음에 하셔도 됩니다.

---

## 마지막 정리 — 한 문장씩

- **fastp / STAR / RSEM** : 안 바뀜. 선생님이 하시던 그대로, 실행파일, 클러스터
- **edgeR** : 안 바뀜. **결과가 모델의 입력 재료(31개 피처)로 들어감**
- **바뀐 것** : edgeR 이후 "손으로 후보 좁히기" → "학습된 랭커로 3,176개 전부 순위 매기기"
- **추가된 것** : GTF 구조, 단백질 서열, GO, DIAMOND 상동성 — edgeR이 볼 수 없던 정보
- **쓴 도구** : Python 라이브러리 (실행파일 아님) + DIAMOND(실행파일 하나)
- **R로도 가능** : `tidymodels` + `xgboost`
