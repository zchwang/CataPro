# Binding Site Mutation Analysis Guide

## 개요

UDH 클러스터링에서 사용한 binding site 잔기에 대해서만 mutation 분석을 수행하는 도구입니다. 각 core position 주변 ±1-3 window까지 포함하여 분석합니다.

## 주요 기능

### 1. Window 기반 잔기 분석
- **Core position (window 0)**: Binding site에 직접 관여하는 잔기
- **Window ±1-3**: Core position 주변 잔기
- 각 잔기가 어느 binding site의 어느 window에 속하는지 자동 라벨링

### 2. Single Mutation 분석
- 각 position에서의 mutation 효과 분석
- Wild-type vs Mutant 비교
- 통계적 유의성 검정 (Mann-Whitney U test)
- Effect size 계산

### 3. Pairwise Interaction 분석
- 두 position 간의 epistasis (상호작용) 분석
- Expected vs Observed double mutant 효과 비교
- 가능한 pair 수 계산 및 feasibility 평가

## 사용 방법

### Step 1: Binding Site Position 정의

```python
from binding_site_mutation_analyzer import BindingSiteMutationAnalyzer

# 방법 1: 리스트로 정의 (자동 라벨링)
binding_sites = [84, 85, 122, 155, 180]

# 방법 2: Dictionary로 정의 (custom 라벨)
binding_sites = {
    84: 'NAD_ribose',
    85: 'NAD_adenine',
    122: 'substrate_binding',
    155: 'catalytic_site',
    180: 'cofactor_binding'
}

# Window size 설정 (1, 2, or 3)
window_size = 3

# Analyzer 초기화
analyzer = BindingSiteMutationAnalyzer(binding_sites, window_size=window_size)
```

### Step 2: 분석 대상 확인

```python
# 분석 요약 정보
summary = analyzer.get_analysis_summary()

print(f"Core binding sites: {summary['n_core_sites']}")
print(f"Window size: ±{summary['window_size']}")
print(f"Total positions: {summary['total_positions']}")
print(f"Positions per site: {summary['positions_per_site']}")

# Pairwise 분석 가능 여부
from itertools import combinations
n_pairs = len(list(combinations(summary['all_positions'], 2)))
print(f"Possible pairwise combinations: {n_pairs}")
```

### Step 3: 데이터 로드 및 Mutation 식별

```python
# CataPro 예측 파일 로드
df = analyzer.load_and_filter_mutations(
    'catapro_predict_score_NADplus.csv',
    aligned_seq_col='aligned_seq'  # 정렬된 서열이 있는 컬럼명
)

# 자동으로 수행되는 작업:
# 1. Reference (wild-type) 서열 식별
# 2. 각 binding site position에서의 mutation 검출
# 3. Window 정보 부여
# 4. Kinetic parameter와 결합
```

### Step 4: Single Mutation 분석

```python
# 각 position의 mutation 효과 분석
single_results = analyzer.analyze_single_mutations(
    kinetic_params=[
        'pred_log10[kcat(s^-1)]',
        'pred_log10[Km(mM)]',
        'pred_log10[kcat/Km(s^-1mM^-1)]'
    ]
)

# 결과 확인
print(single_results.head())

# 유의한 결과만 필터링
significant = single_results[single_results['significant']]
print(f"Significant effects: {len(significant)}")
```

### Step 5: Pairwise Interaction 분석

```python
# Epistasis 분석
pairwise_results = analyzer.analyze_pairwise_interactions(
    kinetic_param='pred_log10[kcat/Km(s^-1mM^-1)]',
    min_samples=5  # 최소 샘플 수
)

# 유의한 상호작용
sig_interactions = pairwise_results[pairwise_results['significant']]
print(f"Significant interactions: {len(sig_interactions)}")

# Top epistatic interactions
top_epistasis = sig_interactions.nlargest(10, 'epistasis')
print(top_epistasis[['pos1', 'pos2', 'label1', 'label2', 'epistasis', 'p_value']])
```

### Step 6: 결과 저장

```python
# 리포트 생성
report_file = analyzer.generate_report(output_dir='./')

# 시각화 생성
analyzer.plot_results(output_dir='./')
```

## 출력 결과

### 1. DataFrame Columns

**mutations_df** (로드된 데이터):
- `pos_{N}`: Position N의 아미노산
- `pos_{N}_mutation`: Mutation 라벨 (예: "A85V" 또는 "A85A_WT")
- `pos_{N}_is_mutated`: Boolean (mutation 여부)
- `pos_{N}_site`: Binding site 이름
- `pos_{N}_window`: Window 번호 (-3 ~ +3, 0=core)
- `pos_{N}_label`: 전체 라벨 (예: "NAD_ribose_w+2")
- `n_mutations_in_sites`: 해당 서열의 총 mutation 수

**single_mutation_results**:
- `position`: 잔기 위치
- `site`: Binding site 이름
- `window`: Window 번호
- `label`: Position 라벨
- `core_position`: 해당하는 core position
- `parameter`: Kinetic parameter
- `n_wt` / `n_mut`: Wild-type / Mutant 샘플 수
- `wt_mean` / `mut_mean`: 평균값
- `delta`: Mutation 효과 (mut_mean - wt_mean)
- `effect_size`: Cohen's d effect size
- `p_value`: Statistical significance
- `significant`: p < 0.05 여부
- `most_common_mutation`: 가장 흔한 mutation
- `n_mutation_types`: Mutation 종류 수

**pairwise_results**:
- `pos1`, `pos2`: 두 위치
- `site1`, `site2`: 각 위치의 binding site
- `window1`, `window2`: 각 위치의 window
- `label1`, `label2`: Position 라벨
- `mean_wt_wt`, `mean_wt_mut`, `mean_mut_wt`, `mean_mut_mut`: 4가지 조합의 평균
- `single_effect_1`, `single_effect_2`: 각 single mutation 효과
- `expected_double`: 예상 double mutant 효과 (additive model)
- `observed_double`: 실제 관찰된 double mutant 효과
- `epistasis`: Epistatic effect (observed - expected)
- `p_value`: Interaction significance
- `significant`: p < 0.05 여부

### 2. 리포트 파일

`binding_site_mutation_report_{timestamp}.txt`:
- 분석 설정 정보
- 데이터 요약
- Single mutation 분석 결과
  - Top effects by parameter
- Pairwise interaction 분석 결과
  - Top positive epistasis
  - Top negative epistasis

### 3. 시각화 파일

**binding_site_single_mutations_{timestamp}.png**:
1. Effect size by position (bar plot)
2. Average effect size by window (bar plot)
3. Average effect size by binding site (horizontal bar)
4. P-value distribution (histogram)

**binding_site_pairwise_interactions_{timestamp}.png**:
1. Epistatic interaction matrix (heatmap)
2. Epistasis distribution (histogram)

## 예제 시나리오

### 시나리오 1: 5개 Core Sites, Window 3

```python
binding_sites = {
    84: 'NAD_ribose',
    85: 'NAD_adenine',
    122: 'substrate_binding',
    155: 'catalytic_site',
    180: 'cofactor_binding'
}

analyzer = BindingSiteMutationAnalyzer(binding_sites, window_size=3)
summary = analyzer.get_analysis_summary()
```

**결과:**
- Core binding sites: 5
- Total positions: 29 (5 × 7 positions per site, with overlaps removed)
- Positions per site: 7 (core + 6 window positions)
- Possible pairwise combinations: 406

**Feasibility:** ✓ 권장 (manageable size)

### 시나리오 2: 10개 Core Sites, Window 3

```python
binding_sites = {
    20: 'site1', 45: 'site2', 68: 'site3', 84: 'site4',
    85: 'site5', 122: 'site6', 155: 'site7', 180: 'site8',
    205: 'site9', 230: 'site10'
}

analyzer = BindingSiteMutationAnalyzer(binding_sites, window_size=3)
```

**결과:**
- Core binding sites: 10
- Total positions: ~70 (assuming no overlaps)
- Possible pairwise combinations: 2,415

**Feasibility:** ✓ 권장 (여전히 관리 가능)

### 시나리오 3: 20개 Core Sites, Window 3

**결과:**
- Total positions: ~140
- Possible pairwise combinations: 9,730

**Feasibility:** ⚠ 크지만 가능 (약간 느릴 수 있음)

### 시나리오 4: 50개 Core Sites, Window 3

**결과:**
- Total positions: ~350
- Possible pairwise combinations: 61,075

**Feasibility:** ⚠ 매우 큰 분석 (시간 소요, 메모리 고려 필요)

## Pairwise 분석 Feasibility Guide

| Core Sites | Total Positions | Pairs | 예상 시간 | 권장 |
|------------|----------------|-------|---------|------|
| 5 | ~29 | 406 | < 1분 | ✓ |
| 10 | ~70 | 2,415 | 1-5분 | ✓ |
| 15 | ~105 | 5,460 | 5-10분 | ✓ |
| 20 | ~140 | 9,730 | 10-20분 | △ |
| 30 | ~210 | 21,945 | 30-60분 | △ |
| 50 | ~350 | 61,075 | 1-2시간 | ⚠ |

**참고:**
- Window size를 줄이면 (3→2→1) position 수가 감소하여 더 빠름
- `min_samples` 값을 높이면 분석되는 pair가 줄어 더 빠름

## 데이터 요구사항

### 필수 컬럼

CataPro prediction 파일에 있어야 하는 컬럼:
1. **aligned_seq**: 정렬된 서열 (MSA 결과)
2. **pred_log10[kcat(s^-1)]**: kcat 예측값 (optional)
3. **pred_log10[Km(mM)]**: Km 예측값 (optional)
4. **pred_log10[kcat/Km(s^-1mM^-1)]**: kcat/Km 예측값 (optional)

### 권장 컬럼

- **type**: 'wild-type' 또는 'mutant' (reference 서열 식별용)
- **fasta_id**: Sequence identifier

## Tips

### 1. Position Indexing

- Python은 0-based indexing
- 하지만 결과는 1-based로 표시됨
- 예: Position 85 (1-based) = Index 84 (0-based)

### 2. Window 해석

- `window = 0`: Core binding site
- `window = +1`: Core position +1 (downstream)
- `window = -1`: Core position -1 (upstream)
- `window = +3`: Core position +3
- `window = -3`: Core position -3

### 3. Epistasis 해석

- **Positive epistasis (> 0)**: 두 mutation의 combined effect가 예상보다 큼
  - Synergistic effect
  - 두 mutation이 함께 있을 때 더 유리

- **Negative epistasis (< 0)**: Combined effect가 예상보다 작음
  - Antagonistic effect
  - 두 mutation이 서로 방해

- **Zero epistasis (≈ 0)**: Additive effect
  - 독립적인 효과

### 4. 성능 최적화

큰 데이터셋의 경우:
1. Window size 줄이기 (3 → 2 → 1)
2. `min_samples` 증가 (5 → 10 → 20)
3. Core sites 수 줄이기 (가장 중요한 것만)
4. Parameter 하나만 분석 (kcat/Km만)

## 문제 해결

### Q1: "Column 'aligned_seq' not found"
**A:** CSV 파일의 정렬된 서열 컬럼 이름을 확인하고 `aligned_seq_col` 파라미터에 정확한 이름 입력

### Q2: "Sequence length mismatch"
**A:** 모든 서열이 동일한 길이의 MSA 결과여야 함. Gap(-)으로 정렬 확인

### Q3: Pairwise 분석이 너무 느림
**A:**
- `min_samples` 값 증가
- Window size 감소
- Core sites 수 감소

### Q4: "Not enough positions for pairwise analysis"
**A:** 최소 2개 이상의 position이 필요. Binding sites 또는 window size 확인

## 참고

이 도구는 다음을 가정합니다:
1. 모든 서열이 동일한 MSA (Multiple Sequence Alignment)로 정렬됨
2. Position 번호는 정렬된 서열의 인덱스
3. Wild-type 또는 consensus 서열을 reference로 사용
4. Mutation은 단일 아미노산 치환

---

**생성 날짜:** 2025-11-04
**도구 버전:** 1.0
