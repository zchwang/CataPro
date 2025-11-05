# Sequence-to-Sequence Comparison Guide

## 🎯 목적

~1000개 서열을 서로 비교하여 **어떤 아미노산이 kinetic parameter를 증가/감소시키는지** 분석합니다.

### 이 분석이 알려주는 것:

1. **Position별 아미노산 효과**
   - Position 85에서 Val을 가진 서열들이 Ala보다 kcat이 높은가?
   - 각 position에서 "최고" vs "최저" 아미노산은?

2. **기질 특이성**
   - NAD에서는 Asp가 좋지만 Glucuronate에서는 Glu가 좋은가?
   - 어느 position이 기질 특이성을 결정하는가?

3. **CataPro 예측 패턴**
   - CataPro가 어떤 sequence를 높은 활성으로 예측하는가?
   - 어떤 residue combination이 최적인가?

## 📊 분석 방식

### Inter-molecular Comparison (서열 간 비교)

```
Sequence A: Position 85 = Val → kcat = 1.5
Sequence B: Position 85 = Ala → kcat = 1.2
Sequence C: Position 85 = Val → kcat = 1.6
Sequence D: Position 85 = Asp → kcat = 0.9

결론: Position 85에서 Val이 가장 높음 (평균 1.55)
      Asp가 가장 낮음 (평균 0.9)
      Delta = 0.65
```

이를 **각 binding site position**과 **NAD + 4 uronates**에 대해 수행

## 🔧 사용 방법

### 1. 필요한 파일 준비

각 기질에 대한 CataPro 예측 파일이 필요합니다:

```
inference/
├── catapro_predict_score_NADplus.csv
├── catapro_predict_score_α-D-Glucuronate.csv
├── catapro_predict_score_α-L-Guluronate.csv
├── catapro_predict_score_β-D-Galacturonate.csv
└── catapro_predict_score_β-D-Mannuronate.csv
```

**각 파일의 필수 컬럼:**
- `fasta_id`: 서열 ID
- `aligned_seq`: MSA 정렬된 서열
- `pred_log10[kcat(s^-1)]`: kcat 예측값
- `pred_log10[Km(mM)]`: Km 예측값
- `pred_log10[kcat/Km(s^-1mM^-1)]`: kcat/Km 예측값

### 2. Binding Site Positions 정의

```python
from sequence_comparison_analyzer import SequenceComparisonAnalyzer

# Clustering에서 사용한 binding site positions (0-indexed)
binding_sites = {
    84: 'NAD_ribose',
    85: 'NAD_adenine',
    122: 'substrate_binding',
    155: 'catalytic_site',
    180: 'cofactor_binding'
}

# Window size 설정
window_size = 3  # ±3 residues around each core position

# Analyzer 생성
analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=window_size)
```

### 3. 데이터 로드

```python
# 기질 파일 정의
substrate_files = {
    'NADplus': 'catapro_predict_score_NADplus.csv',
    'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
    'alpha_L_Guluronate': 'catapro_predict_score_α-L-Guluronate.csv',
    'beta_D_Galacturonate': 'catapro_predict_score_β-D-Galacturonate.csv',
    'beta_D_Mannuronate': 'catapro_predict_score_β-D-Mannuronate.csv'
}

# 데이터 로드
analyzer.load_substrate_data(substrate_files)

# Position 데이터 추출
analyzer.extract_position_data(aligned_seq_col='aligned_seq')
```

### 4. 아미노산 효과 분석

```python
# 각 position에서 아미노산 효과 분석
aa_effects = analyzer.analyze_amino_acid_effects(
    kinetic_params=['log_kcat', 'log_Km', 'log_kcat_Km'],
    min_sequences=5  # 최소 서열 수
)

# 결과 확인
print(aa_effects.head())

# 유의한 결과 필터링
significant = aa_effects[aa_effects['significant']]
print(f"Significant effects: {len(significant)}")
```

### 5. 기질 특이적 잔기 찾기

```python
# 기질에 따라 다른 아미노산을 선호하는 positions 찾기
substrate_specific = analyzer.find_substrate_specific_residues(
    parameter='log_kcat_Km',
    min_delta=0.1  # 최소 효과 크기
)

# 결과 확인
print(substrate_specific.head())
```

### 6. 결과 저장

```python
# 리포트 생성
analyzer.generate_report(output_dir='./')

# 시각화 생성
analyzer.plot_amino_acid_effects(output_dir='./')
```

### 7. 간편 실행

```bash
cd inference
python run_sequence_comparison.py
```

## 📈 결과 해석

### aa_effects DataFrame

각 행은 하나의 **position-substrate-parameter** 조합:

| 컬럼 | 설명 | 예시 |
|------|------|------|
| `position` | 잔기 위치 (1-based) | 85 |
| `site` | Binding site 이름 | NAD_adenine |
| `window` | Window 번호 | 0 (core) |
| `substrate` | 기질 이름 | NADplus |
| `parameter` | Kinetic parameter | log_kcat_Km |
| `best_aa` | 최고 아미노산 | V (Val) |
| `best_mean` | 최고 평균값 | 1.55 |
| `worst_aa` | 최저 아미노산 | D (Asp) |
| `worst_mean` | 최저 평균값 | 0.90 |
| `delta` | 효과 크기 | 0.65 |
| `p_value` | 통계적 유의성 | 0.001 |
| `significant` | p < 0.05 여부 | True |

**해석 예시:**
```
Position 85 (NAD_adenine core)
Substrate: NADplus
Parameter: log_kcat_Km
Best AA: V (mean=1.55)
Worst AA: D (mean=0.90)
Delta: 0.65
P-value: 0.001

→ Position 85에서 Val을 가진 서열들이 Asp보다
  log(kcat/Km)이 0.65 높음 (통계적으로 유의)
→ CataPro는 이 위치에 Val이 있는 서열을 더 활성 높다고 예측
```

### substrate_specific DataFrame

기질에 따라 선호하는 아미노산이 다른 positions:

| 컬럼 | 설명 |
|------|------|
| `position` | 잔기 위치 |
| `n_different_preferences` | 서로 다른 선호 아미노산 수 |
| `max_delta` | 최대 효과 크기 |
| `substrate_preferences` | 각 기질의 선호 아미노산 |

**해석 예시:**
```
Position 122 (substrate_binding core)
Different preferences: 3
Max delta: 0.8

NADplus prefers: A (Ala)
Glucuronate prefers: S (Ser)
Galacturonate prefers: T (Thr)

→ 이 위치는 기질 특이성 결정에 중요
→ Directed evolution 시 기질별로 다른 residue 선택 가능
```

## 📊 출력 파일

### 1. 리포트 (TXT)

`sequence_comparison_report_{timestamp}.txt`:

- 분석 설정 정보
- 데이터 요약
- 각 기질별 결과:
  - Top 10 positions by effect size
  - 유의한 효과 개수
- 기질 특이적 잔기 목록

### 2. 시각화 (PNG)

**aa_effects_{parameter}_{timestamp}.png** (parameter별):
1. Effect size heatmap (position × substrate)
2. Significant positions count by substrate

**substrate_comparison_top_positions_{timestamp}.png**:
- Top 10 positions의 기질별 비교
- 각 bar에 best amino acid 표시

## 💡 사용 예시

### 예시 1: NAD 특이적 잔기 찾기

```python
# NADplus에서 효과 큰 positions
nadplus_effects = aa_effects[
    (aa_effects['substrate'] == 'NADplus') &
    (aa_effects['parameter'] == 'log_kcat_Km') &
    (aa_effects['significant'])
].nlargest(10, 'delta')

for _, row in nadplus_effects.iterrows():
    print(f"Position {row['position']}: {row['best_aa']} > {row['worst_aa']} "
          f"(delta={row['delta']:.2f})")
```

### 예시 2: 특정 position의 기질 비교

```python
# Position 85를 모든 기질에 대해 비교
comparison = analyzer.compare_substrates(position=85, parameter='log_kcat_Km')
print(comparison)
```

출력:
```
     substrate best_aa  best_mean worst_aa  worst_mean  delta  p_value significant
      NADplus       V       1.55        D        0.90   0.65    0.001        True
  Glucuronate       S       1.42        D        0.95   0.47    0.010        True
 Galacturonate      T       1.38        A        0.92   0.46    0.015        True
```

### 예시 3: 최적 서열 설계 힌트

```python
# 각 position의 best amino acid 모으기
best_residues = {}

for substrate in ['NADplus', 'alpha_D_Glucuronate']:
    sub_effects = aa_effects[
        (aa_effects['substrate'] == substrate) &
        (aa_effects['parameter'] == 'log_kcat_Km') &
        (aa_effects['significant'])
    ]

    best_residues[substrate] = {}
    for pos in sub_effects['position'].unique():
        pos_data = sub_effects[sub_effects['position'] == pos]
        if len(pos_data) > 0:
            best_aa = pos_data.iloc[0]['best_aa']
            best_residues[substrate][pos] = best_aa

print("Optimal residues for NADplus:", best_residues['NADplus'])
print("Optimal residues for Glucuronate:", best_residues['alpha_D_Glucuronate'])
```

## 🔍 주의사항

### 1. 데이터 요구사항

- **모든 서열이 동일한 MSA로 정렬**되어 있어야 함
- Position 번호는 **aligned sequence의 index**
- 각 기질에 대해 **동일한 서열 세트** 사용 권장

### 2. 통계적 고려사항

- `min_sequences=5`: 너무 적으면 통계적 power 부족
- `min_sequences=20`: 더 신뢰할 수 있지만 많은 positions 제외될 수 있음
- P-value 보정 (Bonferroni 등) 고려 가능

### 3. 해석 주의

- **상관관계 ≠ 인과관계**
- 이 분석은 CataPro의 **예측 패턴**을 보여줌
- 실제 실험적 검증 필요
- 여러 position의 조합 효과 고려해야 함

## 🆚 다른 분석과의 비교

### vs. Pairwise Epistasis Analysis

| Feature | Sequence Comparison | Pairwise Epistasis |
|---------|---------------------|-------------------|
| 비교 대상 | 서로 다른 서열들 | 한 서열 내 position pairs |
| 질문 | "어떤 AA가 좋은가?" | "두 위치의 상호작용?" |
| 출력 | Position별 best AA | Position pair별 epistasis |
| 용도 | 서열 설계 힌트 | 상호작용 이해 |

### vs. Mutation Analysis

| Feature | Sequence Comparison | Mutation Analysis |
|---------|---------------------|-------------------|
| 기준 | 전체 서열 비교 | WT vs Mutant |
| 질문 | "절대적으로 어떤 AA가?" | "WT 대비 변화는?" |
| 데이터 | 모든 자연 변이 | 특정 mutations |

## 📊 예상 결과 규모

5 core sites, window 3, 5 substrates:

- **Total positions**: ~29
- **Parameters**: 3 (kcat, Km, kcat/Km)
- **Substrates**: 5
- **Total combinations**: ~435

각 combination당:
- Amino acids at position: 5-20 (varies)
- Sequences per AA: 5-200 (varies)

**Output:**
- ~400-500 rows in aa_effects
- ~10-50 substrate-specific positions

## 🚀 다음 단계

분석 후:

1. **상위 positions 검증**
   - 문헌에서 해당 위치의 역할 확인
   - 구조적 중요성 검토

2. **최적 서열 설계**
   - Best amino acids 조합
   - Directed evolution 라이브러리 설계

3. **실험 검증**
   - Site-directed mutagenesis
   - 실제 kinetic parameter 측정

4. **구조 분석**
   - AlphaFold 등으로 구조 예측
   - Best residue가 왜 좋은지 이해

---

**생성 날짜**: 2025-11-04
**도구 버전**: 1.0
