# Position Selection Summary - UDH Sequence Analysis

## 초기 정의 (Initial Definition)

### Core Binding Sites (8개)
PDB 3RFV (Uronate dehydrogenase) 기준:

1. **Position 12**: NAD_binding_1
2. **Position 13**: NAD_binding_2
3. **Position 14**: NAD_binding_3
4. **Position 87**: catalytic_1
5. **Position 110**: catalytic_2
6. **Position 136**: substrate_1
7. **Position 159**: substrate_2
8. **Position 181**: cofactor

### Window 설정
- 각 core site 주변 **±3 residues** (total 6 windows per site)
- Window labels: w-3, w-2, w-1, core, w+1, w+2, w+3
- **이론적 최대 positions**: 8 core × 7 (core + 6 windows) = **56 positions**

---

## 1차 선별: Alignment Coverage 필터링

### 제거 기준
- **Gap coverage > 80%** (min_coverage = 0.2)
- Alignment에서 해당 position이 대부분 gap인 경우

### 제거된 Positions (20개)

#### NAD_binding 관련 (8개)
| Position | Label | Coverage | Gaps | Reason |
|----------|-------|----------|------|--------|
| 9 | NAD_binding_1_w-3 | 0.3% | 5890/5905 | Insertion/deletion region |
| 10 | NAD_binding_2_w-3 | 4.0% | 5670/5905 | Alignment artifact |
| 11 | NAD_binding_3_w-3 | 3.0% | 5730/5905 | Alignment artifact |
| 12 | NAD_binding_3_w-2 | 0.6% | 5870/5905 | Insertion/deletion region |
| 13 | NAD_binding_3_w-1 | 0.6% | 5870/5905 | Insertion/deletion region |
| 14 | NAD_binding_3_core | 0.6% | 5870/5905 | **Core site이지만 gap 과다** |
| 16 | NAD_binding_3_w+2 | 16.9% | 4910/5905 | Low coverage |
| 17 | NAD_binding_3_w+3 | 10.8% | 5265/5905 | Low coverage |

**➡️ NAD_binding sites에서 실제 분석된 position: Position 15 (NAD_binding_3_w+1) 단 1개만 유지**

#### Catalytic_1 관련 (3개)
| Position | Label | Coverage | Gaps | Reason |
|----------|-------|----------|------|--------|
| 84 | catalytic_1_w-3 | 2.1% | 5780/5905 | Insertion/deletion region |
| 85 | catalytic_1_w-2 | 9.1% | 5365/5905 | Low coverage |
| 87 | catalytic_1_core | 2.5% | 5760/5905 | **Core site이지만 gap 과다** |

**➡️ Catalytic_1에서 실제 분석된 positions: 86, 88, 89, 90 (4개)**

#### Catalytic_2 관련 (5개)
| Position | Label | Coverage | Gaps | Reason |
|----------|-------|----------|------|--------|
| 110 | catalytic_2_core | 0.1% | 5900/5905 | **Core site이지만 gap 과다** |
| 111 | catalytic_2_w+1 | 14.9% | 5025/5905 | Low coverage |
| 112 | catalytic_2_w+2 | 19.3% | 4765/5905 | Low coverage |
| 113 | catalytic_2_w+3 | 19.4% | 4760/5905 | Low coverage |

**➡️ Catalytic_2에서 실제 분석된 positions: 107, 108, 109 (3개, 모두 upstream windows)**

#### Substrate_1 관련 (2개)
| Position | Label | Coverage | Gaps | Reason |
|----------|-------|----------|------|--------|
| 133 | substrate_1_w-3 | 0.5% | 5875/5905 | Insertion/deletion region |
| 138 | substrate_1_w+2 | 0.1% | 5900/5905 | Insertion/deletion region |

**➡️ Substrate_1에서 실제 분석된 positions: 134, 135, 136, 137, 139 (5개)**

#### Substrate_2 관련 (3개)
| Position | Label | Coverage | Gaps | Reason |
|----------|-------|----------|------|--------|
| 158 | substrate_2_w-1 | 0.4% | 5880/5905 | Insertion/deletion region |
| 160 | substrate_2_w+1 | 0.1% | 5900/5905 | Insertion/deletion region |
| 161 | substrate_2_w+2 | 0.7% | 5865/5905 | Insertion/deletion region |

**➡️ Substrate_2에서 실제 분석된 positions: 156, 157, 159, 162 (4개)**

#### Cofactor 관련 (0개 제거)
**➡️ Cofactor에서 실제 분석된 positions: 178, 179, 180, 181, 182, 183, 184 (7개 모두 유지)**

---

## 2차 선별: Sample Size 필터링

### 제거 기준
- 각 amino acid의 **n < 5** (최소 5개 sequences 필요)
- n = 1인 amino acids는 자동 제외 (unreliable)

### 결과
- Position별로 분석 가능한 amino acid 개수가 다름
- 예: Position 15 (NAD_binding_3_w+1)
  - Total unique AAs: 8개
  - n=1 (excluded): 7개
  - **n≥5 (usable): 단 1개만** ← 거의 분석 불가능

---

## 최종 분석된 Positions (24개)

### Site별 분포

| Site | Core Position | Windows Analyzed | Total Positions | Coverage |
|------|--------------|------------------|-----------------|----------|
| NAD_binding | 12, 13, 14 | Only w+1 (pos 15) | 1 | 12.5% |
| Catalytic_1 | 87 | w-1, w+1, w+2, w+3 (86, 88, 89, 90) | 4 | 57% |
| Catalytic_2 | 110 | w-3, w-2, w-1 (107, 108, 109) | 3 | 43% |
| Substrate_1 | 136 | w-2, w-1, core, w+1, w+3 (134-137, 139) | 5 | 71% |
| Substrate_2 | 159 | w-3, w-2, core, w+3 (156, 157, 159, 162) | 4 | 57% |
| Cofactor | 181 | All 7 positions (178-184) | 7 | **100%** |

### 주요 발견

#### ✅ 잘 보존된 regions
1. **Cofactor binding site (181)**: 모든 window 유지 (7/7)
2. **Substrate_1 (136)**: 5/7 유지
3. **Catalytic_1 (87)**: 4/7 유지 (core 제외)

#### ❌ 많이 손실된 regions
1. **NAD_binding sites**: 24개 중 1개만 유지 (4%)
   - Core sites (12, 13, 14) 모두 gap 과다로 제거
2. **Catalytic_2 core (110)**: Gap 과다로 제거
3. **Catalytic_1 core (87)**: Gap 과다로 제거

---

## 손실 원인 분석

### 1. Alignment Quality Issues
- **N-terminal region (positions 9-17)**:
  - 많은 sequences에서 N-terminal이 짧거나 missing
  - NAD binding domain의 구조적 다양성

### 2. Insertion/Deletion Hotspots
- Position 87, 110: Core catalytic sites이지만 gap 과다
  - → 다른 UDH variants에서 해당 위치에 insertion/deletion 빈번
  - → Sequence alignment에서 gap으로 처리됨

### 3. Evolutionary Variability
- NAD binding region이 특히 variable
  - 다른 dehydrogenase family members 포함 시 alignment 어려움

---

## 통계 요약

| Category | Count | Percentage |
|----------|-------|------------|
| **Initial defined positions** | ~44 | 100% |
| **Removed by coverage filter** | 20 | 45% |
| **Final analyzed positions** | 24 | 55% |
| | | |
| **Core sites defined** | 8 | 100% |
| **Core sites analyzed** | 3 | 37.5% |
| **Core sites lost** | 5 | 62.5% |

---

## Significant Findings (p<0.05)

### By Site Type
| Site Type | Positions Analyzed | Significant Positions | Hit Rate |
|-----------|-------------------|----------------------|----------|
| Substrate_1 | 5 | 3-4 | 60-80% |
| Substrate_2 | 4 | 3-4 | 75-100% |
| Catalytic_1 | 4 | 3-4 | 75-100% |
| Catalytic_2 | 3 | 1-2 | 33-67% |
| Cofactor | 7 | 4-6 | 57-86% |
| NAD_binding | 1 | 0 | 0% |

**➡️ Substrate binding sites와 Catalytic sites 주변이 가장 중요한 determinants**

---

## Recommendations

### 1. NAD Binding Region 분석 개선
- 현재 거의 분석 불가 (1/24 positions)
- 해결 방법:
  - Sequence curation: N-terminal이 complete한 sequences만 선별
  - Manual alignment correction
  - Alternative alignment tool 시도

### 2. Core Site Gap 문제 해결
- Position 87, 110이 core site임에도 분석 불가
- 해결 방법:
  - Gap 주변 window positions 분석으로 간접 추론
  - Structure-based analysis로 보완

### 3. 현재 분석 가능한 영역에 집중
- **Substrate binding sites (136, 159)**: 가장 완전한 데이터
- **Catalytic windows (86, 88-90, 107-109)**: Core는 빠졌지만 주변 잘 보존
- **Cofactor site (178-184)**: 완벽한 coverage

---

Generated: 2025-11-05
