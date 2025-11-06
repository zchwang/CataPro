# 보고된 UDH 분석

이 디렉토리는 **14개의 실험적으로 보고된 UDH 효소 전체**의 분석을 포함하며, 아미노산 분석, CataPro 예측, 클러스터링 분석을 종합적으로 비교합니다.

## ✨ 최신 업데이트 (2025-11-06)

**14개 효소 완전 분석 완료!**
- ✅ 모든 14개 효소 클러스터링 데이터 포함
- ✅ 새로 추가된 3개 효소: AtUdh, PpUdh, PsUdh
- ✅ 종합 분석 및 시각화 완료

## 주요 스크립트

### `analyze_reported_udhs.py`
메인 분석 스크립트:
- 14개 보고된 UDH 효소의 활성 부위 잔기 분석
- PDB 위치를 정렬 위치에 매핑
- 주요 위치에서 "최적" 아미노산 식별
- AA 점수 계산 및 돌연변이 제안
- 포괄적인 Excel 보고서 생성

### `prepare_reported_udh_for_catapro.py`
보고된 효소 서열에서 CataPro 입력 준비:
- FASTA 파일에서 서열 로드
- CataPro 예측을 위한 형식 지정
- 효소-기질 쌍이 포함된 입력 CSV 생성

### `visualize_reported_udh_rankings.py`
순위 비교 시각화 생성:
- AA 분석 순위 vs CataPro 순위
- 다양한 시각화 유형: 산점도, 막대 차트, 히트맵
- 전체 순위 비교

### `cluster_and_activity_analysis.py`
클러스터링 및 활성 종합 분석 (11개 효소):
- 8개 시나리오에서 클러스터링 결과 분석
- CataPro 예측 활성 프로파일
- 클러스터별 consensus sequence 계산
- 클러스터링 안정성 vs 활성 상관관계
- 기질 특이성 프로파일 비교
- 종합 시각화 (16개 서브플롯)

### `complete_14_enzyme_analysis.py` 🆕 **[최신]**
완전한 14개 효소 종합 분석:
- **14개 효소 전체** 클러스터링 및 활성 분석
- 새로 추가된 3개 효소 (AtUdh, PpUdh, PsUdh) 포함
- 기질별 활성 프로파일 (kcat, Km, kcat/Km)
- 상관관계 분석 (서열 길이, 클러스터링 안정성)
- Top 5 vs Bottom 5 비교
- 종합 시각화 (12개 서브플롯)

## 데이터 파일

### 입력 파일
- `reported_UDHs_complete.fasta` - 14개 보고된 효소의 서열
- `DWC0THCP211-alignment.aln` - 다중 서열 정렬

### 출력 파일
- `reported_udh_analysis.xlsx` - 포괄적인 Excel 보고서:
  - Overall_Summary: AA 점수별 효소 순위
  - Active_Site_Residues: 잔기별 상세 분석
  - Engineering_Suggestions: 돌연변이 권장사항
- `reported_udh_analysis_20251106_002815.xlsx` - 클러스터링 원본 데이터 (11개):
  - Summary_by_UDH: 효소별 클러스터링 통계
  - Reported_UDH_Distribution: 시나리오별 클러스터 할당
  - Full_S1~S8, Reported_S1~S8: 각 시나리오 상세 결과
- `reported_udh_analysis_20251106_010710.xlsx` 🆕 **[최신]** - 클러스터링 완전 데이터 (14개):
  - **14개 효소 모두 포함**
  - Summary_by_UDH, Reported_UDH_Distribution
  - Full_S1~S8, Reported_S1~S8
- `cluster_activity_consensus_analysis.xlsx` 🆕 - 종합 분석 결과:
  - UDH_Summary: 효소별 활성 및 클러스터링 통계
  - Consensus_Sequences: 클러스터별 consensus sequence
  - Clustering_Details: 전체 클러스터링 상세
  - Substrate_Specificity: 기질별 활성 매트릭스
  - Statistics: 상관관계 분석
- `reported_udh_catapro_input.csv` - CataPro 예측을 위한 형식 지정
- `reported_udh_catapro_predictions.csv` - CataPro가 예측한 동역학 매개변수
- `reported_udh_rankings_summary.csv` - 순위 비교
- `sequence_identity_matrix.csv` 🆕 - 서열 동일성 매트릭스
- `enzyme_cluster_stats.csv` 🆕 - 효소별 클러스터 통계
- `cluster_summary_stats.csv` 🆕 - 클러스터 요약 통계

### 시각화
- `reported_udh_rankings_scatter.png` - AA 점수 vs CataPro 산점도
- `reported_udh_rankings_bar_plots.png` - 막대 차트 비교
- `reported_udh_rankings_heatmap.png` - 잔기 패턴 히트맵
- `reported_udh_rankings_overall_ranking.png` - 순위 비교
- `cluster_activity_comprehensive_analysis.png` - 종합 분석 시각화 (11개):
  - 16개 서브플롯 (클러스터링 안정성, 활성, 기질 특이성, 시나리오별 분석)
- `complete_14_enzyme_analysis.png` 🆕 **[최신]** - 완전 분석 시각화 (14개):
  - 12개 서브플롯 (순위, 기질 특이성, 상관관계, Top 5 vs Bottom 5)

## 효소 목록 (14개 전체)

분석된 14개의 보고된 UDH 효소:
1. **AtUdh** - Agrobacterium fabrum (tumefaciens) - 7위, 완벽한 안정성
2. **CsUdh** - Chromohalobacter salixigens - 3위
3. **FpUdh** - Fulvimarina pelagi - 6위
4. **m1Udh** - Metagenomic 1 - 5위
5. **m2Udh** - Metagenomic 2 - **1위 (최고 활성)**
6. **m3Udh** - Metagenomic 3 - 4위
7. **OgUdh** - Oceanicola granulosus - 9위
8. **PfUdh** - Pseudomonas fluorescens - 11위, 완벽한 안정성
9. **PmUdh** - Pseudomonas mendocina - 13위
10. **PnUdh** - Polaromonas naphthalenviorans - 2위
11. **PpUdh** - Pseudomonas putida - 8위
12. **PsUdh** - Pseudomonas syringae - 10위, 완벽한 안정성
13. **SvUdh** - Streptomyces viridochromogenes - 14위, 완벽한 안정성
14. **TbUdh** - Thermobispora bispora - 12위

## 분석된 주요 위치

PDB 2F6Y (BsUdh) 기반 활성 부위 위치:
- K95, K96, K135 - 촉매 라이신
- W174, W175 - 기질 결합
- E140 - 양성자 전달
- Y197, Y198, Y221 - NAD+ 결합
- T130, T131, S132 - 기질 특이성
- L106, L107, L166 - 소수성 상호작용

## 주요 발견사항 (14개 전체)

### 촉매 활성 순위 (CataPro 예측)
평균 log10(kcat/Km) 기준 **Top 5**:
1. **m2Udh** (1.641) - Metagenomic 2, 모든 기질에 최고 활성
2. **PnUdh** (1.474) - Polaromonas naphthalenviorans
3. **CsUdh** (1.462) - Chromohalobacter salixigens
4. **m3Udh** (1.455) - Metagenomic 3
5. **m1Udh** (1.422) - Metagenomic 1

### 새로 추가된 3개 효소 🆕
- **AtUdh** - 7위 (1.339), **8/8 클러스터링 (완벽한 안정성)**
- **PpUdh** - 8위 (1.324), 5/8 클러스터링
- **PsUdh** - 10위 (1.318), **8/8 클러스터링 (완벽한 안정성)**

### 클러스터링 안정성 (14개 기준)
**완벽한 안정성 그룹 (Noise Rate 0%, 8/8 클러스터링):**
- **AtUdh** 🆕 (7위) - 안정성 그룹 중 **최고 활성**
- **PsUdh** 🆕 (10위)
- PfUdh (11위)
- SvUdh (14위)

**가장 불안정:**
- TbUdh (87.5% noise rate)
- FpUdh (75.0% noise rate)

### 클러스터링 vs 활성 (14개 기준)
- 상관계수: **r = 0.403 (p = 0.154)**
- **결론**: 여전히 유의미한 상관관계 없음
- 11개 → 14개로 확장해도 결과 동일

### 기질 특이성
m2Udh가 **모든 5가지 기질**에 대해 최고 활성 (변함없음):
- NAD+: 1.883 (2위 대비 +0.113)
- α-D-Glucuronate: 1.565 (2위 대비 +0.176)
- α-L-Guluronate: 1.607 (2위 대비 +0.188)
- β-D-Galacturonate: 1.618 (2위 대비 +0.180)
- β-D-Mannuronate: 1.533 (2위 대비 +0.181)

### Pseudomonas 그룹 완성 🆕
이제 4개 Pseudomonas 효소 비교 가능:
- PpUdh: 8위 (1.324) - 그룹 내 최고
- PsUdh: 10위 (1.318)
- PfUdh: 11위 (1.237)
- PmUdh: 13위 (1.095)
- **평균**: 1.244 (전체 평균 1.341보다 낮음)

## 사용법

전체 분석 파이프라인 실행:

```bash
# 1. 보고된 UDH 서열 분석
python analyze_reported_udhs.py

# 2. CataPro 예측 준비
python prepare_reported_udh_for_catapro.py

# 3. CataPro 예측 실행
python predict.py -inp_fpath reported_udh_catapro_input.csv \
                  -model_dpath ../models \
                  -batch_size 32 \
                  -device cpu \
                  -out_fpath reported_udh_catapro_predictions.csv

# 4. 순위 시각화
python visualize_reported_udh_rankings.py

# 5. 클러스터링 및 활성 종합 분석 (11개)
python cluster_and_activity_analysis.py

# 6. 완전한 14개 효소 분석 (최신) 🆕
python complete_14_enzyme_analysis.py
```

## 상세 분석 보고서

종합 분석 결과 및 고찰은 다음 문서를 참조:

### 14개 효소 완전 분석 🆕 **[최신 권장]**
- **`COMPLETE_14_ENZYME_REPORT_KO.md`** - 14개 효소 종합 보고서
  - 새로 추가된 3개 효소 상세 분석
  - 완전한 촉매 활성 순위 (14개)
  - Pseudomonas 그룹 완성 분석
  - 효소 공학 전략 (14개 기준)
  - 산업적 응용 가이드

### 11개 효소 초기 분석
- **`ANALYSIS_INSIGHTS_KO.md`** - 11개 효소 분석 보고서
  - 촉매 활성 순위 분석
  - 클러스터링 안정성 vs 활성 상관관계
  - Consensus sequence 분석
  - Binding site 시나리오별 인사이트

## 결론 (14개 전체 기준)

이 분석은 아미노산 분석과 CataPro 예측 간의 기준 비교를 확립하고, **14개 효소 전체** 클러스터링 분석을 통해 다음을 밝혀냈습니다:

1. **m2Udh의 압도적 우월성**: 모든 기질에 대해 최고 활성 (2위 대비 평균 0.168 log units 차이)
2. **새로운 안정성 그룹 발견**: AtUdh가 완벽한 안정성(8/8) + 높은 활성(7위) 조합
3. **Pseudomonas 그룹 완성**: 4개 효소 비교 가능, 안정성 높지만 활성은 중하위권
4. **클러스터링-활성 무관**: 14개로 확장해도 여전히 유의미한 상관관계 없음 (r=0.403, p=0.154)
5. **Metagenomic 효소의 우월성**: 평균 1.506 (Natural source 1.285보다 0.221 높음)

### 산업적 권장사항
- **최고 활성**: m2Udh
- **활성+안정성 균형**: PnUdh (추천 ⭐)
- **안정성 최우선**: AtUdh 또는 PsUdh
- **효소 공학 템플릿**: m2Udh (활성) + AtUdh (안정성) 하이브리드

상세 분석 결과는 **`COMPLETE_14_ENZYME_REPORT_KO.md`** 참조.
