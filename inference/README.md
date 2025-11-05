# CataPro Inference - UDH 효소 분석

이 디렉토리는 전통적인 아미노산 분석과 머신러닝 예측을 비교하는 UDH (Uronate Dehydrogenase) 효소의 종합 분석을 포함합니다.

## 프로젝트 구조

```
inference/
├── mutation_analysis/       # 돌연변이 검증 연구
├── correlation_analysis/    # AA 분석 vs CataPro 상관관계
├── reported_udh_analysis/   # 14개 보고된 효소 분석
├── documentation/          # 가이드 및 문서
└── [core files]           # CataPro 예측 스크립트
```

## 분석 개요

### 1. 보고된 UDH 분석 (`reported_udh_analysis/`)
**목적**: 14개의 실험적으로 특성화된 UDH 효소를 아미노산 빈도 분석을 사용하여 분석

**방법**:
- 1181개의 UDH 서열에 대한 다중 서열 정렬
- 15개 활성 부위 위치에서 "최적" 아미노산의 통계적 식별
- 최적 잔기 일치도를 기반으로 한 점수 매기기 및 순위 결정
- 효소 성능 향상을 위한 돌연변이 제안

**주요 출력**: `reported_udh_analysis.xlsx` - 효소 순위 및 돌연변이 권장사항

### 2. 상관관계 분석 (`correlation_analysis/`)
**목적**: AA 분석 순위와 CataPro ML 예측 비교

**방법**:
- NAD+ 보조인자를 사용한 14개 효소 모두에 대한 CataPro 예측
- Pearson 및 Spearman 상관관계 분석
- 순위 비교 시각화

**주요 발견**: **음의 상관관계** (r = -0.573, p=0.032)
- AA 점수가 높은 효소는 예측 성능이 낮음
- 최적 AA 일치가 많은 효소는 성능이 더 나쁨
- 더 많은 돌연변이가 필요한 효소가 성능이 더 좋음

**함의**: 통계적 아미노산 분석은 역효과를 내는 엔지니어링 제안을 생성할 수 있음

### 3. 돌연변이 검증 (`mutation_analysis/`)
**목적**: 제안된 돌연변이가 실제로 성능을 향상시키는지 실증적으로 테스트

**방법**:
- 280개의 돌연변이 서열 생성 (14개 효소 × 3개 돌연변이 세트 × 5개 기질 + 야생형)
- top5, top10 및 모든 제안된 돌연변이 적용
- 모든 돌연변이에 대해 CataPro를 사용하여 성능 예측
- 돌연변이 대 야생형 성능 비교

**주요 결과**:
- **Top5 돌연변이**: 71.4% 저하 (Δkcat/Km = -0.0476)
- **Top10 돌연변이**: 71.4% 저하 (Δkcat/Km = -0.1032)
- **모든 돌연변이**: 94.3% 저하 (Δkcat/Km = -0.6143)
- **전체**: 79%의 돌연변이가 야생형보다 성능이 나쁨

**결론**: **가설 확인** - 제안된 돌연변이가 실제로 효소 성능을 저하시킴

## 핵심 발견사항

### 주요 발견
전통적인 아미노산 빈도 분석은 효소 엔지니어링에 **역효과**를 내는 돌연변이 제안을 생성합니다:

1. **통계적 빈도 ≠ 최적 성능**: 위치에서 가장 흔한 아미노산이 반드시 더 나은 촉매 효율을 가져오지는 않음

2. **진화적 맥락이 중요함**: 고성능 효소(예: m2Udh)는 합의에서 벗어나지만 특정 기능에 최적화된 고유한 잔기를 가질 수 있음

3. **머신러닝이 복잡성을 포착함**: CataPro의 딥러닝 접근법은 단순한 통계 분석이 놓치는 복잡한 서열-기능 관계를 포착함

4. **엔지니어링 함의**: 효소 엔지니어링 전략은 순진한 빈도 기반 접근법을 피해야 함

### 효소별 사례

**m2Udh** (CataPro #1, AA 분석 #14):
- AA 분석에서 꼴찌 순위 (최저 점수)
- CataPro에서 1위 순위 (최고 예측 성능)
- 모든 제안된 돌연변이가 성능 저하 (-0.01 ~ -0.24)
- **해석**: 이미 최적화됨, 정당한 이유로 "합의"에서 벗어남

**AtUdh** (AA 분석 #1, CataPro #7):
- AA 분석에서 1위 순위 (최고 점수)
- CataPro에서 7위 순위 (중간 성능)
- 모든 제안된 돌연변이가 성능 저하 (-0.11 ~ -0.29)
- **해석**: 높은 합의 일치가 실제 성능을 예측하지 못함

## 핵심 CataPro 파일

### 예측 스크립트
- `predict.py` - 메인 CataPro 예측 스크립트
- `model.py` - 모델 아키텍처 (ProtT5 + MolT5)
- `act_model.py` - 활성화 모델 구성요소
- `utils.py` - 유틸리티 함수

### 분석 지원 스크립트
- `parse_alignment_simple.py` - CLUSTAL 정렬 파싱
- `sequence_comparison_analyzer.py` - 서열 비교 유틸리티
- `binding_site_mutation_analyzer.py` - 돌연변이 분석 도구
- `prepare_analysis_data.py` - 데이터 준비 유틸리티
- `run_sequence_comparison.py` - 서열 비교 실행
- `run_analysis_with_existing_data.py` - 분석 파이프라인

### 샘플 데이터
- `sample_inp.csv` - CataPro 입력 형식 예제
- `catapro_predict_score.csv` - 대규모 예측 결과
- `catapro_test-pred.csv` - 테스트 예측

## 사용법

### 전체 분석 파이프라인 실행

```bash
# 1. 보고된 UDH 서열 분석
cd reported_udh_analysis
python analyze_reported_udhs.py

# 2. CataPro 예측 실행
python prepare_reported_udh_for_catapro.py
python ../predict.py -inp_fpath reported_udh_catapro_input.csv \
                     -model_dpath ../../models \
                     -batch_size 32 \
                     -device cpu \
                     -out_fpath reported_udh_catapro_predictions.csv

# 3. 순위 시각화
python visualize_reported_udh_rankings.py

# 4. 상관관계 분석
cd ../correlation_analysis
python compare_mutation_impact_en.py
python plot_correlation_english.py

# 5. 돌연변이 생성 및 테스트
cd ../mutation_analysis
python create_mutant_sequences.py
python ../predict.py -inp_fpath mutant_sequences_catapro_input.csv \
                     -model_dpath ../../models \
                     -batch_size 32 \
                     -device cpu \
                     -out_fpath mutant_sequences_predictions.csv

# 6. 돌연변이 영향 분석
python analyze_mutant_impact.py
```

### 단일 CataPro 예측 실행

```bash
python predict.py \
    -inp_fpath input.csv \
    -model_dpath ../models \
    -batch_size 32 \
    -device cpu \
    -out_fpath output.csv
```

입력 CSV 형식:
```csv
,Enzyme_id,type,sequence,smiles
0,enzyme1,wild-type,MKLV...,C1=CC=C...
```

## 주요 시사점

1. **연구자를 위한 조언**:
   - 효소 엔지니어링을 위한 통계적 아미노산 분석에 주의
   - 계산 예측을 실증적으로 검증
   - 진화적 맥락 및 구조적 정보 고려

2. **엔지니어를 위한 조언**:
   - 머신러닝 접근법(CataPro 등)은 단순 통계보다 더 많은 복잡성을 포착
   - 직관에 반하는 잔기가 기능적으로 중요할 수 있음
   - 엔지니어링 전략은 여러 정보 소스를 통합해야 함

3. **방법론 개발자를 위한 조언**:
   - 빈도 기반 분석에는 근본적인 한계가 있음
   - 음성 검증 연구는 분야에 가치가 있음
   - ML 모델은 예측 패턴을 통해 기계론적 통찰을 제공할 수 있음

## 참고문헌

- **CataPro**: Wang et al., 단백질 및 분자 언어 모델을 사용한 효소 동역학 매개변수 예측
- **PDB 구조**: 2F6Y (BsUdh)를 활성 부위 위치 참조로 사용
- **데이터셋**: 1181개 UDH 효소 서열의 다중 서열 정렬
- **검증 세트**: 14개의 실험적으로 특성화된 UDH 효소

## 인용

이 분석이나 방법론을 사용하는 경우 다음을 인용해주세요:
- 원본 CataPro 논문 (Wang et al.)
- 이 분석 저장소 (DOI가 있는 경우)

## 연락처

이 분석에 대한 질문:
- 프로젝트 저장소: https://github.com/jahyunlee00299/CataPro
- 연락처: jahyunlee00299@gmail.com
