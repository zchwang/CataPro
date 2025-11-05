# 상관관계 분석

이 디렉토리는 통계적 AA 분석과 실제 효소 성능 간의 관계를 식별하기 위해 아미노산 분석 점수와 CataPro 예측을 비교하는 분석을 포함합니다.

## 주요 스크립트

### `compare_mutation_impact.py`
초기 상관관계 분석 (한글 라벨):
- AA 분석 점수와 CataPro 예측 성능 비교
- Pearson 및 Spearman 상관관계 계산
- 한글 라벨로 시각화 생성

### `compare_mutation_impact_en.py`
추가 분석이 포함된 영어 버전:
- 위와 동일한 분석이지만 영어 라벨 사용
- 더 상세한 통계 보고

### `plot_correlation_english.py`
독립형 시각화 스크립트:
- 깔끔한 영어 라벨 상관관계 플롯 생성
- 다양한 상관관계 측면을 보여주는 4개 패널 시각화

## 데이터 파일

### 입력/출력
- `enzyme_comparison_aa_vs_catapro.csv` - AA 점수 및 CataPro 예측이 포함된 비교 데이터

### 시각화
- `mutation_catapro_correlation.png` - 한글 버전 (폰트 문제로 더 이상 사용 안 함)
- `mutation_catapro_correlation_english.png` - 깔끔한 영어 버전

## 주요 발견사항

### 발견된 음의 상관관계

**AA 점수 vs CataPro kcat/Km**:
- Pearson 상관관계: **-0.573** (p=0.032)
- Spearman 순위 상관관계: **-0.560** (p=0.037)
- **해석**: AA 분석 점수가 높은 효소는 CataPro 성능이 낮음

**최적 AA 일치 vs CataPro kcat/Km**:
- Pearson 상관관계: **-0.505** (p=0.067)
- Spearman 순위 상관관계: **-0.704** (p=0.005)
- **해석**: 최적 AA 일치가 많은 효소는 CataPro 성능이 낮음

**필요한 돌연변이 vs CataPro kcat/Km**:
- Pearson 상관관계: **+0.557** (p=0.038)
- Spearman 순위 상관관계: **+0.701** (p=0.005)
- **해석**: 더 많은 돌연변이가 필요한 효소는 CataPro 성능이 높음

## 시각화 세부사항

4개 패널 플롯은 다음을 보여줍니다:
1. **AA 점수 vs CataPro 성능** - 음의 추세
2. **최적 AA 일치 vs CataPro** - 강한 음의 순위 상관관계
3. **필요한 돌연변이 vs CataPro** - 양의 상관관계
4. **순위 비교** - AA 분석 순위와 CataPro 순위가 크게 다름

## 핵심 결론

**음의 상관관계 발견**: 통계적 아미노산 분석 접근법은 실제 효소 성능 예측과 역상관 관계에 있는 순위를 생성합니다. 이는 다음을 시사합니다:

1. 진화에 의해 자연적으로 최적화된 효소는 통계 분석의 "최적" 잔기와 일치하지 않을 수 있음
2. 고성능 효소(m2Udh 등)는 "최적" 일치가 적음
3. AA 점수를 높이기 위해 돌연변이를 적용하면 성능이 저하될 가능성이 높음

이 발견은 돌연변이 분석 검증 연구로 이어졌습니다 (`../mutation_analysis/` 참조).
