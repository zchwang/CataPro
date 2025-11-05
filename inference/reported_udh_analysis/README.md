# 보고된 UDH 분석

이 디렉토리는 14개의 실험적으로 보고된 UDH 효소의 분석을 포함하며, 아미노산 분석과 CataPro 예측을 비교합니다.

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

## 데이터 파일

### 입력 파일
- `reported_UDHs_complete.fasta` - 14개 보고된 효소의 서열
- `DWC0THCP211-alignment.aln` - 다중 서열 정렬

### 출력 파일
- `reported_udh_analysis.xlsx` - 포괄적인 Excel 보고서:
  - Overall_Summary: AA 점수별 효소 순위
  - Active_Site_Residues: 잔기별 상세 분석
  - Engineering_Suggestions: 돌연변이 권장사항
- `reported_udh_catapro_input.csv` - CataPro 예측을 위한 형식 지정
- `reported_udh_catapro_predictions.csv` - CataPro가 예측한 동역학 매개변수
- `reported_udh_rankings_summary.csv` - 순위 비교

### 시각화
- `reported_udh_rankings_scatter.png` - AA 점수 vs CataPro 산점도
- `reported_udh_rankings_bar_plots.png` - 막대 차트 비교
- `reported_udh_rankings_heatmap.png` - 잔기 패턴 히트맵
- `reported_udh_rankings_overall_ranking.png` - 순위 비교

## 효소 목록

분석된 14개의 보고된 UDH 효소:
1. AtUdh - Agrobacterium tumefaciens
2. SvUdh - Sinorhizobium violaceuinigrum
3. m1Udh - Mutant 1
4. OgUdh - Oceanicola granulosus
5. m2Udh - Mutant 2 (CataPro 최고 성능)
6. TbUdh - Thermoanaerobacter brockii
7. PsUdh - Pseudomonas syringae
8. TtUdh - Thermus thermophilus
9. RmUdh - Rhizobium meliloti
10. BsUdh - Bacillus subtilis
11. PsaUdh - Pseudomonas aeruginosa
12. SvdUdh - Streptomyces violaceoruber
13. MtUdh - Mycobacterium tuberculosis
14. PnUdh - Pseudomonas nitroreducens

## 분석된 주요 위치

PDB 2F6Y (BsUdh) 기반 활성 부위 위치:
- K95, K96, K135 - 촉매 라이신
- W174, W175 - 기질 결합
- E140 - 양성자 전달
- Y197, Y198, Y221 - NAD+ 결합
- T130, T131, S132 - 기질 특이성
- L106, L107, L166 - 소수성 상호작용

## 주요 발견사항

### AA 분석 순위
AA 점수별 상위 성능:
1. AtUdh (AA 점수: 94.05)
2. SvUdh (93.29)
3. m1Udh (92.06)

### CataPro 순위
예측된 kcat/Km별 상위 성능:
1. m2Udh (예측값: 4.19)
2. PnUdh (4.12)
3. TbUdh (3.94)

### 불일치
**핵심**: 두 순위 방법은 음의 상관관계 (-0.573)를 보여, AA 분석 기반 순위가 실제 효소 성능을 예측하지 못함을 시사합니다.

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
```

## 결론

이 분석은 아미노산 분석과 CataPro 예측 간의 기준 비교를 확립하여, 나중에 돌연변이 테스트를 통해 확인된 음의 상관관계를 밝혀냈습니다 (`../mutation_analysis/` 참조).
