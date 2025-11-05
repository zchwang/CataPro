# Sequence Comparison Analysis Setup Guide

## 현재 상태 요약

### ✅ 완료된 것

1. **Alignment 파일 파싱 성공**
   - 파일: `DWC0THCP211-alignment.aln`
   - 서열 수: **1,181**
   - Alignment 길이: **711**
   - 형식: CLUSTAL

2. **PDB 구조 파일**
   - 파일: `pdb3rfv_chainA.pdb`
   - Chain A: 265 residues
   - 서열: AMKRLLVTGAAGQLGRVMRERLAPMAEILRLADLSPLDPAGPNEECVQCD...

3. **분석 도구 준비 완료**
   - `sequence_comparison_analyzer.py` - 메인 분석 도구
   - `parse_alignment_simple.py` - Alignment 파싱
   - `run_sequence_comparison.py` - 실행 스크립트
   - 가이드 문서들

### ❓ 필요한 것

**기질별로 나누어진 CataPro 예측 파일:**

현재 있는 파일:
- `catapro_predict_score.csv` (전체, 5,905 sequences)

필요한 파일:
- `catapro_predict_score_NADplus.csv` (1,181 sequences)
- `catapro_predict_score_α-D-Glucuronate.csv` (1,181 sequences)
- `catapro_predict_score_α-L-Guluronate.csv` (1,181 sequences)
- `catapro_predict_score_β-D-Galacturonate.csv` (1,181 sequences)
- `catapro_predict_score_β-D-Mannuronate.csv` (1,181 sequences)

## 다음 단계

### Option 1: 전체 파일에서 SMILES로 구분

`catapro_predict_score.csv`를 SMILES 컬럼으로 기질별로 분리:

```python
import pandas as pd

# Load full data
df = pd.read_csv('catapro_predict_score.csv')

# SMILES for each substrate
smiles_map = {
    'NADplus': '[specific SMILES for NAD+]',
    'alpha_D_Glucuronate': '[specific SMILES]',
    # ... etc
}

# Split by SMILES
for substrate, smiles in smiles_map.items():
    substrate_df = df[df['smiles'] == smiles]
    substrate_df.to_csv(f'catapro_predict_score_{substrate}.csv', index=False)
```

### Option 2: 이미 분리된 파일이 있는 경우

파일 이름이 다를 수 있습니다. 확인:
```bash
ls -la inference/*.csv
```

### Option 3: Alignment ID와 CataPro 매칭

Alignment의 fasta_id와 CataPro의 fasta_id가 다를 수 있습니다:

**Alignment IDs 형식:**
```
lcl|Query_10001
lcl|Query_10002
...
```

**CataPro IDs 형식:**
```
sp|Q7CRQ0|URODH_AGRFC
sp|Q888H1|URODH_PSESM
...
```

ID 매핑 파일이 필요할 수 있습니다.

## Binding Site Positions

### PDB에서 추출 (자동)

PDB 파일에서 ligand 주변 residues를 자동으로 찾을 수 있습니다:

```python
# Biopython이 설치되어 있다면
from prepare_analysis_data import DataPreparator

prep = DataPreparator()
prep.load_pdb_structure('pdb3rfv_chainA.pdb')

# NAD 주변 5Å 이내 residues
binding_sites = prep.extract_binding_sites_from_pdb('NAD', distance_cutoff=5.0)
```

### 수동 정의 (권장)

Clustering에서 사용한 positions를 직접 입력:

```python
binding_sites = {
    84: 'NAD_ribose',      # 0-indexed (PDB 85)
    85: 'NAD_adenine',     # 0-indexed (PDB 86)
    122: 'substrate',      # 0-indexed (PDB 123)
    155: 'catalytic',      # 0-indexed (PDB 156)
    180: 'cofactor'        # 0-indexed (PDB 181)
}
```

## 완전한 실행 예시

### 1. 데이터 준비

```python
# Parse alignment
from parse_alignment_simple import parse_clustal_alignment

alignment_df = parse_clustal_alignment('DWC0THCP211-alignment.aln')
# -> 1,181 sequences with aligned_seq column
```

### 2. CataPro 데이터와 병합

```python
# Option A: If you have substrate-specific files
import pandas as pd

substrate_files = {
    'NADplus': 'catapro_predict_score_NADplus.csv',
    'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
    # ...
}

for substrate, file in substrate_files.items():
    catapro_df = pd.read_csv(file)
    merged = alignment_df.merge(catapro_df, on='fasta_id', how='inner')
    merged.to_csv(f'catapro_with_alignment_{substrate}.csv', index=False)
```

```python
# Option B: If IDs don't match
# Need mapping file: Query_10001 -> sp|Q7CRQ0|URODH_AGRFC
mapping = pd.read_csv('id_mapping.csv')  # fasta_id_alignment, fasta_id_catapro

merged = alignment_df.merge(mapping, left_on='fasta_id', right_on='fasta_id_alignment')
merged = merged.merge(catapro_df, left_on='fasta_id_catapro', right_on='fasta_id')
```

### 3. Sequence Comparison 분석 실행

```python
from sequence_comparison_analyzer import SequenceComparisonAnalyzer

# Define binding sites
binding_sites = {
    84: 'NAD_ribose',
    85: 'NAD_adenine',
    122: 'substrate',
    155: 'catalytic',
    180: 'cofactor'
}

# Initialize
analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=3)

# Load substrate files (now with aligned_seq column!)
substrate_files = {
    'NADplus': 'catapro_with_alignment_NADplus.csv',
    'alpha_D_Glucuronate': 'catapro_with_alignment_alpha_D_Glucuronate.csv',
    'alpha_L_Guluronate': 'catapro_with_alignment_alpha_L_Guluronate.csv',
    'beta_D_Galacturonate': 'catapro_with_alignment_beta_D_Galacturonate.csv',
    'beta_D_Mannuronate': 'catapro_with_alignment_beta_D_Mannuronate.csv'
}

analyzer.load_substrate_data(substrate_files)

# Extract positions
analyzer.extract_position_data(aligned_seq_col='aligned_seq')

# Analyze
aa_effects = analyzer.analyze_amino_acid_effects(min_sequences=5)
substrate_specific = analyzer.find_substrate_specific_residues()

# Save results
analyzer.generate_report()
analyzer.plot_amino_acid_effects()
```

## 체크리스트

- [ ] Alignment 파일 확인 (`DWC0THCP211-alignment.aln`) ✓
- [ ] PDB 파일 확인 (`pdb3rfv_chainA.pdb`) ✓
- [ ] CataPro 예측 파일 확인
  - [ ] NADplus
  - [ ] α-D-Glucuronate
  - [ ] α-L-Guluronate
  - [ ] β-D-Galacturonate
  - [ ] β-D-Mannuronate
- [ ] Alignment ID와 CataPro ID 매칭 확인
- [ ] Binding site positions 정의
- [ ] 데이터 병합 완료
- [ ] 분석 실행

## 문제 해결

### Q1: Alignment ID와 CataPro ID가 다름

**확인:**
```python
aln_df = parse_clustal_alignment('DWC0THCP211-alignment.aln')
print("Alignment IDs:", aln_df['fasta_id'].head())

catapro_df = pd.read_csv('catapro_predict_score.csv')
print("CataPro IDs:", catapro_df['fasta_id'].head())
```

**해결:**
- ID 매핑 파일 생성
- 또는 공통 identifier 찾기 (예: position in file)

### Q2: Binding site positions를 모름

**방법 1: PDB에서 추출**
- Ligand (NAD, substrate) 주변 residues

**방법 2: Clustering 결과 사용**
- Clustering에서 사용한 positions

**방법 3: 문헌 참고**
- UDH family의 알려진 binding sites

### Q3: SMILES로 기질 구분 방법

```python
import pandas as pd

df = pd.read_csv('catapro_predict_score.csv')

# Unique SMILES
unique_smiles = df['smiles'].unique()
print(f"Found {len(unique_smiles)} unique SMILES")

# Check counts
smiles_counts = df['smiles'].value_counts()
print(smiles_counts)

# Expected:
# SMILES_1    1181    # NADplus
# SMILES_2    1181    # α-D-Glucuronate
# SMILES_3    1181    # α-L-Guluronate
# SMILES_4    1181    # β-D-Galacturonate
# SMILES_5    1181    # β-D-Mannuronate
```

## 현재 파일 구조

```
inference/
├── DWC0THCP211-alignment.aln          # ✓ Alignment (1,181 seqs)
├── pdb3rfv_chainA.pdb                  # ✓ Structure
├── catapro_predict_score.csv           # ✓ All predictions (5,905)
│
├── sequence_comparison_analyzer.py     # ✓ Main analysis tool
├── parse_alignment_simple.py           # ✓ Alignment parser
├── run_sequence_comparison.py          # ✓ Run script
│
├── SEQUENCE_COMPARISON_GUIDE.md        # ✓ Usage guide
└── SETUP_GUIDE.md                      # ✓ This file
```

## 다음에 필요한 것

1. **기질별 CataPro 파일 생성** (SMILES로 분리)
2. **Alignment와 CataPro 병합** (ID 매칭)
3. **Binding site positions 정의** (clustering 또는 PDB)

이 3가지가 준비되면 바로 분석 실행 가능!

---

**작성일**: 2025-11-04
