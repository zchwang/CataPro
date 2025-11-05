# CataPro Inference - UDH Enzyme Analysis

This directory contains comprehensive analysis of UDH (Uronate Dehydrogenase) enzymes, comparing traditional amino acid analysis with machine learning predictions.

## Project Structure

```
inference/
├── mutation_analysis/       # Mutation validation study
├── correlation_analysis/    # AA analysis vs CataPro correlation
├── reported_udh_analysis/   # Analysis of 14 reported enzymes
├── documentation/          # Guides and documentation
├── mutation_analysis/      # Cluster analysis results (if present)
└── [core files]           # CataPro prediction scripts
```

## Analysis Overview

### 1. Reported UDH Analysis (`reported_udh_analysis/`)
**Objective**: Analyze 14 experimentally characterized UDH enzymes using amino acid frequency analysis

**Methods**:
- Multiple sequence alignment of 1181 UDH sequences
- Statistical identification of "optimal" amino acids at 15 active site positions
- Scoring and ranking based on matches to optimal residues
- Mutation suggestions to improve enzyme performance

**Key Output**: `reported_udh_analysis.xlsx` with enzyme rankings and mutation recommendations

### 2. Correlation Analysis (`correlation_analysis/`)
**Objective**: Compare AA analysis rankings with CataPro ML predictions

**Methods**:
- CataPro predictions for all 14 enzymes with NAD+ cofactor
- Pearson and Spearman correlation analysis
- Rank comparison visualization

**Key Finding**: **NEGATIVE CORRELATION** (r = -0.573, p=0.032)
- Enzymes with higher AA scores have LOWER predicted performance
- Enzymes with more optimal AA matches perform WORSE
- Enzymes needing more mutations perform BETTER

**Implication**: Statistical amino acid analysis may produce counter-productive engineering suggestions

### 3. Mutation Validation (`mutation_analysis/`)
**Objective**: Empirically test whether suggested mutations actually improve performance

**Methods**:
- Created 280 mutant sequences (14 enzymes × 3 mutation sets × 5 substrates + wild-types)
- Applied top5, top10, and all suggested mutations
- Predicted performance using CataPro for all mutants
- Compared mutant vs wild-type performance

**Key Results**:
- **Top5 mutations**: 71.4% degraded (Δkcat/Km = -0.0476)
- **Top10 mutations**: 71.4% degraded (Δkcat/Km = -0.1032)
- **All mutations**: 94.3% degraded (Δkcat/Km = -0.6143)
- **Overall**: 79% of mutants performed worse than wild-type

**Conclusion**: **HYPOTHESIS CONFIRMED** - Suggested mutations actually degrade enzyme performance

## Critical Findings

### Main Discovery
Traditional amino acid frequency analysis produces mutation suggestions that are **counter-productive** for enzyme engineering:

1. **Statistical frequency ≠ optimal performance**: Most common amino acids at a position do not necessarily lead to better catalytic efficiency

2. **Evolutionary context matters**: High-performing enzymes (e.g., m2Udh) may have unique residues that deviate from the consensus but are optimized for their specific function

3. **Machine learning captures complexity**: CataPro's deep learning approach captures complex sequence-function relationships that simple statistical analysis misses

4. **Engineering implications**: Enzyme engineering strategies should avoid naive frequency-based approaches

### Per-Enzyme Examples

**m2Udh** (CataPro #1, AA analysis #14):
- Ranked LAST by AA analysis (lowest score)
- Ranked FIRST by CataPro (best predicted performance)
- All suggested mutations degraded performance (-0.01 to -0.24)
- **Interpretation**: Already optimized, deviates from "consensus" for good reasons

**AtUdh** (AA analysis #1, CataPro #7):
- Ranked FIRST by AA analysis (highest score)
- Ranked 7th by CataPro (middle performance)
- All suggested mutations degraded performance (-0.11 to -0.29)
- **Interpretation**: High consensus match doesn't predict actual performance

## Core CataPro Files

### Prediction Scripts
- `predict.py` - Main CataPro prediction script
- `model.py` - Model architecture (ProtT5 + MolT5)
- `act_model.py` - Activation model components
- `utils.py` - Utility functions

### Analysis Support Scripts
- `parse_alignment_simple.py` - Parse CLUSTAL alignments
- `sequence_comparison_analyzer.py` - Sequence comparison utilities
- `binding_site_mutation_analyzer.py` - Mutation analysis tools
- `prepare_analysis_data.py` - Data preparation utilities
- `run_sequence_comparison.py` - Run sequence comparisons
- `run_analysis_with_existing_data.py` - Analysis pipeline

### Sample Data
- `sample_inp.csv` - Example input format for CataPro
- `catapro_predict_score.csv` - Large-scale prediction results
- `catapro_test-pred.csv` - Test predictions

## Usage

### Run Complete Analysis Pipeline

```bash
# 1. Analyze reported UDH sequences
cd reported_udh_analysis
python analyze_reported_udhs.py

# 2. Run CataPro predictions
python prepare_reported_udh_for_catapro.py
python ../predict.py -inp_fpath reported_udh_catapro_input.csv \
                     -model_dpath ../../models \
                     -batch_size 32 \
                     -device cpu \
                     -out_fpath reported_udh_catapro_predictions.csv

# 3. Visualize rankings
python visualize_reported_udh_rankings.py

# 4. Correlation analysis
cd ../correlation_analysis
python compare_mutation_impact_en.py
python plot_correlation_english.py

# 5. Create and test mutants
cd ../mutation_analysis
python create_mutant_sequences.py
python ../predict.py -inp_fpath mutant_sequences_catapro_input.csv \
                     -model_dpath ../../models \
                     -batch_size 32 \
                     -device cpu \
                     -out_fpath mutant_sequences_predictions.csv

# 6. Analyze mutation impact
python analyze_mutant_impact.py
```

### Run Single CataPro Prediction

```bash
python predict.py \
    -inp_fpath input.csv \
    -model_dpath ../models \
    -batch_size 32 \
    -device cpu \
    -out_fpath output.csv
```

Input CSV format:
```csv
,Enzyme_id,type,sequence,smiles
0,enzyme1,wild-type,MKLV...,C1=CC=C...
```

## Key Takeaways

1. **For Researchers**:
   - Be cautious with statistical amino acid analysis for enzyme engineering
   - Validate computational predictions empirically
   - Consider evolutionary context and structural information

2. **For Engineers**:
   - Machine learning approaches (like CataPro) capture more complexity than simple statistics
   - Counter-intuitive residues may be functionally important
   - Engineering strategies should integrate multiple information sources

3. **For Method Developers**:
   - Frequency-based analysis has fundamental limitations
   - Negative validation studies are valuable for the field
   - ML models can provide mechanistic insights through prediction patterns

## References

- **CataPro**: Wang et al., enzyme kinetic parameter prediction using protein and molecule language models
- **PDB Structure**: 2F6Y (BsUdh) used as reference for active site positions
- **Dataset**: Multiple sequence alignment of 1181 UDH enzyme sequences
- **Validation Set**: 14 experimentally characterized UDH enzymes

## Citation

If you use this analysis or methodology, please cite:
- Original CataPro publication (Wang et al.)
- This analysis repository with DOI (if available)

## Contact

For questions about this analysis:
- Project Repository: [GitHub URL]
- Contact: [Your contact information]
