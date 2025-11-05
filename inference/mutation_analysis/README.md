# Mutation Analysis

This directory contains scripts and results for testing the hypothesis that amino acid analysis-suggested mutations improve enzyme performance.

## Key Scripts

### `create_mutant_sequences.py`
Creates mutant enzyme sequences by applying suggested mutations from AA analysis:
- Loads original sequences from `reported_UDHs_complete.fasta`
- Loads mutation suggestions from `reported_udh_analysis.xlsx`
- Creates three mutation sets: top5, top10, all
- Generates sequences for 14 enzymes × 5 substrates × 4 types (WT + 3 mutant sets)

### `analyze_mutant_impact.py`
Analyzes the performance difference between wild-type and mutant enzymes:
- Compares CataPro predictions for wild-type vs mutants
- Calculates Δkcat/Km for each mutant
- Generates statistical summaries and visualizations

## Data Files

### Input Files
- `mutant_sequences_catapro_input.csv` - 280 sequences for CataPro prediction (70 WT + 210 mutants)
- `mutant_sequences_metadata.csv` - Metadata tracking enzyme, substrate, mutation type

### Output Files
- `mutant_sequences_predictions.csv` - CataPro predictions for all 280 sequences
- `mutant_sequences_summary.csv` - Summary of mutations applied per enzyme
- `mutant_vs_wildtype_comparison.csv` - Detailed comparison (210 mutant-WT pairs)

### Visualizations
- `mutant_impact_analysis.png` - 4-panel visualization showing:
  - Distribution of performance changes
  - Wild-type vs Mutant scatter plot
  - Average change by enzyme
  - Performance change vs number of mutations

## Key Findings

**CRITICAL RESULT**: Mutations suggested by amino acid analysis actually **DEGRADE** enzyme performance in 79% of cases.

### Overall Performance Changes
- **Top5 mutations**: Average Δkcat/Km = -0.0476 (71.4% degraded)
- **Top10 mutations**: Average Δkcat/Km = -0.1032 (71.4% degraded)
- **All mutations**: Average Δkcat/Km = -0.6143 (94.3% degraded)

### Per-Enzyme Results
- **m2Udh** (CataPro #1): All mutations degraded performance (-0.01 to -0.24)
- **AtUdh** (AA #1): All mutations degraded performance (-0.11 to -0.29)
- **SvUdh**: Only enzyme with slight improvement from "all" mutations (+0.03)

## Conclusion

The empirical validation confirms the negative correlation hypothesis: applying mutations to reach "optimal" amino acids as identified by statistical analysis actually decreases enzyme performance on average.
