# Correlation Analysis

This directory contains analysis comparing amino acid analysis scores with CataPro predictions to identify relationships between statistical AA analysis and actual enzyme performance.

## Key Scripts

### `compare_mutation_impact.py`
Initial correlation analysis (Korean labels):
- Compares AA analysis scores with CataPro predicted performance
- Calculates Pearson and Spearman correlations
- Generates visualization with Korean labels

### `compare_mutation_impact_en.py`
English version with additional analysis:
- Same analysis as above but with English labels
- More detailed statistical reporting

### `plot_correlation_english.py`
Standalone visualization script:
- Generates clean English-labeled correlation plots
- 4-panel visualization showing different correlation aspects

## Data Files

### Input/Output
- `enzyme_comparison_aa_vs_catapro.csv` - Comparison data with AA scores and CataPro predictions

### Visualizations
- `mutation_catapro_correlation.png` - Korean version (deprecated, font issues)
- `mutation_catapro_correlation_english.png` - Clean English version

## Key Findings

### Negative Correlations Found

**AA Score vs CataPro kcat/Km**:
- Pearson correlation: **-0.573** (p=0.032)
- Spearman rank correlation: **-0.560** (p=0.037)
- **Interpretation**: Enzymes with HIGHER AA analysis scores have LOWER CataPro performance

**Optimal AA Matches vs CataPro kcat/Km**:
- Pearson correlation: **-0.505** (p=0.067)
- Spearman rank correlation: **-0.704** (p=0.005)
- **Interpretation**: Enzymes with MORE optimal AA matches have LOWER CataPro performance

**Mutations Needed vs CataPro kcat/Km**:
- Pearson correlation: **+0.557** (p=0.038)
- Spearman rank correlation: **+0.701** (p=0.005)
- **Interpretation**: Enzymes needing MORE mutations have HIGHER CataPro performance

## Visualization Details

The 4-panel plot shows:
1. **AA Score vs CataPro Performance** - Negative trend
2. **Optimal AA Matches vs CataPro** - Strong negative rank correlation
3. **Mutations Needed vs CataPro** - Positive correlation
4. **Ranking Comparison** - AA analysis ranks vs CataPro ranks differ significantly

## Critical Conclusion

**NEGATIVE CORRELATION DETECTED**: The statistical amino acid analysis approach produces rankings that are inversely correlated with actual enzyme performance predictions. This suggests that:

1. Enzymes naturally optimized by evolution may not match "optimal" residues from statistical analysis
2. High-performing enzymes (like m2Udh) have fewer "optimal" matches
3. Applying mutations to increase AA score likely degrades performance

This finding led to the mutation analysis validation study (see `../mutation_analysis/`).
