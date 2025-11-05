# Documentation

This directory contains documentation and guides for the UDH enzyme analysis project.

## Files

### Analysis Summaries

#### `ANALYSIS_COMPARISON_SUMMARY.md` (English)
Comprehensive summary of the mutation analysis validation study:
- Background on amino acid analysis approach
- CataPro predictions methodology
- Correlation analysis results
- Mutation validation results
- Critical findings and conclusions

#### `분석_비교_요약_한글.md` (Korean)
Korean translation of the analysis comparison summary for local team collaboration.

### Methodology Guides

#### `BINDING_SITE_ANALYSIS_GUIDE.md`
Guide for binding site mutation analysis:
- Methodology for identifying key residues
- Statistical approaches for optimal amino acid identification
- Mutation suggestion algorithm
- Scoring system for enzyme engineering

#### `SEQUENCE_COMPARISON_GUIDE.md`
Guide for sequence comparison analysis:
- Multiple sequence alignment processing
- PDB position mapping methodology
- Active site residue identification
- Comparison across enzyme variants

#### `POSITION_SELECTION_SUMMARY.md`
Summary of position selection criteria:
- Active site positions analyzed
- Rationale for position selection
- Structural basis for each position
- Position importance ranking

#### `SETUP_GUIDE.md`
Setup instructions for the analysis pipeline:
- Environment setup
- Required dependencies
- Data preparation steps
- Running the analysis pipeline

## Key Concepts

### Amino Acid Analysis
Statistical approach to identify "optimal" amino acids at active site positions:
1. Analyze 1181 UDH sequences from multiple sequence alignment
2. Identify most common amino acids at each position
3. Score enzymes based on matches to optimal residues
4. Suggest mutations to improve enzyme performance

### CataPro Predictions
Machine learning approach using protein and molecule language models:
- Predicts kcat, Km, and kcat/Km from sequence and substrate
- Trained on experimental enzyme kinetics data
- Uses ProtT5 for protein encoding and MolT5 for substrate encoding

### Critical Discovery

**NEGATIVE CORRELATION FOUND**: The amino acid analysis approach produces rankings that are inversely correlated (r = -0.573, p=0.032) with CataPro predictions.

**EMPIRICAL VALIDATION**: Creating actual mutant sequences with suggested mutations and testing them confirmed that 79% of mutants show decreased performance compared to wild-type.

### Implications

1. **Statistical analysis limitations**: Most common ≠ optimal for performance
2. **Evolutionary context matters**: High-performing enzymes may have unique residues
3. **Machine learning advantage**: CataPro captures complex sequence-function relationships
4. **Engineering strategy**: Avoid naive amino acid frequency-based approaches

## Usage

These documents provide context and methodology for understanding the analysis results in:
- `../reported_udh_analysis/` - Initial enzyme analysis
- `../correlation_analysis/` - Correlation discovery
- `../mutation_analysis/` - Empirical validation

## References

- BsUdh structure (PDB: 2F6Y) used as reference for active site positions
- Multiple sequence alignment of 1181 UDH sequences
- 14 experimentally characterized UDH enzymes for validation
