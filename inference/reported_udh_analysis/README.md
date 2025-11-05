# Reported UDH Analysis

This directory contains analysis of 14 experimentally reported UDH enzymes, comparing amino acid analysis with CataPro predictions.

## Key Scripts

### `analyze_reported_udhs.py`
Main analysis script that:
- Analyzes active site residues for 14 reported UDH enzymes
- Maps PDB positions to alignment positions
- Identifies "optimal" amino acids at key positions
- Calculates AA scores and suggests mutations
- Generates comprehensive Excel report

### `prepare_reported_udh_for_catapro.py`
Prepares CataPro input from reported enzyme sequences:
- Loads sequences from FASTA file
- Formats for CataPro prediction
- Creates input CSV with enzyme-substrate pairs

### `visualize_reported_udh_rankings.py`
Creates visualizations comparing rankings:
- AA analysis rankings vs CataPro rankings
- Multiple visualization types: scatter plots, bar charts, heatmaps
- Overall ranking comparison

## Data Files

### Input Files
- `reported_UDHs_complete.fasta` - Sequences of 14 reported enzymes
- `DWC0THCP211-alignment.aln` - Multiple sequence alignment

### Output Files
- `reported_udh_analysis.xlsx` - Comprehensive Excel report with:
  - Overall_Summary: Enzyme rankings by AA score
  - Active_Site_Residues: Detailed residue-by-residue analysis
  - Engineering_Suggestions: Mutation recommendations
- `reported_udh_catapro_input.csv` - Formatted for CataPro prediction
- `reported_udh_catapro_predictions.csv` - CataPro predicted kinetic parameters
- `reported_udh_rankings_summary.csv` - Comparison of rankings

### Visualizations
- `reported_udh_rankings_scatter.png` - Scatter plot of AA score vs CataPro
- `reported_udh_rankings_bar_plots.png` - Bar chart comparison
- `reported_udh_rankings_heatmap.png` - Heatmap of residue patterns
- `reported_udh_rankings_overall_ranking.png` - Ranking comparison

## Enzyme List

The 14 reported UDH enzymes analyzed:
1. AtUdh - Agrobacterium tumefaciens
2. SvUdh - Sinorhizobium violaceuinigrum
3. m1Udh - Mutant 1
4. OgUdh - Oceanicola granulosus
5. m2Udh - Mutant 2 (best CataPro performer)
6. TbUdh - Thermoanaerobacter brockii
7. PsUdh - Pseudomonas syringae
8. TtUdh - Thermus thermophilus
9. RmUdh - Rhizobium meliloti
10. BsUdh - Bacillus subtilis
11. PsaUdh - Pseudomonas aeruginosa
12. SvdUdh - Streptomyces violaceoruber
13. MtUdh - Mycobacterium tuberculosis
14. PnUdh - Pseudomonas nitroreducens

## Key Positions Analyzed

Active site positions based on PDB 2F6Y (BsUdh):
- K95, K96, K135 - Catalytic lysines
- W174, W175 - Substrate binding
- E140 - Proton transfer
- Y197, Y198, Y221 - NAD+ binding
- T130, T131, S132 - Substrate specificity
- L106, L107, L166 - Hydrophobic interactions

## Key Findings

### AA Analysis Rankings
Top performers by AA score:
1. AtUdh (AA score: 94.05)
2. SvUdh (93.29)
3. m1Udh (92.06)

### CataPro Rankings
Top performers by predicted kcat/Km:
1. m2Udh (predicted: 4.19)
2. PnUdh (4.12)
3. TbUdh (3.94)

### Discrepancy
**CRITICAL**: The two ranking methods show negative correlation (-0.573), suggesting that AA analysis-based rankings do not predict actual enzyme performance.

## Usage

To run the full analysis pipeline:

```bash
# 1. Analyze reported UDH sequences
python analyze_reported_udhs.py

# 2. Prepare for CataPro prediction
python prepare_reported_udh_for_catapro.py

# 3. Run CataPro predictions
python predict.py -inp_fpath reported_udh_catapro_input.csv \
                  -model_dpath ../models \
                  -batch_size 32 \
                  -device cpu \
                  -out_fpath reported_udh_catapro_predictions.csv

# 4. Visualize rankings
python visualize_reported_udh_rankings.py
```

## Conclusion

This analysis established the baseline comparison between amino acid analysis and CataPro predictions, revealing the negative correlation that was later confirmed through mutation testing (see `../mutation_analysis/`).
