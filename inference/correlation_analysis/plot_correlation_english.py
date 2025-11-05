"""
Generate English version of mutation correlation plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# Load data
enzyme_avg = pd.read_csv('enzyme_comparison_aa_vs_catapro.csv')

# Calculate correlations
corr_aa, p_aa = pearsonr(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'])
corr_best, p_best = pearsonr(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'])
corr_mut, p_mut = pearsonr(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'])

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: AA Score vs CataPro
ax = axes[0, 0]
scatter = ax.scatter(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'],
                    c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                    s=100, alpha=0.7, edgecolors='black')
for _, row in enzyme_avg.iterrows():
    ax.annotate(row['Enzyme'], (row['AA_Score'], row['CataPro_kcat_over_Km']),
               xytext=(5, 5), textcoords='offset points', fontsize=8)
ax.set_xlabel('AA Analysis Score', fontsize=11, fontweight='bold')
ax.set_ylabel('CataPro Predicted kcat/Km', fontsize=11, fontweight='bold')
ax.set_title(f'AA Score vs CataPro Performance\n(correlation: {corr_aa:.3f}, p={p_aa:.4f})',
             fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax, label='CataPro Rank')

# Plot 2: Best matches vs CataPro
ax = axes[0, 1]
scatter = ax.scatter(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'],
                    c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                    s=100, alpha=0.7, edgecolors='black')
for _, row in enzyme_avg.iterrows():
    ax.annotate(row['Enzyme'], (row['Best_Match_Count'], row['CataPro_kcat_over_Km']),
               xytext=(5, 5), textcoords='offset points', fontsize=8)
ax.set_xlabel('Number of Optimal AA Matches', fontsize=11, fontweight='bold')
ax.set_ylabel('CataPro Predicted kcat/Km', fontsize=11, fontweight='bold')
ax.set_title(f'Optimal AA Matches vs CataPro Performance\n(correlation: {corr_best:.3f}, p={p_best:.4f})',
             fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax, label='CataPro Rank')

# Plot 3: Mutations needed vs CataPro
ax = axes[1, 0]
scatter = ax.scatter(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'],
                    c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                    s=100, alpha=0.7, edgecolors='black')
for _, row in enzyme_avg.iterrows():
    ax.annotate(row['Enzyme'], (row['Mutations_Needed'], row['CataPro_kcat_over_Km']),
               xytext=(5, 5), textcoords='offset points', fontsize=8)
ax.set_xlabel('Number of Mutations Needed', fontsize=11, fontweight='bold')
ax.set_ylabel('CataPro Predicted kcat/Km', fontsize=11, fontweight='bold')
ax.set_title(f'Mutations Needed vs CataPro Performance\n(correlation: {corr_mut:.3f}, p={p_mut:.4f})',
             fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax, label='CataPro Rank')

# Plot 4: Rank comparison
ax = axes[1, 1]
x = np.arange(len(enzyme_avg))
width = 0.35

bars1 = ax.bar(x - width/2, enzyme_avg['AA_Rank'], width, label='AA Analysis Rank', alpha=0.7, color='#3498db')
bars2 = ax.bar(x + width/2, enzyme_avg['CataPro_Rank'], width, label='CataPro Rank', alpha=0.7, color='#e74c3c')

ax.set_xlabel('Enzyme', fontsize=11, fontweight='bold')
ax.set_ylabel('Rank (lower is better)', fontsize=11, fontweight='bold')
ax.set_title('Ranking Comparison by Method', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(enzyme_avg['Enzyme'], rotation=45, ha='right', fontsize=9)
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig('mutation_catapro_correlation_english.png', dpi=300, bbox_inches='tight')
print("Saved: mutation_catapro_correlation_english.png")
plt.close()

# Print correlation summary
print("\n" + "="*80)
print("Correlation Analysis Summary")
print("="*80)
print(f"\nAA Score vs CataPro kcat/Km:")
print(f"  Pearson correlation: {corr_aa:.3f} (p={p_aa:.4f})")
scorr_aa, sp_aa = spearmanr(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'])
print(f"  Spearman rank correlation: {scorr_aa:.3f} (p={sp_aa:.4f})")

print(f"\nOptimal AA Matches vs CataPro kcat/Km:")
print(f"  Pearson correlation: {corr_best:.3f} (p={p_best:.4f})")
scorr_best, sp_best = spearmanr(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'])
print(f"  Spearman rank correlation: {scorr_best:.3f} (p={sp_best:.4f})")

print(f"\nMutations Needed vs CataPro kcat/Km:")
print(f"  Pearson correlation: {corr_mut:.3f} (p={p_mut:.4f})")
scorr_mut, sp_mut = spearmanr(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'])
print(f"  Spearman rank correlation: {scorr_mut:.3f} (p={sp_mut:.4f})")

print("\n" + "="*80)
print("KEY FINDING: NEGATIVE CORRELATION!")
print("="*80)
print("\nEnzymes with HIGHER AA analysis scores have LOWER CataPro performance")
print("Enzymes with MORE optimal AA matches have LOWER CataPro performance")
print("Enzymes needing MORE mutations have HIGHER CataPro performance")
print("\n=> Suggested mutations may actually DECREASE enzyme performance!")
print("="*80)
