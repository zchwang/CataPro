"""
Analyze the impact of mutations on enzyme performance
Compare wild-type vs mutants predicted by CataPro
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_mutant_impact():
    """Analyze mutant predictions vs wild-type"""

    print("="*80)
    print("Analyzing Mutant Impact: Wild-type vs Engineered Variants")
    print("="*80)

    # Load predictions
    print("\n1. Loading predictions...")
    pred_df = pd.read_csv('mutant_sequences_predictions.csv', index_col=0)
    metadata_df = pd.read_csv('mutant_sequences_metadata.csv')

    # Rename columns
    pred_df = pred_df.rename(columns={
        'pred_log10[kcat(s^-1)]': 'kcat',
        'pred_log10[Km(mM)]': 'Km',
        'pred_log10[kcat/Km(s^-1mM^-1)]': 'kcat_over_Km'
    })

    # Fix fasta_id: CataPro appends "_wild-type" or "_mutant" to the enzyme_id
    # We need to remove it to match with metadata
    pred_df['Enzyme_id_clean'] = pred_df['fasta_id'].str.replace('_wild-type$', '', regex=True).str.replace('_mutant$', '', regex=True)

    # Merge with metadata
    df = pred_df.merge(metadata_df[['Enzyme_id', 'original_enzyme', 'mutant_type',
                                     'substrate', 'num_mutations']],
                       left_on='Enzyme_id_clean', right_on='Enzyme_id', how='left')

    print(f"   Total predictions: {len(df)}")

    # Separate wild-type and mutants
    wt_df = df[df['mutant_type'] == 'wildtype'].copy()
    mut_df = df[df['mutant_type'] != 'wildtype'].copy()

    print(f"   Wild-type: {len(wt_df)}")
    print(f"   Mutants: {len(mut_df)}")

    # === 2. Calculate changes ===
    print("\n2. Calculating performance changes...")

    comparison_data = []

    for _, mut_row in mut_df.iterrows():
        enzyme = mut_row['original_enzyme']
        substrate = mut_row['substrate']
        mutant_type = mut_row['mutant_type']

        # Find corresponding wild-type
        wt_row = wt_df[(wt_df['original_enzyme'] == enzyme) &
                       (wt_df['substrate'] == substrate)]

        if len(wt_row) == 0:
            continue

        wt_row = wt_row.iloc[0]

        # Calculate change
        delta_kcat_over_Km = mut_row['kcat_over_Km'] - wt_row['kcat_over_Km']
        delta_kcat = mut_row['kcat'] - wt_row['kcat']
        delta_Km = mut_row['Km'] - wt_row['Km']

        comparison_data.append({
            'Enzyme': enzyme,
            'Substrate': substrate,
            'Mutant_Type': mutant_type,
            'Num_Mutations': mut_row['num_mutations'],
            'WT_kcat_over_Km': wt_row['kcat_over_Km'],
            'Mutant_kcat_over_Km': mut_row['kcat_over_Km'],
            'Delta_kcat_over_Km': delta_kcat_over_Km,
            'WT_kcat': wt_row['kcat'],
            'Mutant_kcat': mut_row['kcat'],
            'Delta_kcat': delta_kcat,
            'WT_Km': wt_row['Km'],
            'Mutant_Km': mut_row['Km'],
            'Delta_Km': delta_Km,
            'Improved': delta_kcat_over_Km > 0
        })

    comparison_df = pd.DataFrame(comparison_data)

    # Save detailed comparison
    comparison_df.to_csv('mutant_vs_wildtype_comparison.csv', index=False)
    print(f"   Saved detailed comparison to: mutant_vs_wildtype_comparison.csv")

    # === 3. Summary statistics ===
    print("\n" + "="*80)
    print("3. SUMMARY STATISTICS")
    print("="*80)

    for mutant_type in ['top5', 'top10', 'all']:
        type_df = comparison_df[comparison_df['Mutant_Type'] == mutant_type]

        if len(type_df) == 0:
            continue

        print(f"\n{mutant_type.upper()} Mutations:")
        print("-"*80)

        avg_delta = type_df['Delta_kcat_over_Km'].mean()
        improved_pct = (type_df['Improved'].sum() / len(type_df)) * 100
        degraded_pct = 100 - improved_pct

        print(f"  Average change in kcat/Km: {avg_delta:+.4f}")
        print(f"  Improved: {type_df['Improved'].sum()}/{len(type_df)} ({improved_pct:.1f}%)")
        print(f"  Degraded: {len(type_df) - type_df['Improved'].sum()}/{len(type_df)} ({degraded_pct:.1f}%)")

        # Best and worst cases
        best = type_df.nlargest(3, 'Delta_kcat_over_Km')
        worst = type_df.nsmallest(3, 'Delta_kcat_over_Km')

        print(f"\n  Best improvements:")
        for i, (_, row) in enumerate(best.iterrows(), 1):
            print(f"    {i}. {row['Enzyme']} ({row['Substrate']}): {row['Delta_kcat_over_Km']:+.4f}")

        print(f"\n  Worst degradations:")
        for i, (_, row) in enumerate(worst.iterrows(), 1):
            print(f"    {i}. {row['Enzyme']} ({row['Substrate']}): {row['Delta_kcat_over_Km']:+.4f}")

    # === 4. Per-enzyme analysis ===
    print("\n" + "="*80)
    print("4. PER-ENZYME ANALYSIS")
    print("="*80)

    # Focus on key enzymes
    key_enzymes = ['m2Udh', 'AtUdh', 'SvUdh', 'PnUdh']

    for enzyme in key_enzymes:
        enzyme_df = comparison_df[comparison_df['Enzyme'] == enzyme]

        if len(enzyme_df) == 0:
            continue

        print(f"\n{enzyme}:")
        print("-"*80)

        for mutant_type in ['top5', 'top10', 'all']:
            type_df = enzyme_df[enzyme_df['Mutant_Type'] == mutant_type]

            if len(type_df) == 0:
                continue

            avg_delta = type_df['Delta_kcat_over_Km'].mean()
            improved = (type_df['Delta_kcat_over_Km'] > 0).sum()

            print(f"  {mutant_type:6s}: Δkcat/Km = {avg_delta:+.4f}, Improved: {improved}/{len(type_df)}")

    # === 5. Visualizations ===
    print("\n" + "="*80)
    print("5. Generating visualizations...")
    print("="*80)

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Distribution of changes
    ax = axes[0, 0]
    for mutant_type, color in [('top5', '#3498db'), ('top10', '#e74c3c'), ('all', '#2ecc71')]:
        type_df = comparison_df[comparison_df['Mutant_Type'] == mutant_type]
        ax.hist(type_df['Delta_kcat_over_Km'], bins=30, alpha=0.5, label=mutant_type, color=color)

    ax.axvline(x=0, color='black', linestyle='--', linewidth=2, label='No change')
    ax.set_xlabel('Change in kcat/Km (log10)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Distribution of Performance Changes', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Plot 2: Wild-type vs Mutant (all mutations)
    ax = axes[0, 1]
    all_mut_df = comparison_df[comparison_df['Mutant_Type'] == 'all']

    scatter = ax.scatter(all_mut_df['WT_kcat_over_Km'], all_mut_df['Mutant_kcat_over_Km'],
                        c=all_mut_df['Delta_kcat_over_Km'], cmap='RdYlGn',
                        s=50, alpha=0.6, edgecolors='black', linewidth=0.5)

    # Draw diagonal (no change line)
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0, label='No change')

    ax.set_xlabel('Wild-type kcat/Km', fontsize=11, fontweight='bold')
    ax.set_ylabel('Mutant kcat/Km (all mutations)', fontsize=11, fontweight='bold')
    ax.set_title('Wild-type vs Mutant Performance', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Change in kcat/Km')

    # Plot 3: Change by enzyme
    ax = axes[1, 0]

    # Average change per enzyme for 'all' mutations
    all_mut_df = comparison_df[comparison_df['Mutant_Type'] == 'all']
    enzyme_avg = all_mut_df.groupby('Enzyme')['Delta_kcat_over_Km'].mean().sort_values()

    colors = ['red' if x < 0 else 'green' for x in enzyme_avg.values]
    bars = ax.barh(range(len(enzyme_avg)), enzyme_avg.values, color=colors, alpha=0.7)

    ax.axvline(x=0, color='black', linestyle='--', linewidth=2)
    ax.set_yticks(range(len(enzyme_avg)))
    ax.set_yticklabels(enzyme_avg.index)
    ax.set_xlabel('Average Change in kcat/Km', fontsize=11, fontweight='bold')
    ax.set_ylabel('Enzyme', fontsize=11, fontweight='bold')
    ax.set_title('Average Performance Change by Enzyme\n(All Mutations Applied)', fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    # Plot 4: Number of mutations vs change
    ax = axes[1, 1]

    for mutant_type, color, marker in [('top5', '#3498db', 'o'), ('top10', '#e74c3c', 's'), ('all', '#2ecc71', '^')]:
        type_df = comparison_df[comparison_df['Mutant_Type'] == mutant_type]
        ax.scatter(type_df['Num_Mutations'], type_df['Delta_kcat_over_Km'],
                  alpha=0.5, label=mutant_type, color=color, marker=marker, s=50)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=2)
    ax.set_xlabel('Number of Mutations', fontsize=11, fontweight='bold')
    ax.set_ylabel('Change in kcat/Km', fontsize=11, fontweight='bold')
    ax.set_title('Performance Change vs Number of Mutations', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('mutant_impact_analysis.png', dpi=300, bbox_inches='tight')
    print("   Saved: mutant_impact_analysis.png")
    plt.close()

    # === 6. KEY FINDINGS ===
    print("\n" + "="*80)
    print("6. KEY FINDINGS")
    print("="*80)

    overall_avg = comparison_df.groupby('Mutant_Type')['Delta_kcat_over_Km'].mean()

    print(f"\nOverall average performance changes:")
    for mutant_type in ['top5', 'top10', 'all']:
        if mutant_type in overall_avg.index:
            avg = overall_avg[mutant_type]
            sign = "+" if avg > 0 else ""
            print(f"  {mutant_type:6s}: {sign}{avg:.4f}")

    # Statistical test
    improved_count = (comparison_df['Delta_kcat_over_Km'] > 0).sum()
    degraded_count = (comparison_df['Delta_kcat_over_Km'] < 0).sum()
    total_count = len(comparison_df)

    print(f"\nOverall results:")
    print(f"  Improved: {improved_count}/{total_count} ({improved_count/total_count*100:.1f}%)")
    print(f"  Degraded: {degraded_count}/{total_count} ({degraded_count/total_count*100:.1f}%)")

    if overall_avg['all'] < 0:
        print("\n🔴 CRITICAL FINDING:")
        print("   Suggested mutations DECREASE enzyme performance on average!")
        print("   The amino acid analysis-based engineering strategy is COUNTER-PRODUCTIVE.")
    else:
        print("\n✓ Mutations improve performance on average")

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)

    return comparison_df


if __name__ == "__main__":
    comparison_df = analyze_mutant_impact()
