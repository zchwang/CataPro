"""
Visualize CataPro predictions for 14 reported UDH enzymes
Show rankings based on predicted kcat, Km, and kcat/Km
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

def parse_enzyme_substrate(fasta_id):
    """Parse fasta_id to get enzyme name and substrate
    Format: EnzymeName_Substrate_wild-type
    """
    # Remove '_wild-type' suffix
    fasta_id_clean = fasta_id.replace('_wild-type', '')

    # Known substrates
    substrates = ['NADplus', 'alpha_D_Glucuronate', 'alpha_L_Guluronate',
                  'beta_D_Galacturonate', 'beta_D_Mannuronate']

    # Find which substrate this is
    for substrate in substrates:
        if fasta_id_clean.endswith('_' + substrate):
            enzyme = fasta_id_clean[:-len('_' + substrate)]
            return enzyme, substrate

    # Fallback
    return fasta_id_clean, 'unknown'

def create_rankings_visualization(predictions_file, output_prefix='reported_udh_rankings'):
    """
    Create comprehensive visualization of enzyme rankings

    Parameters:
    -----------
    predictions_file : str
        Path to CataPro predictions CSV
    output_prefix : str
        Prefix for output files
    """

    # Load predictions
    print("Loading CataPro predictions...")
    df = pd.read_csv(predictions_file, index_col=0)

    # Rename columns for easier access
    df = df.rename(columns={
        'pred_log10[kcat(s^-1)]': 'kcat',
        'pred_log10[Km(mM)]': 'Km',
        'pred_log10[kcat/Km(s^-1mM^-1)]': 'kcat_over_Km'
    })

    # Parse enzyme and substrate from fasta_id
    df[['Enzyme', 'Substrate']] = df['fasta_id'].apply(
        lambda x: pd.Series(parse_enzyme_substrate(x))
    )

    print(f"Loaded {len(df)} predictions for {df['Enzyme'].nunique()} enzymes across {df['Substrate'].nunique()} substrates")

    # Substrate name mapping for display
    substrate_names = {
        'NADplus': 'NAD+',
        'alpha_D_Glucuronate': 'α-D-Glucuronate',
        'alpha_L_Guluronate': 'α-L-Guluronate',
        'beta_D_Galacturonate': 'β-D-Galacturonate',
        'beta_D_Mannuronate': 'β-D-Mannuronate'
    }

    df['Substrate_Display'] = df['Substrate'].map(substrate_names).fillna(df['Substrate'])

    # Get unique enzymes and substrates
    enzymes = sorted(df['Enzyme'].unique())
    substrates = ['NADplus', 'alpha_D_Glucuronate', 'alpha_L_Guluronate',
                  'beta_D_Galacturonate', 'beta_D_Mannuronate']

    # === 1. Heatmap: kcat/Km for all enzymes and substrates ===
    print("\nCreating kcat/Km heatmap...")
    fig, ax = plt.subplots(figsize=(12, 10))

    # Pivot to create matrix
    heatmap_data = df.pivot(index='Enzyme', columns='Substrate_Display', values='kcat_over_Km')

    # Reorder columns
    col_order = [substrate_names[s] for s in substrates if s in df['Substrate'].values]
    heatmap_data = heatmap_data[col_order]

    # Create heatmap
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='RdYlGn',
                cbar_kws={'label': 'log10(kcat/Km)'}, ax=ax)

    ax.set_title('Predicted Catalytic Efficiency (kcat/Km) for Reported UDH Enzymes',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Substrate', fontsize=12, fontweight='bold')
    ax.set_ylabel('Enzyme', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_heatmap.png")
    plt.close()

    # === 2. Bar plots: Rankings by substrate ===
    print("\nCreating ranking bar plots...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, substrate in enumerate(substrates):
        substrate_df = df[df['Substrate'] == substrate].copy()
        substrate_df = substrate_df.sort_values('kcat_over_Km', ascending=False)

        ax = axes[idx]
        bars = ax.barh(range(len(substrate_df)), substrate_df['kcat_over_Km'])

        # Color bars: top 3 green, bottom 3 red
        colors = ['green'] * 3 + ['yellow'] * (len(substrate_df) - 6) + ['red'] * 3
        for bar, color in zip(bars, colors):
            bar.set_color(color)
            bar.set_alpha(0.7)

        ax.set_yticks(range(len(substrate_df)))
        ax.set_yticklabels(substrate_df['Enzyme'], fontsize=9)
        ax.set_xlabel('log10(kcat/Km)', fontsize=10, fontweight='bold')
        ax.set_title(substrate_names[substrate], fontsize=11, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)

        # Add rank numbers
        for i, (idx_val, row) in enumerate(substrate_df.iterrows()):
            ax.text(row['kcat_over_Km'] + 0.05, i, f"#{i+1}",
                   va='center', fontsize=8, fontweight='bold')

    # Hide unused subplot
    axes[5].axis('off')

    fig.suptitle('Enzyme Rankings by Substrate (Predicted kcat/Km)',
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_bar_plots.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_bar_plots.png")
    plt.close()

    # === 3. Scatter plots: kcat vs Km ===
    print("\nCreating kcat vs Km scatter plots...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, substrate in enumerate(substrates):
        substrate_df = df[df['Substrate'] == substrate].copy()

        ax = axes[idx]
        scatter = ax.scatter(substrate_df['Km'], substrate_df['kcat'],
                           c=substrate_df['kcat_over_Km'], cmap='RdYlGn',
                           s=100, alpha=0.7, edgecolors='black')

        # Add enzyme labels
        for _, row in substrate_df.iterrows():
            ax.annotate(row['Enzyme'], (row['Km'], row['kcat']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=7, alpha=0.8)

        ax.set_xlabel('log10(Km) [M]', fontsize=10, fontweight='bold')
        ax.set_ylabel('log10(kcat) [s⁻¹]', fontsize=10, fontweight='bold')
        ax.set_title(substrate_names[substrate], fontsize=11, fontweight='bold')
        ax.grid(alpha=0.3)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('log10(kcat/Km)', fontsize=8)

    # Hide unused subplot
    axes[5].axis('off')

    fig.suptitle('kcat vs Km for Reported UDH Enzymes',
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_scatter.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_scatter.png")
    plt.close()

    # === 4. Summary statistics table ===
    print("\nCreating summary statistics...")
    summary_data = []

    for enzyme in enzymes:
        enzyme_df = df[df['Enzyme'] == enzyme]

        summary_data.append({
            'Enzyme': enzyme,
            'Mean_kcat_over_Km': enzyme_df['kcat_over_Km'].mean(),
            'Std_kcat_over_Km': enzyme_df['kcat_over_Km'].std(),
            'Mean_kcat': enzyme_df['kcat'].mean(),
            'Mean_Km': enzyme_df['Km'].mean(),
            'Best_Substrate': substrate_names[enzyme_df.loc[enzyme_df['kcat_over_Km'].idxmax(), 'Substrate']],
            'Best_kcat_over_Km': enzyme_df['kcat_over_Km'].max()
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('Mean_kcat_over_Km', ascending=False)

    # Save to CSV
    output_csv = f'{output_prefix}_summary.csv'
    summary_df.to_csv(output_csv, index=False)
    print(f"\nSaved: {output_csv}")

    # === 5. Overall ranking bar plot ===
    print("\nCreating overall ranking plot...")
    fig, ax = plt.subplots(figsize=(10, 8))

    bars = ax.barh(range(len(summary_df)), summary_df['Mean_kcat_over_Km'])

    # Color gradient
    colors = plt.cm.RdYlGn(np.linspace(0.3, 0.9, len(summary_df))[::-1])
    for bar, color in zip(bars, colors):
        bar.set_color(color)

    ax.set_yticks(range(len(summary_df)))
    ax.set_yticklabels(summary_df['Enzyme'], fontsize=10)
    ax.set_xlabel('Mean log10(kcat/Km) across all substrates', fontsize=11, fontweight='bold')
    ax.set_title('Overall Enzyme Rankings (Average across all substrates)',
                 fontsize=13, fontweight='bold', pad=20)
    ax.grid(axis='x', alpha=0.3)

    # Add rank and value labels
    for i, (idx_val, row) in enumerate(summary_df.iterrows()):
        ax.text(row['Mean_kcat_over_Km'] + 0.02, i,
               f"#{i+1} ({row['Mean_kcat_over_Km']:.2f})",
               va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_overall_ranking.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_overall_ranking.png")
    plt.close()

    # === 6. Print top rankings ===
    print("\n" + "="*80)
    print("TOP 5 ENZYMES (Average kcat/Km across all substrates)")
    print("="*80)
    for i, (_, row) in enumerate(summary_df.head().iterrows(), 1):
        print(f"\n{i}. {row['Enzyme']}")
        print(f"   Mean kcat/Km: {row['Mean_kcat_over_Km']:.3f} ± {row['Std_kcat_over_Km']:.3f}")
        print(f"   Best substrate: {row['Best_Substrate']} (kcat/Km = {row['Best_kcat_over_Km']:.3f})")

    # Print rankings by substrate
    print("\n" + "="*80)
    print("TOP 3 ENZYMES BY SUBSTRATE")
    print("="*80)
    for substrate in substrates:
        print(f"\n{substrate_names[substrate]}:")
        substrate_df = df[df['Substrate'] == substrate].sort_values('kcat_over_Km', ascending=False)
        for i, (_, row) in enumerate(substrate_df.head(3).iterrows(), 1):
            print(f"  {i}. {row['Enzyme']}: {row['kcat_over_Km']:.3f}")

    print("\n" + "="*80)
    print("Visualization complete!")
    print("="*80)

if __name__ == "__main__":
    predictions_file = 'reported_udh_catapro_predictions.csv'
    create_rankings_visualization(predictions_file)
