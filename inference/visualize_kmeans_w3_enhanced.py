"""
Visualize k-means clustering results with window=3 for kcat/Km values by substrate
- Enhanced modern color scheme with larger fonts for better readability
- Violin plots for distribution
- No noise clusters (k-means specific)
- Highlight 14 reported UDH data points
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns
from matplotlib.gridspec import GridSpec
from datetime import datetime
from matplotlib.patches import PathPatch
import matplotlib.colors as mcolors

# Set style for modern, clean appearance with larger fonts
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
plt.rcParams['font.size'] = 14  # Larger base font
plt.rcParams['figure.dpi'] = 100
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.15
plt.rcParams['axes.facecolor'] = '#F8F9FA'
plt.rcParams['figure.facecolor'] = 'white'

# Define modern color palette
VIOLIN_COLORS = {
    'NAD+': {'face': '#FF6B6B', 'edge': '#C92A2A', 'median': '#862E2E'},
    'α-D-Glucuronate': {'face': '#4ECDC4', 'edge': '#15AAA0', 'median': '#0B7C75'},
    'α-L-Guluronate': {'face': '#45B7D1', 'edge': '#2E86AB', 'median': '#1E5F8E'},
    'β-D-Galacturonate': {'face': '#96CEB4', 'edge': '#69B098', 'median': '#4A8B70'},
    'β-D-Mannuronate': {'face': '#DDA0DD', 'edge': '#BA55D3', 'median': '#8B008B'}
}

# Highlight colors for reported UDHs
HIGHLIGHT_COLOR = '#FFD700'  # Gold
HIGHLIGHT_EDGE = '#FF8C00'   # Dark orange
ANNOTATION_COLOR = '#2C3E50'  # Dark blue-gray

# Define substrate mapping
SUBSTRATE_MAPPING = {
    '[C@@H]1([C@@H]([C@@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O': 'α-D-Glucuronate',
    '[C@@H]1([C@H]([C@H](O[C@H]([C@@H]1O)O)C(=O)[O-])O)O': 'α-L-Guluronate',
    '[C@@H]1([C@@H]([C@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O': 'β-D-Galacturonate',
    '[C@@H]1([C@@H]([C@H](O[C@@H]([C@@H]1O)O)C(=O)[O-])O)O': 'β-D-Mannuronate',
    'C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N': 'NAD+'
}

SCENARIO_NAMES = {
    1: 'All binding sites',
    2: 'NAD binding sites',
    3: 'Substrate binding sites',
    4: 'Direct binding sites',
    5: 'Proximal binding sites',
    6: 'No binding site info',
    7: 'Substrate direct only',
    8: 'Substrate proximal only'
}

SUBSTRATE_ORDER = ['NAD+', 'α-D-Glucuronate', 'α-L-Guluronate',
                   'β-D-Galacturonate', 'β-D-Mannuronate']

# Target window size for analysis
TARGET_WINDOW_SIZE = 3

def load_data_for_scenario(scenario_num, excel_file, window_size=TARGET_WINDOW_SIZE):
    """Load and merge data for a specific scenario and window size"""

    # Load catapro predictions
    df_pred = pd.read_csv('catapro_predict_score.csv')
    df_pred['substrate'] = df_pred['smiles'].map(SUBSTRATE_MAPPING)

    # Load cluster assignments for this scenario and window
    # New file format uses S{scenario}_W{window} naming
    sheet_name = f'S{scenario_num}_W{window_size}'

    try:
        df_clusters = pd.read_excel(excel_file, sheet_name=sheet_name)
    except ValueError:
        print(f"  Warning: Sheet '{sheet_name}' not found")
        return pd.DataFrame(), pd.DataFrame()

    print(f"  Using window size = {window_size}, {len(df_clusters)} sequences")

    df_clusters = df_clusters[['fasta_id', 'Cluster_ID']].copy()
    df_clusters.columns = ['fasta_id', 'cluster']

    # Load reported UDH information for this scenario
    df_reported = pd.read_excel(excel_file, sheet_name='Reported_UDH_Distribution')

    # Filter by scenario and window size
    df_reported = df_reported[
        (df_reported['Scenario'] == scenario_num) &
        (df_reported['Window_Size'] == window_size)
    ].copy()

    df_reported = df_reported[['UDH_Name', 'FASTA_ID', 'Cluster']].copy()
    df_reported.columns = ['udh_name', 'fasta_id', 'cluster']

    # Merge predictions with cluster assignments
    df_merged = df_pred.merge(df_clusters, on='fasta_id', how='left')

    # Mark reported UDHs
    df_merged['is_reported'] = df_merged['fasta_id'].isin(df_reported['fasta_id'])

    # Add UDH names for reported enzymes
    df_merged = df_merged.merge(
        df_reported[['fasta_id', 'udh_name']],
        on='fasta_id',
        how='left'
    )

    return df_merged, df_reported

def prepare_cluster_statistics(df, substrate):
    """Calculate statistics for each cluster for a specific substrate"""

    # Filter for specific substrate
    df_sub = df[df['substrate'] == substrate].copy()

    # Group by cluster
    cluster_stats = df_sub.groupby('cluster').agg({
        'pred_log10[kcat/Km(s^-1mM^-1)]': ['mean', 'median', 'std', 'count']
    }).reset_index()

    cluster_stats.columns = ['cluster', 'mean', 'median', 'std', 'count']

    # For k-means, all clusters are regular (no noise cluster -1)
    # Sort clusters by median (descending) for better visualization
    cluster_stats = cluster_stats.sort_values('median', ascending=False)

    return cluster_stats, df_sub

def plot_substrate_cluster_violin(ax, df, df_sub, cluster_stats, scenario_num, substrate, is_top_row=False):
    """Plot cluster distribution with enhanced violin plot for a specific substrate/scenario"""

    # Get color scheme for this substrate
    colors = VIOLIN_COLORS.get(substrate, {'face': '#95A5A6', 'edge': '#7F8C8D', 'median': '#34495E'})

    # Get sorted cluster order
    cluster_order = cluster_stats['cluster'].tolist()
    cluster_positions = {cluster: i for i, cluster in enumerate(cluster_order)}

    # Prepare data for violin plot
    violin_data = []
    positions = []

    for cluster in cluster_order:
        cluster_data = df_sub[df_sub['cluster'] == cluster]['pred_log10[kcat/Km(s^-1mM^-1)]'].values
        if len(cluster_data) > 0:
            violin_data.append(cluster_data)
            positions.append(cluster_positions[cluster])

    # Draw violin plot with enhanced styling
    if len(violin_data) > 0:
        vp = ax.violinplot(violin_data, positions=positions, widths=0.7,
                          showmeans=False, showmedians=True, showextrema=True)

        # Style violin bodies with gradient effect
        for i, pc in enumerate(vp['bodies']):
            # Create gradient effect by setting alpha and colors
            pc.set_facecolor(colors['face'])
            pc.set_alpha(0.65)
            pc.set_edgecolor(colors['edge'])
            pc.set_linewidth(2)

            # Add subtle gradient effect using path effects
            from matplotlib import patheffects
            pc.set_path_effects([patheffects.SimplePatchShadow(offset=(1, -1),
                                                               shadow_rgbFace=colors['edge'],
                                                               alpha=0.3)])

        # Style median lines with bold color
        vp['cmedians'].set_color(colors['median'])
        vp['cmedians'].set_linewidth(3)
        vp['cmedians'].set_alpha(0.9)

        # Style extrema with matching colors
        vp['cmaxes'].set_color(colors['edge'])
        vp['cmaxes'].set_linewidth(2)
        vp['cmaxes'].set_alpha(0.8)
        vp['cmins'].set_color(colors['edge'])
        vp['cmins'].set_linewidth(2)
        vp['cmins'].set_alpha(0.8)
        vp['cbars'].set_color(colors['edge'])
        vp['cbars'].set_linewidth(1.5)
        vp['cbars'].set_alpha(0.6)

    # Plot reported UDH points with enhanced styling
    for cluster in cluster_order:
        cluster_data = df_sub[df_sub['cluster'] == cluster]
        reported_points = cluster_data[cluster_data['is_reported'] == True]

        if len(reported_points) > 0:
            x_pos = cluster_positions[cluster]

            for _, row in reported_points.iterrows():
                # Add jitter for visibility
                jitter = np.random.normal(0, 0.08)

                # Create glow effect for reported points
                for size, alpha in [(120, 0.2), (100, 0.3), (80, 1.0)]:
                    ax.scatter(x_pos + jitter,
                              row['pred_log10[kcat/Km(s^-1mM^-1)]'],
                              alpha=alpha, s=size, color=HIGHLIGHT_COLOR,
                              edgecolors='none' if size > 80 else HIGHLIGHT_EDGE,
                              linewidth=0 if size > 80 else 2.5,
                              marker='*',
                              zorder=10 + (80 - size) // 20)

                # Add UDH name with improved styling
                ax.annotate(row['udh_name'],
                           xy=(x_pos + jitter, row['pred_log10[kcat/Km(s^-1mM^-1)]']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=10, color=ANNOTATION_COLOR, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3',
                                   facecolor='white', alpha=0.95,
                                   edgecolor=HIGHLIGHT_EDGE, linewidth=2))

    # Create labels with improved styling (k-means specific - no noise cluster)
    cluster_labels = []
    for cluster in cluster_order:
        stats = cluster_stats[cluster_stats['cluster'] == cluster].iloc[0]
        count_val = int(stats['count'])
        cluster_labels.append(f'K{int(cluster)}\n({count_val})')

    # Customize plot with modern styling
    ax.set_xticks(range(len(cluster_order)))
    ax.set_xticklabels(cluster_labels, fontsize=11, rotation=45, ha='right')

    if is_top_row:
        ax.set_title(f'Scenario {scenario_num}: {SCENARIO_NAMES[scenario_num]}',
                    fontsize=16, fontweight='bold', pad=10,
                    bbox=dict(boxstyle='round,pad=0.6',
                             facecolor='white', alpha=0.8,
                             edgecolor='gray', linewidth=1.5))

    ax.set_ylabel('log₁₀(kcat/Km) [s⁻¹mM⁻¹]', fontsize=14, fontweight='bold')

    # Enhanced grid styling
    ax.grid(True, alpha=0.2, linestyle='--', axis='y', color='gray')
    ax.grid(True, alpha=0.1, linestyle=':', axis='x', color='gray')
    ax.set_axisbelow(True)

    # Add subtle background gradient
    ax.axhspan(-0.6, 2.1, alpha=0.02, color='blue', zorder=0)

    # Set consistent y-axis range with minor ticks
    ax.set_ylim(-0.6, 2.1)
    ax.set_yticks(np.arange(-0.5, 2.1, 0.5))
    ax.set_yticks(np.arange(-0.6, 2.1, 0.1), minor=True)
    ax.tick_params(axis='y', which='minor', length=2, color='gray')

def create_scenario_visualization(scenario_num, excel_file, output_dir):
    """Create visualization for one scenario with all 5 substrates"""

    print(f"\nProcessing Scenario {scenario_num}: {SCENARIO_NAMES[scenario_num]}")

    try:
        # Load data for this scenario
        df, df_reported = load_data_for_scenario(scenario_num, excel_file, TARGET_WINDOW_SIZE)

        # Check if we have data
        if df.empty or df['cluster'].isna().all():
            print(f"  Warning: No data found for window size {TARGET_WINDOW_SIZE}")
            return False

        # Create figure with enhanced layout
        fig = plt.figure(figsize=(28, 34))
        gs = GridSpec(5, 1, figure=fig, hspace=0.45, left=0.1, right=0.94, top=0.92, bottom=0.05)

        axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]

        # Plot each substrate
        for idx, substrate in enumerate(SUBSTRATE_ORDER):
            ax = axes[idx]
            cluster_stats, df_sub = prepare_cluster_statistics(df, substrate)

            print(f"  {substrate}: {len(cluster_stats)} k-means clusters")

            is_top_row = (idx == 0)
            plot_substrate_cluster_violin(ax, df, df_sub, cluster_stats,
                                         scenario_num, substrate, is_top_row)

            # Add substrate label with enhanced styling
            substrate_color = VIOLIN_COLORS.get(substrate, {'face': '#95A5A6'})['edge']
            ax.text(-0.09, 0.5, substrate, transform=ax.transAxes,
                   fontsize=18, fontweight='bold', rotation=90,
                   verticalalignment='center', horizontalalignment='right',
                   color=substrate_color,
                   bbox=dict(boxstyle='round,pad=0.4',
                            facecolor='white', alpha=0.9,
                            edgecolor=substrate_color, linewidth=2.5))

        # Add main title with enhanced styling
        fig.suptitle(f'K-means Clustering (Window={TARGET_WINDOW_SIZE}) - Scenario {scenario_num}: {SCENARIO_NAMES[scenario_num]}\n' +
                    f'Cluster-wise Catalytic Efficiency Distribution by Substrate',
                    fontsize=24, fontweight='bold', y=0.97,
                    bbox=dict(boxstyle='round,pad=0.8',
                             facecolor='#F0F8FF', alpha=0.9,
                             edgecolor='#4169E1', linewidth=3))

        # Add enhanced legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor='#4ECDC4', markersize=14,
                   markeredgecolor='#15AAA0', markeredgewidth=2,
                   label='K-means Cluster Distribution', alpha=0.7),
            Line2D([0], [0], marker='_', color='#0B7C75',
                   markersize=12, linewidth=4,
                   label='Median Value', markeredgewidth=0),
            Line2D([0], [0], marker='*', color='w', markerfacecolor=HIGHLIGHT_COLOR,
                   markersize=14, label='14 Reported UDHs',
                   markeredgecolor=HIGHLIGHT_EDGE, markeredgewidth=2.5)
        ]

        legend = fig.legend(handles=legend_elements, loc='upper right',
                           bbox_to_anchor=(0.98, 0.96), fontsize=15,
                           framealpha=0.95, fancybox=True, shadow=True,
                           title='Legend', title_fontsize=16,
                           edgecolor='gray', facecolor='white')
        legend.get_frame().set_linewidth(2.5)

        # Save with high quality
        output_path = output_dir / f'kmeans_w{TARGET_WINDOW_SIZE}_scenario_{scenario_num}_by_substrate_violin.png'
        fig.savefig(output_path, dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"  Saved: {output_path}")

        plt.close(fig)
        return True

    except Exception as e:
        print(f"  Error processing scenario {scenario_num}: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main execution"""

    print("="*80)
    print(f"K-MEANS CLUSTERING VISUALIZATION (Window Size = {TARGET_WINDOW_SIZE})")
    print("Enhanced Violin Plots with Modern Color Scheme")
    print("="*80)

    # Define the Excel file to use
    excel_file = 'kmeans_14udh_analysis_20251106_160853.xlsx'

    # Check if file exists
    if not Path(excel_file).exists():
        print(f"\nError: File '{excel_file}' not found.")
        return

    # Create output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = Path(f'kmeans_window{TARGET_WINDOW_SIZE}_violin_{timestamp}')
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir}")

    # Process each scenario (skip scenario 6 if missing)
    success_count = 0
    for scenario_num in [1, 2, 3, 4, 5, 7, 8]:  # Note: Scenario 6 might be special
        try:
            if create_scenario_visualization(scenario_num, excel_file, output_dir):
                success_count += 1
        except Exception as e:
            print(f"Error processing Scenario {scenario_num}: {e}")
            continue

    if success_count > 0:
        print(f"\n{'='*80}")
        print(f"Successfully created {success_count} visualizations")
        print(f"All visualizations saved to: {output_dir}")
        print(f"{'='*80}\n")
    else:
        print(f"\n{'='*80}")
        print(f"No visualizations were created.")
        print(f"Please check the data availability.")
        print(f"{'='*80}\n")

if __name__ == "__main__":
    main()