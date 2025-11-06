"""
Visualize kcat/Km values by substrate across all scenarios
- One file per substrate
- 8 scenarios in 2 rows × 4 columns layout
- Violin plots for distribution with enhanced modern color scheme
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
plt.rcParams['font.size'] = 14  # Increased from 10
plt.rcParams['figure.dpi'] = 100
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 14  # Axis labels
plt.rcParams['axes.titlesize'] = 16  # Subplot titles
plt.rcParams['xtick.labelsize'] = 12  # X-axis tick labels
plt.rcParams['ytick.labelsize'] = 12  # Y-axis tick labels
plt.rcParams['legend.fontsize'] = 14  # Legend text
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.15
plt.rcParams['axes.facecolor'] = '#F8F9FA'
plt.rcParams['figure.facecolor'] = 'white'

# Define modern color palette
# Using a professional gradient color scheme
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

def load_data_for_scenario(scenario_num):
    """Load and merge data for a specific scenario"""

    # Load catapro predictions
    df_pred = pd.read_csv('catapro_predict_score.csv')
    df_pred['substrate'] = df_pred['smiles'].map(SUBSTRATE_MAPPING)

    # Load cluster assignments for this scenario
    sheet_name = f'Full_S{scenario_num}'
    df_clusters = pd.read_excel('reported_udh_analysis_20251106_010710.xlsx',
                                 sheet_name=sheet_name)
    df_clusters = df_clusters[['fasta_id', 'Cluster_ID']].copy()
    df_clusters.columns = ['fasta_id', 'cluster']

    # Load reported UDH information for this scenario
    df_reported = pd.read_excel('reported_udh_analysis_20251106_010710.xlsx',
                                sheet_name='Reported_UDH_Distribution')
    df_reported = df_reported[df_reported['Scenario'] == scenario_num].copy()
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

    # Separate noise and regular clusters
    noise_mask = cluster_stats['cluster'] == -1
    noise_stats = cluster_stats[noise_mask].copy()
    regular_stats = cluster_stats[~noise_mask].copy()

    # Sort regular clusters by median (descending)
    regular_stats = regular_stats.sort_values('median', ascending=False)

    # Combine: regular clusters first, then noise at the end
    sorted_stats = pd.concat([regular_stats, noise_stats], ignore_index=True)

    return sorted_stats, df_sub

def create_gradient_violin(ax, data, position, width, colors):
    """Create a gradient-filled violin plot"""
    from matplotlib.collections import PolyCollection

    # Calculate violin shape using kernel density estimation
    from scipy import stats

    if len(data) > 1:
        kde = stats.gaussian_kde(data, bw_method='scott')
        y_range = np.linspace(data.min(), data.max(), 100)
        density = kde(y_range)
        density = density / density.max() * width / 2
    else:
        y_range = [data[0]]
        density = [width / 2]

    # Create vertices for the violin shape
    vertices_left = [(position - d, y) for d, y in zip(density, y_range)]
    vertices_right = [(position + d, y) for d, y in zip(density[::-1], y_range[::-1])]
    vertices = vertices_left + vertices_right

    # Create gradient effect
    gradient = np.linspace(0, 1, len(vertices))

    # Create polygon with gradient
    poly = plt.Polygon(vertices, facecolor=colors['face'],
                       edgecolor=colors['edge'], linewidth=2,
                       alpha=0.7, zorder=5)
    ax.add_patch(poly)

    return poly

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

    # Create labels with improved styling
    cluster_labels = []
    for cluster in cluster_order:
        stats = cluster_stats[cluster_stats['cluster'] == cluster].iloc[0]
        count_val = int(stats['count'])
        if cluster == -1:
            cluster_labels.append(f'Noise\n({count_val})')
        else:
            cluster_labels.append(f'C{int(cluster)}\n({count_val})')

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

def create_scenario_visualization(scenario_num, output_dir):
    """Create visualization for one scenario with all 5 substrates"""

    print(f"\nProcessing Scenario {scenario_num}: {SCENARIO_NAMES[scenario_num]}")

    # Load data for this scenario
    df, df_reported = load_data_for_scenario(scenario_num)

    # Create figure with enhanced layout
    fig = plt.figure(figsize=(28, 34))
    gs = GridSpec(5, 1, figure=fig, hspace=0.45, left=0.1, right=0.94, top=0.92, bottom=0.05)

    axes = [fig.add_subplot(gs[i, 0]) for i in range(5)]

    # Plot each substrate
    for idx, substrate in enumerate(SUBSTRATE_ORDER):
        ax = axes[idx]
        cluster_stats, df_sub = prepare_cluster_statistics(df, substrate)

        print(f"  {substrate}: {len(cluster_stats)} clusters")

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
    fig.suptitle(f'Scenario {scenario_num}: {SCENARIO_NAMES[scenario_num]}\n' +
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
               label='Cluster Distribution', alpha=0.7),
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
    output_path = output_dir / f'scenario_{scenario_num}_by_substrate_violin.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {output_path}")

    plt.close(fig)

def main():
    """Main execution"""

    print("="*80)
    print("ENHANCED SCENARIO-WISE VISUALIZATION WITH MODERN VIOLIN PLOTS")
    print("="*80)

    # Create output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = Path(f'scenario_violin_{timestamp}')
    output_dir.mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir}")

    # Process each scenario
    for scenario_num in range(1, 9):
        try:
            create_scenario_visualization(scenario_num, output_dir)
        except Exception as e:
            print(f"Error processing Scenario {scenario_num}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{'='*80}")
    print(f"All enhanced visualizations saved to: {output_dir}")
    print(f"{'='*80}\n")

if __name__ == "__main__":
    main()