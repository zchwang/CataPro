"""
Sequence-to-Sequence Comparison for Binding Site Residues
Compares ~1000 sequences to identify which amino acids at binding sites
increase/decrease kinetic parameters for NAD and 4 uronates
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')


def calculate_sequence_identity(seq1, seq2):
    """
    Calculate pairwise sequence identity between two aligned sequences

    Parameters:
    -----------
    seq1, seq2 : str
        Aligned sequences (may contain gaps '-')

    Returns:
    --------
    float : Identity percentage (0-100)
    """
    if len(seq1) != len(seq2):
        return 0.0

    matches = 0
    valid_positions = 0

    for a, b in zip(seq1, seq2):
        # Skip positions where either sequence has a gap
        if a != '-' and b != '-':
            valid_positions += 1
            if a == b:
                matches += 1

    if valid_positions == 0:
        return 0.0

    return (matches / valid_positions) * 100


def calculate_average_pairwise_identity(sequences, max_pairs=1000):
    """
    Calculate average pairwise sequence identity for a group of sequences

    Parameters:
    -----------
    sequences : list of str
        List of aligned sequences
    max_pairs : int
        Maximum number of pairs to sample (default: 1000)
        If total pairs > max_pairs, randomly sample

    Returns:
    --------
    dict with:
        - mean_identity: Average pairwise identity
        - min_identity: Minimum pairwise identity
        - max_identity: Maximum pairwise identity
        - n_pairs: Number of pairs compared
    """
    if len(sequences) < 2:
        return {
            'mean_identity': 100.0,
            'min_identity': 100.0,
            'max_identity': 100.0,
            'n_pairs': 0
        }

    identities = []

    # Get all possible pairs
    all_pairs = list(combinations(sequences, 2))
    n_total_pairs = len(all_pairs)

    # If too many pairs, randomly sample
    if n_total_pairs > max_pairs:
        sampled_pairs = np.random.choice(len(all_pairs), size=max_pairs, replace=False)
        pairs_to_compare = [all_pairs[i] for i in sampled_pairs]
    else:
        pairs_to_compare = all_pairs

    # Calculate identities
    for seq1, seq2 in pairs_to_compare:
        identity = calculate_sequence_identity(seq1, seq2)
        identities.append(identity)

    return {
        'mean_identity': np.mean(identities),
        'min_identity': np.min(identities),
        'max_identity': np.max(identities),
        'n_pairs': len(identities),
        'sampled': n_total_pairs > max_pairs
    }


class SequenceComparisonAnalyzer:
    """
    Analyzes which amino acids at specific positions lead to higher/lower activity
    """

    def __init__(self, binding_site_positions, window_size=3):
        """
        Parameters:
        -----------
        binding_site_positions : list or dict
            Binding site positions (0-indexed)
        window_size : int
            Window around each position (1-3)
        """
        if isinstance(binding_site_positions, dict):
            self.binding_sites = binding_site_positions
        else:
            self.binding_sites = {pos: f'site_{i+1}' for i, pos in enumerate(binding_site_positions)}

        self.window_size = window_size
        self.positions_with_window = self._expand_positions()
        self.substrates = []
        self.results = {}

    def _expand_positions(self):
        """Expand binding site positions to include window"""
        position_info = {}

        for core_pos, site_name in self.binding_sites.items():
            # Core position
            position_info[core_pos] = {
                'site': site_name,
                'window': 0,
                'core_position': core_pos,
                'label': f'{site_name}_core'
            }

            # Window positions
            for window in range(1, self.window_size + 1):
                upstream_pos = core_pos - window
                position_info[upstream_pos] = {
                    'site': site_name,
                    'window': -window,
                    'core_position': core_pos,
                    'label': f'{site_name}_w-{window}'
                }

                downstream_pos = core_pos + window
                position_info[downstream_pos] = {
                    'site': site_name,
                    'window': window,
                    'core_position': core_pos,
                    'label': f'{site_name}_w+{window}'
                }

        return position_info

    def load_substrate_data(self, substrate_files):
        """
        Load CataPro prediction files for multiple substrates

        Parameters:
        -----------
        substrate_files : dict
            {substrate_name: file_path}
            Example: {
                'NADplus': 'catapro_predict_score_NADplus.csv',
                'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
                ...
            }
        """
        self.substrate_data = {}
        self.substrates = list(substrate_files.keys())

        print("Loading substrate data...")
        for substrate, filepath in substrate_files.items():
            if not os.path.exists(filepath):
                print(f"Warning: {filepath} not found, skipping {substrate}")
                continue

            df = pd.read_csv(filepath)

            # Clean column names
            if 'Unnamed: 0' in df.columns:
                df = df.drop('Unnamed: 0', axis=1)

            self.substrate_data[substrate] = df
            print(f"  {substrate}: {len(df)} sequences")

        print(f"\nLoaded {len(self.substrate_data)} substrates")
        return self.substrate_data

    def extract_position_data(self, aligned_seq_col='aligned_seq', min_coverage=0.2):
        """
        Extract amino acid information at binding site positions for all sequences

        Parameters:
        -----------
        aligned_seq_col : str
            Column name containing aligned sequences
        min_coverage : float
            Minimum fraction of sequences that must have non-gap at this position
            (default: 0.2 = at least 20% of sequences must have amino acid)

        Returns a combined dataframe with:
        - fasta_id
        - substrate
        - amino acids at each position
        - kinetic parameters
        """
        # First pass: count gaps at each position to filter out insertion-heavy positions
        print(f"\nFiltering positions by coverage (min_coverage={min_coverage})...")

        first_substrate = list(self.substrate_data.values())[0]
        n_sequences = len(first_substrate)

        positions_to_keep = {}
        positions_removed = []

        for pos, info in sorted(self.positions_with_window.items()):
            gap_count = 0
            total_count = 0

            for df in self.substrate_data.values():
                for idx, row in df.iterrows():
                    seq = row[aligned_seq_col]
                    if pos >= 0 and pos < len(seq):
                        total_count += 1
                        if seq[pos] == '-':
                            gap_count += 1

            if total_count > 0:
                coverage = 1 - (gap_count / total_count)
                if coverage >= min_coverage:
                    positions_to_keep[pos] = info
                else:
                    positions_removed.append({
                        'position': pos + 1,
                        'label': info['label'],
                        'coverage': coverage,
                        'gap_count': gap_count,
                        'total_count': total_count
                    })

        if positions_removed:
            print(f"\nRemoved {len(positions_removed)} positions with low coverage:")
            for p in positions_removed:
                print(f"  Pos {p['position']} ({p['label']}): {p['coverage']*100:.1f}% coverage "
                      f"({p['gap_count']}/{p['total_count']} gaps)")

        print(f"\nKeeping {len(positions_to_keep)} positions with coverage >= {min_coverage*100:.0f}%")

        # Update positions_with_window to only include filtered positions
        self.positions_with_window = positions_to_keep

        # Second pass: extract data
        all_data = []

        for substrate, df in self.substrate_data.items():
            if aligned_seq_col not in df.columns:
                print(f"Warning: '{aligned_seq_col}' not found in {substrate}, skipping")
                continue

            for idx, row in df.iterrows():
                seq = row[aligned_seq_col]

                seq_data = {
                    'fasta_id': row.get('fasta_id', f'seq_{idx}'),
                    'substrate': substrate,
                }

                # Extract amino acids at binding site positions (only filtered ones)
                for pos, info in sorted(self.positions_with_window.items()):
                    if pos >= 0 and pos < len(seq):
                        aa = seq[pos]
                        if aa != '-':
                            seq_data[f'pos_{pos+1}'] = aa
                            seq_data[f'pos_{pos+1}_site'] = info['site']
                            seq_data[f'pos_{pos+1}_window'] = info['window']
                        else:
                            seq_data[f'pos_{pos+1}'] = None

                # Add kinetic parameters
                if 'pred_log10[kcat(s^-1)]' in df.columns:
                    seq_data['log_kcat'] = row['pred_log10[kcat(s^-1)]']
                if 'pred_log10[Km(mM)]' in df.columns:
                    seq_data['log_Km'] = row['pred_log10[Km(mM)]']
                if 'pred_log10[kcat/Km(s^-1mM^-1)]' in df.columns:
                    seq_data['log_kcat_Km'] = row['pred_log10[kcat/Km(s^-1mM^-1)]']

                all_data.append(seq_data)

        self.combined_data = pd.DataFrame(all_data)

        # Also store position_info for later visualization
        self.position_info = {}
        for pos, info in sorted(self.positions_with_window.items()):
            self.position_info[pos] = info

        print(f"\nCombined data: {len(self.combined_data)} sequences × {len(self.substrates)} substrates")
        return self.combined_data

    def analyze_amino_acid_effects(self, kinetic_params=['log_kcat', 'log_Km', 'log_kcat_Km'],
                                   min_sequences=5):
        """
        For each position and each substrate, compare sequences with different amino acids

        Parameters:
        -----------
        kinetic_params : list
            Which kinetic parameters to analyze
        min_sequences : int
            Minimum number of sequences required for each amino acid
        """
        print("\nAnalyzing amino acid effects at each position...")

        results = []
        total_positions = len(self.positions_with_window)

        for i, (pos, info) in enumerate(sorted(self.positions_with_window.items())):
            if (i + 1) % 10 == 0:
                print(f"  Processing position {i+1}/{total_positions}")

            pos_col = f'pos_{pos+1}'

            if pos_col not in self.combined_data.columns:
                continue

            # Analyze for each substrate
            for substrate in self.substrates:
                substrate_data = self.combined_data[self.combined_data['substrate'] == substrate].copy()

                if len(substrate_data) == 0:
                    continue

                # Get amino acid distribution at this position
                aa_counts = substrate_data[pos_col].value_counts()

                # Filter amino acids: exclude n=1 (unreliable) and keep only n >= min_sequences
                valid_aas = aa_counts[(aa_counts > 1) & (aa_counts >= min_sequences)].index.tolist()

                if len(valid_aas) < 2:
                    continue

                # Analyze each kinetic parameter
                for param in kinetic_params:
                    if param not in substrate_data.columns:
                        continue

                    # Get mean values for each amino acid
                    aa_means = {}
                    aa_stds = {}
                    aa_counts_dict = {}

                    for aa in valid_aas:
                        aa_data = substrate_data[substrate_data[pos_col] == aa][param].dropna()
                        if len(aa_data) >= min_sequences:
                            aa_means[aa] = aa_data.mean()
                            aa_stds[aa] = aa_data.std()
                            aa_counts_dict[aa] = len(aa_data)

                    if len(aa_means) < 2:
                        continue

                    # Find best and worst amino acids
                    best_aa = max(aa_means, key=aa_means.get)
                    worst_aa = min(aa_means, key=aa_means.get)

                    # Statistical test between best and worst
                    best_data = substrate_data[substrate_data[pos_col] == best_aa][param].dropna()
                    worst_data = substrate_data[substrate_data[pos_col] == worst_aa][param].dropna()

                    if len(best_data) >= min_sequences and len(worst_data) >= min_sequences:
                        stat, pval = stats.mannwhitneyu(best_data, worst_data, alternative='two-sided')
                        effect_size = (best_data.mean() - worst_data.mean()) / substrate_data[param].std()
                    else:
                        pval = np.nan
                        effect_size = np.nan

                    # Calculate sequence identity for best and worst AA groups
                    # Need to get aligned sequences from original substrate_data
                    best_fasta_ids = substrate_data[substrate_data[pos_col] == best_aa]['fasta_id'].tolist()
                    worst_fasta_ids = substrate_data[substrate_data[pos_col] == worst_aa]['fasta_id'].tolist()

                    # Get aligned sequences from original data
                    substrate_original = self.substrate_data[substrate]
                    best_sequences = substrate_original[substrate_original['fasta_id'].isin(best_fasta_ids)]['aligned_seq'].tolist()
                    worst_sequences = substrate_original[substrate_original['fasta_id'].isin(worst_fasta_ids)]['aligned_seq'].tolist()

                    best_identity = calculate_average_pairwise_identity(best_sequences)
                    worst_identity = calculate_average_pairwise_identity(worst_sequences)

                    result = {
                        'position': pos + 1,
                        'site': info['site'],
                        'window': info['window'],
                        'label': info['label'],
                        'substrate': substrate,
                        'parameter': param,
                        'n_amino_acids': len(aa_means),
                        'best_aa': best_aa,
                        'best_mean': aa_means[best_aa],
                        'best_std': aa_stds[best_aa],
                        'best_count': aa_counts_dict[best_aa],
                        'best_seq_identity_mean': best_identity['mean_identity'],
                        'best_seq_identity_min': best_identity['min_identity'],
                        'best_seq_identity_max': best_identity['max_identity'],
                        'worst_aa': worst_aa,
                        'worst_mean': aa_means[worst_aa],
                        'worst_std': aa_stds[worst_aa],
                        'worst_count': aa_counts_dict[worst_aa],
                        'worst_seq_identity_mean': worst_identity['mean_identity'],
                        'worst_seq_identity_min': worst_identity['min_identity'],
                        'worst_seq_identity_max': worst_identity['max_identity'],
                        'delta': aa_means[best_aa] - aa_means[worst_aa],
                        'effect_size': effect_size,
                        'p_value': pval,
                        'significant': pval < 0.05 if not np.isnan(pval) else False,
                        'all_aa_means': str(aa_means),
                        'all_aa_counts': str(aa_counts_dict)
                    }

                    results.append(result)

        self.aa_effects = pd.DataFrame(results)
        print(f"\nAnalysis complete: {len(self.aa_effects)} position-substrate-parameter combinations")
        return self.aa_effects

    def compare_substrates(self, position, parameter='log_kcat_Km'):
        """
        Compare amino acid preferences across substrates for a specific position

        Returns which amino acids are preferred for each substrate
        """
        pos_data = self.aa_effects[
            (self.aa_effects['position'] == position) &
            (self.aa_effects['parameter'] == parameter)
        ]

        if len(pos_data) == 0:
            print(f"No data for position {position}, {parameter}")
            return None

        comparison = pos_data[['substrate', 'best_aa', 'best_mean', 'worst_aa', 'worst_mean',
                              'delta', 'p_value', 'significant']].copy()

        return comparison

    def find_substrate_specific_residues(self, parameter='log_kcat_Km', min_delta=0.1):
        """
        Find positions where different substrates prefer different amino acids
        (substrate-specific determinants)
        """
        print(f"\nFinding substrate-specific residues for {parameter}...")

        substrate_specific = []

        positions = self.aa_effects['position'].unique()

        for pos in positions:
            pos_data = self.aa_effects[
                (self.aa_effects['position'] == pos) &
                (self.aa_effects['parameter'] == parameter)
            ]

            if len(pos_data) < 2:
                continue

            # Check if different substrates prefer different amino acids
            best_aas = pos_data['best_aa'].unique()

            if len(best_aas) > 1 and pos_data['delta'].max() > min_delta:
                # This position shows substrate specificity
                info = pos_data.iloc[0]

                substrate_preferences = {}
                for _, row in pos_data.iterrows():
                    substrate_preferences[row['substrate']] = {
                        'best_aa': row['best_aa'],
                        'mean': row['best_mean'],
                        'delta': row['delta']
                    }

                substrate_specific.append({
                    'position': pos,
                    'site': info['site'],
                    'window': info['window'],
                    'label': info['label'],
                    'parameter': parameter,
                    'n_substrates': len(pos_data),
                    'n_different_preferences': len(best_aas),
                    'max_delta': pos_data['delta'].max(),
                    'substrate_preferences': str(substrate_preferences)
                })

        self.substrate_specific = pd.DataFrame(substrate_specific)
        print(f"Found {len(self.substrate_specific)} substrate-specific positions")
        return self.substrate_specific

    def generate_report(self, output_dir='./'):
        """Generate comprehensive report"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = os.path.join(output_dir, f"sequence_comparison_report_{timestamp}.txt")

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("SEQUENCE COMPARISON ANALYSIS REPORT\n")
            f.write("="*80 + "\n\n")

            # Configuration
            f.write("ANALYSIS CONFIGURATION\n")
            f.write("-"*80 + "\n")
            f.write(f"Core binding sites: {len(self.binding_sites)}\n")
            f.write(f"Window size: ±{self.window_size}\n")
            f.write(f"Total positions analyzed: {len(self.positions_with_window)}\n")
            f.write(f"Substrates: {', '.join(self.substrates)}\n\n")

            # Data summary
            if hasattr(self, 'combined_data'):
                f.write("DATA SUMMARY\n")
                f.write("-"*80 + "\n")
                for substrate in self.substrates:
                    n_seqs = len(self.combined_data[self.combined_data['substrate'] == substrate])
                    f.write(f"  {substrate}: {n_seqs} sequences\n")
                f.write(f"\nTotal: {len(self.combined_data)} sequence-substrate combinations\n\n")

            # Amino acid effects
            if hasattr(self, 'aa_effects'):
                f.write("AMINO ACID EFFECTS\n")
                f.write("="*80 + "\n\n")

                for param in self.aa_effects['parameter'].unique():
                    f.write(f"\n{param}\n")
                    f.write("-"*80 + "\n\n")

                    param_data = self.aa_effects[self.aa_effects['parameter'] == param]

                    for substrate in self.substrates:
                        sub_data = param_data[param_data['substrate'] == substrate]
                        sig_data = sub_data[sub_data['significant']]

                        f.write(f"\n{substrate}:\n")
                        f.write(f"  Total positions: {len(sub_data)}\n")
                        f.write(f"  Significant effects (p<0.05): {len(sig_data)}\n")

                        if len(sig_data) > 0:
                            f.write(f"\n  Top 10 positions by effect size:\n")
                            top_effects = sig_data.nlargest(10, 'delta')[
                                ['position', 'label', 'best_aa', 'best_mean', 'worst_aa',
                                 'worst_mean', 'delta', 'p_value']
                            ]
                            f.write(top_effects.to_string(index=False))
                            f.write("\n")

            # Substrate-specific residues
            if hasattr(self, 'substrate_specific'):
                f.write("\n" + "="*80 + "\n")
                f.write("SUBSTRATE-SPECIFIC RESIDUES\n")
                f.write("="*80 + "\n\n")
                f.write(f"Positions showing different preferences across substrates: {len(self.substrate_specific)}\n\n")

                if len(self.substrate_specific) > 0:
                    top_specific = self.substrate_specific.nlargest(20, 'max_delta')[
                        ['position', 'label', 'n_different_preferences', 'max_delta']
                    ]
                    f.write(top_specific.to_string(index=False))
                    f.write("\n")

            f.write("\n" + "="*80 + "\n")

        print(f"\nReport saved to: {report_file}")
        return report_file

    def plot_amino_acid_effects(self, output_dir='./'):
        """Generate visualizations"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        if not hasattr(self, 'aa_effects'):
            print("No amino acid effects data to plot")
            return

        # Plot 1: Heatmap of best amino acids by position and substrate
        for param in self.aa_effects['parameter'].unique():
            fig, axes = plt.subplots(1, 2, figsize=(18, 10))
            fig.suptitle(f'Amino Acid Effects on {param}', fontsize=16, fontweight='bold')

            param_data = self.aa_effects[self.aa_effects['parameter'] == param]

            # Heatmap of effect sizes
            ax1 = axes[0]
            pivot = param_data.pivot_table(
                values='delta',
                index='position',
                columns='substrate',
                aggfunc='mean'
            )

            sns.heatmap(pivot, cmap='RdYlGn', center=0, ax=ax1,
                       cbar_kws={'label': 'Effect Size (best - worst AA)'})
            ax1.set_title('Effect Size: Best vs Worst Amino Acid')
            ax1.set_xlabel('Substrate')
            ax1.set_ylabel('Position')

            # Significance count
            ax2 = axes[1]
            sig_counts = param_data[param_data['significant']].groupby('substrate').size()
            sig_counts.plot(kind='bar', ax=ax2, color='skyblue', edgecolor='black')
            ax2.set_title('Number of Significant Positions')
            ax2.set_xlabel('Substrate')
            ax2.set_ylabel('Count')
            ax2.grid(True, alpha=0.3)

            plt.tight_layout()
            filename = f"aa_effects_{param}_{timestamp}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
            plt.close()

        # Plot 2: Substrate comparison for top positions
        param = 'log_kcat_Km'
        if param in self.aa_effects['parameter'].unique():
            param_data = self.aa_effects[self.aa_effects['parameter'] == param]
            top_positions = param_data.groupby('position')['delta'].mean().nlargest(10).index

            fig, axes = plt.subplots(2, 5, figsize=(20, 8))
            fig.suptitle(f'Top 10 Positions: Substrate Comparison ({param})',
                        fontsize=16, fontweight='bold')

            for idx, pos in enumerate(top_positions):
                ax = axes[idx // 5, idx % 5]

                pos_data = param_data[param_data['position'] == pos]

                x = range(len(pos_data))
                colors = plt.cm.Set3(range(len(pos_data)))

                bars = ax.bar(x, pos_data['delta'], color=colors, edgecolor='black')
                ax.set_xticks(x)
                ax.set_xticklabels(pos_data['substrate'], rotation=45, ha='right', fontsize=8)
                ax.set_ylabel('Effect Size')
                ax.set_title(f'Position {pos}', fontweight='bold')
                ax.grid(True, alpha=0.3, axis='y')

                # Add amino acid labels
                for i, (_, row) in enumerate(pos_data.iterrows()):
                    ax.text(i, row['delta']/2, f"{row['best_aa']}",
                           ha='center', va='center', fontweight='bold', fontsize=10)

            plt.tight_layout()
            filename = f"substrate_comparison_top_positions_{timestamp}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
            plt.close()

    def plot_sample_sizes_by_position(self, output_dir='./'):
        """
        Plot amino acid sample size (n) distribution for each position
        Shows which amino acids are present and how many sequences have each AA
        """
        if self.combined_data is None:
            print("No position data. Run extract_position_data() first.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Get positions from position_info
        positions = sorted(self.position_info.keys())

        # Create figure with subplots for each position
        n_positions = len(positions)
        n_cols = 5
        n_rows = (n_positions + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(25, 5*n_rows))
        fig.suptitle('Amino Acid Sample Sizes by Position', fontsize=16, fontweight='bold')

        if n_rows == 1:
            axes = axes.reshape(1, -1)

        for idx, pos in enumerate(positions):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col]

            pos_col = f'pos_{pos+1}'  # Column names are 1-indexed
            info = self.position_info[pos]

            # Count amino acids across all substrates at this position
            if pos_col not in self.combined_data.columns:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax.set_title(f"Pos {pos+1}\n{info['label']}")
                continue

            aa_counts_all = self.combined_data[pos_col].value_counts().sort_index()

            if len(aa_counts_all) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax.set_title(f"Pos {pos+1}\n{info['label']}")
                continue

            # Create bar plot
            colors = ['red' if count == 1 else 'lightblue' if count < 5 else 'steelblue'
                     for count in aa_counts_all.values]

            bars = ax.bar(range(len(aa_counts_all)), aa_counts_all.values,
                         color=colors, edgecolor='black', alpha=0.7)

            ax.set_xticks(range(len(aa_counts_all)))
            ax.set_xticklabels(aa_counts_all.index, fontsize=8)
            ax.set_ylabel('Count (n)', fontsize=9)
            ax.set_title(f"Pos {pos+1}: {info['label']}\n({info['site']}, w={info['window']})",
                        fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')

            # Add value labels on bars
            for i, (aa, count) in enumerate(aa_counts_all.items()):
                if count > 0:
                    ax.text(i, count + max(aa_counts_all.values)*0.02, str(count),
                           ha='center', va='bottom', fontsize=7)

            # Rotate x labels if too many
            if len(aa_counts_all) > 10:
                ax.tick_params(axis='x', rotation=45)

        # Hide empty subplots
        for idx in range(n_positions, n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].axis('off')

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', edgecolor='black', alpha=0.7, label='n=1 (excluded)'),
            Patch(facecolor='lightblue', edgecolor='black', alpha=0.7, label='n=2-4 (low)'),
            Patch(facecolor='steelblue', edgecolor='black', alpha=0.7, label='n≥5 (usable)')
        ]
        fig.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout(rect=[0, 0, 1, 0.98])
        filename = f"sample_sizes_by_position_{timestamp}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved: {filename}")
        plt.close()

        # Also create a summary table
        summary_data = []
        for pos in positions:
            pos_col = f'pos_{pos+1}'  # Column names are 1-indexed
            info = self.position_info[pos]

            if pos_col not in self.combined_data.columns:
                continue

            aa_counts = self.combined_data[pos_col].value_counts()

            n_total = len(self.combined_data[pos_col].dropna())
            n_unique_aa = len(aa_counts)
            n_singleton = (aa_counts == 1).sum()
            n_usable = (aa_counts >= 5).sum()
            max_n = aa_counts.max() if len(aa_counts) > 0 else 0

            summary_data.append({
                'position': pos + 1,
                'label': info['label'],
                'site': info['site'],
                'window': info['window'],
                'total_seqs': n_total,
                'unique_AAs': n_unique_aa,
                'n=1 (excluded)': n_singleton,
                'n≥5 (usable)': n_usable,
                'max_n': max_n
            })

        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(output_dir, f"sample_size_summary_{timestamp}.csv")
        summary_df.to_csv(summary_file, index=False)
        print(f"Saved summary: sample_size_summary_{timestamp}.csv")

        return summary_df


def main():
    """Example usage"""
    print("="*80)
    print("SEQUENCE COMPARISON ANALYZER")
    print("="*80)

    # Define binding sites (0-indexed)
    binding_sites = {
        84: 'NAD_ribose',
        85: 'NAD_adenine',
        122: 'substrate_binding',
        155: 'catalytic_site',
        180: 'cofactor_binding'
    }

    # Initialize analyzer
    analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=3)

    # Define substrate files
    substrate_files = {
        'NADplus': 'catapro_predict_score_NADplus.csv',
        'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
        'alpha_L_Guluronate': 'catapro_predict_score_α-L-Guluronate.csv',
        'beta_D_Galacturonate': 'catapro_predict_score_β-D-Galacturonate.csv',
        'beta_D_Mannuronate': 'catapro_predict_score_β-D-Mannuronate.csv'
    }

    # Load data
    analyzer.load_substrate_data(substrate_files)

    # Extract position data
    analyzer.extract_position_data(aligned_seq_col='aligned_seq')

    # Analyze amino acid effects
    analyzer.analyze_amino_acid_effects(min_sequences=5)

    # Find substrate-specific residues
    analyzer.find_substrate_specific_residues(parameter='log_kcat_Km', min_delta=0.1)

    # Generate outputs
    analyzer.generate_report()
    analyzer.plot_amino_acid_effects()

    # Example: Compare specific position across substrates
    print("\n" + "="*80)
    print("Example: Position 85 comparison across substrates")
    print("="*80)
    comparison = analyzer.compare_substrates(position=85, parameter='log_kcat_Km')
    if comparison is not None:
        print(comparison.to_string(index=False))

    print("\n" + "="*80)


if __name__ == "__main__":
    main()
