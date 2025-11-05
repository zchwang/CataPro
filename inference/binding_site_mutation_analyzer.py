"""
Binding Site-Focused Mutation Analysis
Analyzes mutations only at binding site residues with window context
"""

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os
import warnings
warnings.filterwarnings('ignore')


class BindingSiteMutationAnalyzer:
    """
    Analyzes mutations focusing on binding site residues
    """

    def __init__(self, binding_site_positions, window_size=3):
        """
        Parameters:
        -----------
        binding_site_positions : list or dict
            List of binding site positions, or dict mapping position to description
            Example: [85, 123, 156] or {85: 'NAD_binding', 123: 'substrate', 156: 'catalytic'}
        window_size : int
            Window size around each binding site position (1, 2, or 3)
        """
        if isinstance(binding_site_positions, dict):
            self.binding_sites = binding_site_positions
        else:
            self.binding_sites = {pos: f'site_{i+1}' for i, pos in enumerate(binding_site_positions)}

        self.window_size = window_size
        self.positions_with_window = self._expand_positions()

    def _expand_positions(self):
        """Expand binding site positions to include window"""
        position_info = {}

        for core_pos, site_name in self.binding_sites.items():
            # Core position (window 0)
            position_info[core_pos] = {
                'site': site_name,
                'window': 0,
                'core_position': core_pos,
                'label': f'{site_name}_core'
            }

            # Window positions
            for window in range(1, self.window_size + 1):
                # Upstream
                upstream_pos = core_pos - window
                position_info[upstream_pos] = {
                    'site': site_name,
                    'window': -window,
                    'core_position': core_pos,
                    'label': f'{site_name}_w-{window}'
                }

                # Downstream
                downstream_pos = core_pos + window
                position_info[downstream_pos] = {
                    'site': site_name,
                    'window': window,
                    'core_position': core_pos,
                    'label': f'{site_name}_w+{window}'
                }

        return position_info

    def get_analysis_summary(self):
        """Get summary of positions to be analyzed"""
        summary = {
            'n_core_sites': len(self.binding_sites),
            'window_size': self.window_size,
            'total_positions': len(self.positions_with_window),
            'positions_per_site': 1 + (2 * self.window_size),
            'core_positions': list(self.binding_sites.keys()),
            'all_positions': sorted(self.positions_with_window.keys())
        }
        return summary

    def load_and_filter_mutations(self, catapro_file, aligned_seq_col='aligned_seq'):
        """
        Load CataPro data and identify mutations at binding site positions

        Parameters:
        -----------
        catapro_file : str
            Path to CataPro prediction file with sequences
        aligned_seq_col : str
            Column name containing aligned sequences
        """
        print(f"Loading data from {catapro_file}...")
        df = pd.read_csv(catapro_file)

        if aligned_seq_col not in df.columns:
            raise ValueError(f"Column '{aligned_seq_col}' not found. Available: {df.columns.tolist()}")

        # Find reference (wild-type) sequence
        if 'type' in df.columns:
            wt_sequences = df[df['type'] == 'wild-type'][aligned_seq_col]
            if len(wt_sequences) > 0:
                reference_seq = wt_sequences.iloc[0]
            else:
                # Use consensus
                reference_seq = self._get_consensus_sequence(df[aligned_seq_col])
        else:
            reference_seq = self._get_consensus_sequence(df[aligned_seq_col])

        print(f"Reference sequence length: {len(reference_seq)}")
        print(f"Analyzing {len(self.positions_with_window)} positions (core + windows)")

        # Extract mutations at binding sites
        mutations_data = []

        for idx, row in df.iterrows():
            seq = row[aligned_seq_col]

            if len(seq) != len(reference_seq):
                print(f"Warning: Sequence {idx} length mismatch, skipping")
                continue

            # Check each binding site position
            seq_mutations = {}
            for pos, info in self.positions_with_window.items():
                if pos < 0 or pos >= len(reference_seq):
                    continue

                ref_aa = reference_seq[pos]
                mut_aa = seq[pos]

                if ref_aa != '-' and mut_aa != '-':
                    is_mutation = (ref_aa != mut_aa)
                    mutation_label = f"{ref_aa}{pos+1}{mut_aa}" if is_mutation else f"{ref_aa}{pos+1}{ref_aa}_WT"

                    seq_mutations[f"pos_{pos+1}"] = mut_aa
                    seq_mutations[f"pos_{pos+1}_mutation"] = mutation_label
                    seq_mutations[f"pos_{pos+1}_is_mutated"] = is_mutation
                    seq_mutations[f"pos_{pos+1}_site"] = info['site']
                    seq_mutations[f"pos_{pos+1}_window"] = info['window']
                    seq_mutations[f"pos_{pos+1}_label"] = info['label']

            # Combine with kinetic data
            row_data = row.to_dict()
            row_data.update(seq_mutations)
            row_data['n_mutations_in_sites'] = sum(v for k, v in seq_mutations.items() if k.endswith('_is_mutated'))
            mutations_data.append(row_data)

        self.mutations_df = pd.DataFrame(mutations_data)
        self.reference_seq = reference_seq

        print(f"\nAnalysis Summary:")
        print(f"Total sequences: {len(self.mutations_df)}")
        print(f"Sequences with mutations in binding sites: {(self.mutations_df['n_mutations_in_sites'] > 0).sum()}")
        print(f"Average mutations per sequence: {self.mutations_df['n_mutations_in_sites'].mean():.2f}")

        return self.mutations_df

    def _get_consensus_sequence(self, sequences):
        """Get consensus sequence from multiple sequences"""
        seq_array = np.array([list(seq) for seq in sequences])
        consensus = []
        for i in range(seq_array.shape[1]):
            col = seq_array[:, i]
            col = col[col != '-']
            if len(col) > 0:
                consensus.append(pd.Series(col).mode()[0])
            else:
                consensus.append('-')
        return ''.join(consensus)

    def analyze_single_mutations(self, kinetic_params=['pred_log10[kcat(s^-1)]',
                                                       'pred_log10[Km(mM)]',
                                                       'pred_log10[kcat/Km(s^-1mM^-1)]']):
        """
        Analyze effect of single mutations at each position
        """
        results = []

        for pos, info in sorted(self.positions_with_window.items()):
            pos_col = f"pos_{pos+1}"
            if pos_col not in self.mutations_df.columns:
                continue

            mutation_col = f"{pos_col}_mutation"
            is_mut_col = f"{pos_col}_is_mutated"

            if is_mut_col not in self.mutations_df.columns:
                continue

            # Compare mutated vs wild-type at this position
            wt_data = self.mutations_df[~self.mutations_df[is_mut_col]]
            mut_data = self.mutations_df[self.mutations_df[is_mut_col]]

            if len(mut_data) == 0:
                continue

            for param in kinetic_params:
                if param not in self.mutations_df.columns:
                    continue

                wt_values = wt_data[param].dropna()
                mut_values = mut_data[param].dropna()

                if len(wt_values) < 2 or len(mut_values) < 2:
                    continue

                # Statistical test
                stat, pval = stats.mannwhitneyu(wt_values, mut_values, alternative='two-sided')
                effect_size = (mut_values.mean() - wt_values.mean()) / wt_values.std()

                result = {
                    'position': pos + 1,
                    'site': info['site'],
                    'window': info['window'],
                    'label': info['label'],
                    'core_position': info['core_position'] + 1,
                    'parameter': param,
                    'n_wt': len(wt_values),
                    'n_mut': len(mut_values),
                    'wt_mean': wt_values.mean(),
                    'wt_std': wt_values.std(),
                    'mut_mean': mut_values.mean(),
                    'mut_std': mut_values.std(),
                    'delta': mut_values.mean() - wt_values.mean(),
                    'effect_size': effect_size,
                    'p_value': pval,
                    'significant': pval < 0.05
                }

                # Mutation details
                mutation_types = mut_data[mutation_col].value_counts()
                result['most_common_mutation'] = mutation_types.index[0] if len(mutation_types) > 0 else 'N/A'
                result['n_mutation_types'] = len(mutation_types)

                results.append(result)

        self.single_mutation_results = pd.DataFrame(results)
        return self.single_mutation_results

    def analyze_pairwise_interactions(self, kinetic_param='pred_log10[kcat/Km(s^-1mM^-1)]',
                                     min_samples=5):
        """
        Analyze pairwise interactions between binding site positions

        Parameters:
        -----------
        kinetic_param : str
            Which kinetic parameter to analyze
        min_samples : int
            Minimum number of samples required for each mutation combination
        """
        print("\nAnalyzing pairwise interactions...")

        # Get all position columns
        position_cols = [f"pos_{pos+1}_is_mutated" for pos in sorted(self.positions_with_window.keys())
                        if f"pos_{pos+1}_is_mutated" in self.mutations_df.columns]

        if len(position_cols) < 2:
            print("Not enough positions for pairwise analysis")
            return None

        pairwise_results = []
        total_pairs = len(list(combinations(position_cols, 2)))

        print(f"Analyzing {total_pairs} position pairs...")

        for i, (pos1_col, pos2_col) in enumerate(combinations(position_cols, 2)):
            if (i + 1) % 100 == 0:
                print(f"  Processing pair {i+1}/{total_pairs}")

            pos1 = int(pos1_col.split('_')[1])
            pos2 = int(pos2_col.split('_')[1])

            # Get position info
            info1 = self.positions_with_window.get(pos1 - 1, {})
            info2 = self.positions_with_window.get(pos2 - 1, {})

            # Create 4 groups: WT-WT, WT-Mut, Mut-WT, Mut-Mut
            wt_wt = self.mutations_df[(~self.mutations_df[pos1_col]) & (~self.mutations_df[pos2_col])][kinetic_param].dropna()
            wt_mut = self.mutations_df[(~self.mutations_df[pos1_col]) & (self.mutations_df[pos2_col])][kinetic_param].dropna()
            mut_wt = self.mutations_df[(self.mutations_df[pos1_col]) & (~self.mutations_df[pos2_col])][kinetic_param].dropna()
            mut_mut = self.mutations_df[(self.mutations_df[pos1_col]) & (self.mutations_df[pos2_col])][kinetic_param].dropna()

            # Check minimum samples
            if len(wt_wt) < min_samples or len(mut_mut) < min_samples:
                continue

            # Calculate interaction effect
            # Epistasis: observed double mutant effect vs expected (additive)
            if len(wt_mut) >= min_samples and len(mut_wt) >= min_samples:
                single_effect_1 = mut_wt.mean() - wt_wt.mean()
                single_effect_2 = wt_mut.mean() - wt_wt.mean()
                expected_double = wt_wt.mean() + single_effect_1 + single_effect_2
                observed_double = mut_mut.mean()
                epistasis = observed_double - expected_double

                # Statistical test for interaction
                try:
                    # ANOVA-like test
                    all_values = np.concatenate([wt_wt, wt_mut, mut_wt, mut_mut])
                    all_groups = ['WT-WT']*len(wt_wt) + ['WT-Mut']*len(wt_mut) + ['Mut-WT']*len(mut_wt) + ['Mut-Mut']*len(mut_mut)

                    _, p_interaction = stats.kruskal(wt_wt, wt_mut, mut_wt, mut_mut)

                except:
                    p_interaction = np.nan

                result = {
                    'pos1': pos1,
                    'pos2': pos2,
                    'site1': info1.get('site', 'unknown'),
                    'site2': info2.get('site', 'unknown'),
                    'window1': info1.get('window', 0),
                    'window2': info2.get('window', 0),
                    'label1': info1.get('label', f'pos{pos1}'),
                    'label2': info2.get('label', f'pos{pos2}'),
                    'n_wt_wt': len(wt_wt),
                    'n_wt_mut': len(wt_mut),
                    'n_mut_wt': len(mut_wt),
                    'n_mut_mut': len(mut_mut),
                    'mean_wt_wt': wt_wt.mean(),
                    'mean_wt_mut': wt_mut.mean(),
                    'mean_mut_wt': mut_wt.mean(),
                    'mean_mut_mut': mut_mut.mean(),
                    'single_effect_1': single_effect_1,
                    'single_effect_2': single_effect_2,
                    'expected_double': expected_double,
                    'observed_double': observed_double,
                    'epistasis': epistasis,
                    'p_value': p_interaction,
                    'significant': p_interaction < 0.05 if not np.isnan(p_interaction) else False
                }

                pairwise_results.append(result)

        self.pairwise_results = pd.DataFrame(pairwise_results)

        print(f"\nPairwise Analysis Complete:")
        print(f"Total pairs analyzed: {len(self.pairwise_results)}")
        if len(self.pairwise_results) > 0:
            print(f"Significant interactions (p<0.05): {self.pairwise_results['significant'].sum()}")
            print(f"Average epistasis: {self.pairwise_results['epistasis'].mean():.4f}")

        return self.pairwise_results

    def generate_report(self, output_dir='./'):
        """Generate comprehensive report"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = os.path.join(output_dir, f"binding_site_mutation_report_{timestamp}.txt")

        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("BINDING SITE MUTATION ANALYSIS REPORT\n")
            f.write("="*80 + "\n\n")

            # Summary
            summary = self.get_analysis_summary()
            f.write("ANALYSIS CONFIGURATION\n")
            f.write("-"*80 + "\n")
            f.write(f"Number of core binding sites: {summary['n_core_sites']}\n")
            f.write(f"Window size: ±{summary['window_size']} residues\n")
            f.write(f"Total positions analyzed: {summary['total_positions']}\n")
            f.write(f"Positions per site: {summary['positions_per_site']}\n\n")

            f.write("Core Binding Site Positions:\n")
            for pos, name in self.binding_sites.items():
                f.write(f"  Position {pos+1}: {name}\n")
            f.write("\n")

            # Data summary
            if hasattr(self, 'mutations_df'):
                f.write("DATA SUMMARY\n")
                f.write("-"*80 + "\n")
                f.write(f"Total sequences: {len(self.mutations_df)}\n")
                f.write(f"Sequences with mutations: {(self.mutations_df['n_mutations_in_sites'] > 0).sum()}\n")
                f.write(f"Average mutations per sequence: {self.mutations_df['n_mutations_in_sites'].mean():.2f}\n")
                f.write(f"Max mutations in one sequence: {self.mutations_df['n_mutations_in_sites'].max()}\n\n")

            # Single mutation results
            if hasattr(self, 'single_mutation_results'):
                f.write("SINGLE MUTATION ANALYSIS\n")
                f.write("-"*80 + "\n")
                sig_results = self.single_mutation_results[self.single_mutation_results['significant']]
                f.write(f"Total position-parameter combinations: {len(self.single_mutation_results)}\n")
                f.write(f"Significant effects (p<0.05): {len(sig_results)}\n\n")

                # Top effects by parameter
                for param in self.single_mutation_results['parameter'].unique():
                    param_data = sig_results[sig_results['parameter'] == param]
                    if len(param_data) > 0:
                        f.write(f"\nTop 10 effects for {param}:\n")
                        top_effects = param_data.nlargest(10, 'effect_size')[
                            ['position', 'site', 'window', 'label', 'delta', 'effect_size', 'p_value']
                        ]
                        f.write(top_effects.to_string(index=False))
                        f.write("\n")

            # Pairwise results
            if hasattr(self, 'pairwise_results') and len(self.pairwise_results) > 0:
                f.write("\n" + "="*80 + "\n")
                f.write("PAIRWISE INTERACTION ANALYSIS\n")
                f.write("="*80 + "\n\n")
                f.write(f"Total pairs analyzed: {len(self.pairwise_results)}\n")
                sig_pairs = self.pairwise_results[self.pairwise_results['significant']]
                f.write(f"Significant interactions (p<0.05): {len(sig_pairs)}\n\n")

                if len(sig_pairs) > 0:
                    f.write("Top 20 Epistatic Interactions:\n")
                    f.write("-"*80 + "\n")
                    top_epistasis = sig_pairs.nlargest(20, 'epistasis')[
                        ['pos1', 'pos2', 'label1', 'label2', 'epistasis', 'p_value']
                    ]
                    f.write(top_epistasis.to_string(index=False))
                    f.write("\n\n")

                    f.write("Top 20 Negative Epistatic Interactions:\n")
                    f.write("-"*80 + "\n")
                    bottom_epistasis = sig_pairs.nsmallest(20, 'epistasis')[
                        ['pos1', 'pos2', 'label1', 'label2', 'epistasis', 'p_value']
                    ]
                    f.write(bottom_epistasis.to_string(index=False))

            f.write("\n" + "="*80 + "\n")

        print(f"\nReport saved to: {report_file}")
        return report_file

    def plot_results(self, output_dir='./'):
        """Generate visualizations"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Plot 1: Single mutation effects
        if hasattr(self, 'single_mutation_results'):
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle('Binding Site Mutation Effects', fontsize=16, fontweight='bold')

            # Effect size by position
            ax1 = axes[0, 0]
            sig_data = self.single_mutation_results[self.single_mutation_results['significant']]
            if len(sig_data) > 0:
                pivot = sig_data.pivot_table(values='effect_size', index='position',
                                            columns='parameter', aggfunc='mean')
                pivot.plot(kind='bar', ax=ax1)
                ax1.set_title('Effect Size by Position')
                ax1.set_xlabel('Position')
                ax1.set_ylabel('Effect Size')
                ax1.legend(title='Parameter', bbox_to_anchor=(1.05, 1))
                ax1.grid(True, alpha=0.3)

            # Window analysis
            ax2 = axes[0, 1]
            window_effects = self.single_mutation_results.groupby('window')['effect_size'].mean()
            window_effects.plot(kind='bar', ax=ax2, color='skyblue')
            ax2.set_title('Average Effect Size by Window')
            ax2.set_xlabel('Window (0=core, ±1-3)')
            ax2.set_ylabel('Average Effect Size')
            ax2.grid(True, alpha=0.3)

            # Site comparison
            ax3 = axes[1, 0]
            site_effects = self.single_mutation_results.groupby('site')['effect_size'].mean().sort_values()
            site_effects.plot(kind='barh', ax=ax3, color='lightcoral')
            ax3.set_title('Average Effect Size by Binding Site')
            ax3.set_xlabel('Average Effect Size')
            ax3.grid(True, alpha=0.3)

            # P-value distribution
            ax4 = axes[1, 1]
            ax4.hist(self.single_mutation_results['p_value'], bins=50,
                    color='plum', edgecolor='black', alpha=0.7)
            ax4.axvline(0.05, color='red', linestyle='--', linewidth=2, label='p=0.05')
            ax4.set_title('P-value Distribution')
            ax4.set_xlabel('P-value')
            ax4.set_ylabel('Frequency')
            ax4.legend()
            ax4.grid(True, alpha=0.3)

            plt.tight_layout()
            filename = f"binding_site_single_mutations_{timestamp}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
            plt.close()

        # Plot 2: Pairwise interactions
        if hasattr(self, 'pairwise_results') and len(self.pairwise_results) > 0:
            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            fig.suptitle('Pairwise Epistatic Interactions', fontsize=16, fontweight='bold')

            # Epistasis heatmap (top interactions)
            ax1 = axes[0]
            sig_pairs = self.pairwise_results[self.pairwise_results['significant']]
            if len(sig_pairs) > 0:
                # Create matrix
                positions = sorted(set(sig_pairs['pos1'].tolist() + sig_pairs['pos2'].tolist()))
                matrix = pd.DataFrame(0.0, index=positions, columns=positions)

                for _, row in sig_pairs.iterrows():
                    matrix.loc[row['pos1'], row['pos2']] = row['epistasis']
                    matrix.loc[row['pos2'], row['pos1']] = row['epistasis']

                sns.heatmap(matrix, cmap='RdBu_r', center=0, ax=ax1,
                           cbar_kws={'label': 'Epistasis'})
                ax1.set_title('Epistatic Interaction Matrix')
                ax1.set_xlabel('Position')
                ax1.set_ylabel('Position')

            # Epistasis distribution
            ax2 = axes[1]
            ax2.hist(self.pairwise_results['epistasis'], bins=50,
                    color='skyblue', edgecolor='black', alpha=0.7)
            ax2.axvline(0, color='red', linestyle='--', linewidth=2)
            ax2.set_title('Epistasis Distribution')
            ax2.set_xlabel('Epistasis')
            ax2.set_ylabel('Frequency')
            ax2.grid(True, alpha=0.3)

            plt.tight_layout()
            filename = f"binding_site_pairwise_interactions_{timestamp}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")
            plt.close()


def main():
    """Example usage"""
    print("="*80)
    print("BINDING SITE MUTATION ANALYZER")
    print("="*80)

    # Example: Define binding sites
    # These should be actual binding site positions from your clustering analysis
    binding_sites = {
        84: 'NAD_ribose',      # Example positions (0-indexed)
        85: 'NAD_adenine',
        122: 'substrate_binding',
        155: 'catalytic_site',
        180: 'cofactor_binding'
    }

    window_size = 3  # Analyze ±3 residues around each site

    # Initialize analyzer
    analyzer = BindingSiteMutationAnalyzer(binding_sites, window_size=window_size)

    # Print summary
    summary = analyzer.get_analysis_summary()
    print(f"\nAnalysis Configuration:")
    print(f"  Core binding sites: {summary['n_core_sites']}")
    print(f"  Window size: ±{summary['window_size']}")
    print(f"  Total positions to analyze: {summary['total_positions']}")
    print(f"  Positions per site: {summary['positions_per_site']}")
    print(f"\nCore positions: {summary['core_positions']}")
    print(f"\nAll positions (including windows): {summary['all_positions']}")

    # Check if pairwise analysis is feasible
    n_pairs = len(list(combinations(summary['all_positions'], 2)))
    print(f"\nPairwise Analysis:")
    print(f"  Possible pairs: {n_pairs:,}")
    if n_pairs < 10000:
        print(f"  Feasibility: Recommended (manageable size)")
    else:
        print(f"  Feasibility: Large - may be slow")

    # Example with actual data (uncomment when you have the file)
    """
    # Load data
    df = analyzer.load_and_filter_mutations('catapro_predict_score_NADplus.csv',
                                            aligned_seq_col='aligned_seq')

    # Single mutation analysis
    single_results = analyzer.analyze_single_mutations()
    print(f"\nSingle mutation results: {len(single_results)} position-parameter combinations")

    # Pairwise analysis
    pairwise_results = analyzer.analyze_pairwise_interactions(min_samples=5)

    # Generate report and plots
    analyzer.generate_report()
    analyzer.plot_results()
    """

    print("\n" + "="*80)
    print("To use with your data, uncomment the example code in main() and provide:")
    print("  1. Actual binding site positions from clustering")
    print("  2. CataPro prediction file with aligned sequences")
    print("="*80)


if __name__ == "__main__":
    main()
