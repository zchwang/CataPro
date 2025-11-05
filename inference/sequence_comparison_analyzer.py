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


def create_pdb_to_alignment_mapping(reference_aligned_seq):
    """
    Create mapping from PDB residue numbers to alignment positions

    PDB structure uses sequential residue numbering (1, 2, 3, ...)
    But in MSA, the reference sequence may have gaps, so position != residue number

    Parameters:
    -----------
    reference_aligned_seq : str
        Aligned reference sequence with gaps

    Returns:
    --------
    dict : {pdb_residue_number: alignment_position}

    Example:
    --------
    ref_seq = "---MA-TG"  # alignment length = 8
    mapping = {1: 3, 2: 4, 3: 6, 4: 7}  # M=pos3, A=pos4, T=pos6, G=pos7
    """
    pdb_to_aln = {}
    pdb_residue = 1

    for aln_pos, char in enumerate(reference_aligned_seq):
        if char != '-':  # Skip gaps
            pdb_to_aln[pdb_residue] = aln_pos
            pdb_residue += 1

    return pdb_to_aln


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
            Binding site positions (0-indexed PDB residue numbers)
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
        self.pdb_to_aln = None  # Will be set after loading data

    def set_reference_sequence(self, reference_aligned_seq):
        """
        Create PDB to alignment position mapping from reference sequence

        Parameters:
        -----------
        reference_aligned_seq : str
            Aligned reference sequence with gaps (from alignment file)
        """
        print("\nCreating PDB → Alignment position mapping...")
        self.pdb_to_aln = create_pdb_to_alignment_mapping(reference_aligned_seq)

        # Show some key mappings for verification
        print(f"  Total PDB residues: {len(self.pdb_to_aln)}")
        print(f"  Alignment length: {len(reference_aligned_seq)}")

        # Check key positions
        test_positions = [1, 87, 110, 136, 159, 181]
        for pdb_pos in test_positions:
            if pdb_pos in self.pdb_to_aln:
                aln_pos = self.pdb_to_aln[pdb_pos]
                aa = reference_aligned_seq[aln_pos] if aln_pos < len(reference_aligned_seq) else '?'
                print(f"  PDB residue {pdb_pos:3d} → Alignment position {aln_pos:3d} ({aa})")

        return self.pdb_to_aln

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
        if self.pdb_to_aln is None:
            raise ValueError("Must call set_reference_sequence() before extract_position_data()")

        # First pass: count gaps at each position to filter out insertion-heavy positions
        print(f"\nFiltering positions by coverage (min_coverage={min_coverage})...")

        first_substrate = list(self.substrate_data.values())[0]
        n_sequences = len(first_substrate)

        positions_to_keep = {}
        positions_removed = []

        for pos, info in sorted(self.positions_with_window.items()):
            # pos is 0-indexed PDB position, need to convert to 1-indexed for mapping
            pdb_residue = pos + 1

            if pdb_residue not in self.pdb_to_aln:
                positions_removed.append({
                    'position': pdb_residue,
                    'label': info['label'],
                    'coverage': 0.0,
                    'gap_count': 'N/A',
                    'total_count': 'N/A',
                    'reason': 'PDB position not in reference sequence'
                })
                continue

            aln_pos = self.pdb_to_aln[pdb_residue]
            gap_count = 0
            total_count = 0

            for df in self.substrate_data.values():
                for idx, row in df.iterrows():
                    seq = row[aligned_seq_col]
                    if aln_pos >= 0 and aln_pos < len(seq):
                        total_count += 1
                        if seq[aln_pos] == '-':
                            gap_count += 1

            if total_count > 0:
                coverage = 1 - (gap_count / total_count)
                if coverage >= min_coverage:
                    positions_to_keep[pos] = info
                else:
                    positions_removed.append({
                        'position': pdb_residue,
                        'label': info['label'],
                        'coverage': coverage,
                        'gap_count': gap_count,
                        'total_count': total_count
                    })

        if positions_removed:
            print(f"\nRemoved {len(positions_removed)} positions with low coverage:")
            for p in positions_removed:
                if 'reason' in p:
                    print(f"  PDB {p['position']} ({p['label']}): {p['reason']}")
                else:
                    print(f"  PDB {p['position']} ({p['label']}): {p['coverage']*100:.1f}% coverage "
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
                    # pos is 0-indexed PDB position, convert to alignment position
                    pdb_residue = pos + 1
                    aln_pos = self.pdb_to_aln[pdb_residue]

                    if aln_pos >= 0 and aln_pos < len(seq):
                        aa = seq[aln_pos]
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

                    # Find best amino acid
                    best_aa = max(aa_means, key=aa_means.get)
                    best_data = substrate_data[substrate_data[pos_col] == best_aa][param].dropna()

                    # Get substrate original data for identity calculations
                    substrate_original = self.substrate_data[substrate]

                    # Calculate best AA identity
                    best_fasta_ids = substrate_data[substrate_data[pos_col] == best_aa]['fasta_id'].tolist()
                    best_sequences = substrate_original[substrate_original['fasta_id'].isin(best_fasta_ids)].drop_duplicates(subset='fasta_id')['aligned_seq'].tolist()
                    best_identity = calculate_average_pairwise_identity(best_sequences)
                    best_unique_count = len(best_sequences)

                    # Compare BEST AA with ALL OTHER AAs (not just worst)
                    # Store pairwise comparisons
                    aa_pairwise_pvalues = {}

                    for other_aa in valid_aas:
                        if other_aa == best_aa:
                            aa_pairwise_pvalues[other_aa] = 0.0  # Best vs itself
                            continue

                        other_data = substrate_data[substrate_data[pos_col] == other_aa][param].dropna()

                        if len(best_data) >= min_sequences and len(other_data) >= min_sequences:
                            try:
                                stat, pval = stats.mannwhitneyu(best_data, other_data, alternative='two-sided')
                                # Convert numpy types to Python types for JSON serialization
                                aa_pairwise_pvalues[other_aa] = float(pval) if not np.isnan(pval) else None
                            except:
                                aa_pairwise_pvalues[other_aa] = None
                        else:
                            aa_pairwise_pvalues[other_aa] = None

                    # Find worst amino acid for legacy compatibility
                    worst_aa = min(aa_means, key=aa_means.get)
                    worst_data = substrate_data[substrate_data[pos_col] == worst_aa][param].dropna()

                    # Calculate worst AA identity
                    worst_fasta_ids = substrate_data[substrate_data[pos_col] == worst_aa]['fasta_id'].tolist()
                    worst_sequences = substrate_original[substrate_original['fasta_id'].isin(worst_fasta_ids)].drop_duplicates(subset='fasta_id')['aligned_seq'].tolist()
                    worst_identity = calculate_average_pairwise_identity(worst_sequences)
                    worst_unique_count = len(worst_sequences)

                    # Statistical test between best and worst (for backward compatibility)
                    if len(best_data) >= min_sequences and len(worst_data) >= min_sequences:
                        stat, pval = stats.mannwhitneyu(best_data, worst_data, alternative='two-sided')
                        effect_size = (best_data.mean() - worst_data.mean()) / substrate_data[param].std()
                    else:
                        pval = np.nan
                        effect_size = np.nan

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
                        'best_unique_seqs': best_unique_count,
                        'best_seq_identity_mean': best_identity['mean_identity'],
                        'best_seq_identity_min': best_identity['min_identity'],
                        'best_seq_identity_max': best_identity['max_identity'],
                        'worst_aa': worst_aa,
                        'worst_mean': aa_means[worst_aa],
                        'worst_std': aa_stds[worst_aa],
                        'worst_count': aa_counts_dict[worst_aa],
                        'worst_unique_seqs': worst_unique_count,
                        'worst_seq_identity_mean': worst_identity['mean_identity'],
                        'worst_seq_identity_min': worst_identity['min_identity'],
                        'worst_seq_identity_max': worst_identity['max_identity'],
                        'delta': aa_means[best_aa] - aa_means[worst_aa],
                        'effect_size': effect_size,
                        'p_value': pval,
                        'significant': pval < 0.05 if not np.isnan(pval) else False,
                        'all_aa_means': str(aa_means),
                        'all_aa_counts': str(aa_counts_dict),
                        'all_aa_pvalues_vs_best': str(aa_pairwise_pvalues)  # NEW: p-values for all AAs vs Best
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

    def generate_excel_report(self, output_dir='./'):
        """
        Generate comprehensive Excel report with active site residues for each sequence

        Returns:
        --------
        str : Path to saved Excel file
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        excel_file = os.path.join(output_dir, f"active_site_residues_{timestamp}.xlsx")

        print("\nGenerating Excel report...")

        # Create Excel writer
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:

            # Sheet 1: Sequence Active Site Residues
            print("  Creating Sheet 1: Sequence Active Site Residues...")

            # Prepare data with position labels in header
            sequence_data = []

            for substrate in self.substrates:
                substrate_df = self.substrate_data[substrate]
                unique_df = substrate_df.drop_duplicates(subset='fasta_id')

                for idx, row in unique_df.iterrows():
                    seq_record = {
                        'fasta_id': row['fasta_id'],
                        'substrate': substrate,
                    }

                    # Add amino acids at each active site position
                    for pos, info in sorted(self.positions_with_window.items()):
                        pdb_residue = pos + 1
                        if pdb_residue in self.pdb_to_aln:
                            aln_pos = self.pdb_to_aln[pdb_residue]
                            seq = row['aligned_seq']
                            if aln_pos < len(seq):
                                aa = seq[aln_pos]
                                # Column name: "P12 (NAD_binding_1_core)"
                                col_name = f"P{pdb_residue} ({info['label']})"
                                seq_record[col_name] = aa if aa != '-' else None

                    # Add kinetic parameters
                    if 'pred_log10[kcat(s^-1)]' in substrate_df.columns:
                        seq_record['log_kcat'] = row['pred_log10[kcat(s^-1)]']
                    if 'pred_log10[Km(mM)]' in substrate_df.columns:
                        seq_record['log_Km'] = row['pred_log10[Km(mM)]']
                    if 'pred_log10[kcat/Km(s^-1mM^-1)]' in substrate_df.columns:
                        seq_record['log_kcat_Km'] = row['pred_log10[kcat/Km(s^-1mM^-1)]']

                    sequence_data.append(seq_record)

            sequence_df = pd.DataFrame(sequence_data)
            sequence_df.to_excel(writer, sheet_name='Active_Site_Residues', index=False)

            # Sheet 2: Position Information with Best/Worst per substrate
            print("  Creating Sheet 2: Position Information...")

            position_info_data = []
            for pos, info in sorted(self.positions_with_window.items()):
                pdb_residue = pos + 1

                base_info = {
                    'position': pdb_residue,
                    'label': info['label'],
                    'site': info['site'],
                    'window': info['window'],
                }

                # Add Best/Worst AA for each substrate from aa_effects
                if hasattr(self, 'aa_effects') and self.aa_effects is not None:
                    for substrate in self.substrates:
                        pos_effects = self.aa_effects[
                            (self.aa_effects['position'] == pdb_residue) &
                            (self.aa_effects['substrate'] == substrate) &
                            (self.aa_effects['parameter'] == 'log_kcat_Km')
                        ]

                        # Get significant effects first
                        pos_effects_sig = pos_effects[pos_effects['significant']]
                        if len(pos_effects_sig) > 0:
                            row = pos_effects_sig.iloc[0]
                        elif len(pos_effects) > 0:
                            row = pos_effects.iloc[0]
                        else:
                            continue

                        base_info[f'{substrate}_best_AA'] = row['best_aa']
                        base_info[f'{substrate}_best_mean'] = round(row['best_mean'], 3)
                        base_info[f'{substrate}_worst_AA'] = row['worst_aa']
                        base_info[f'{substrate}_worst_mean'] = round(row['worst_mean'], 3)
                        base_info[f'{substrate}_delta'] = round(row['delta'], 3)
                        base_info[f'{substrate}_p_value'] = row['p_value']
                        base_info[f'{substrate}_significant'] = row['significant']

                position_info_data.append(base_info)

            position_info_df = pd.DataFrame(position_info_data)
            position_info_df.to_excel(writer, sheet_name='Position_Information', index=False)

            # Sheet 3: Amino Acid Effects (full data)
            print("  Creating Sheet 3: Amino Acid Effects...")
            if hasattr(self, 'aa_effects') and self.aa_effects is not None:
                # Select relevant columns (check which columns actually exist)
                available_cols = ['position', 'label', 'site', 'substrate', 'parameter',
                                 'best_aa', 'best_mean', 'best_unique_seqs',
                                 'worst_aa', 'worst_mean', 'worst_unique_seqs',
                                 'delta', 'p_value', 'significant']

                # Only select columns that exist
                cols_to_export = [col for col in available_cols if col in self.aa_effects.columns]
                aa_effects_export = self.aa_effects[cols_to_export].copy()
                aa_effects_export.to_excel(writer, sheet_name='AA_Effects', index=False)

            # Sheet 4: Substrate-Specific Residues
            print("  Creating Sheet 4: Substrate-Specific Residues...")
            if hasattr(self, 'substrate_specific') and len(self.substrate_specific) > 0:
                substrate_specific_export = self.substrate_specific[[
                    'position', 'label', 'site', 'window', 'parameter',
                    'n_substrates', 'n_different_preferences', 'max_delta'
                ]].copy()
                substrate_specific_export.to_excel(writer, sheet_name='Substrate_Specific', index=False)

        print(f"\nExcel report saved: {excel_file}")
        return excel_file

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

            # Display image
            try:
                import subprocess
                # subprocess.Popen(['start', filepath], shell=True)
            except:
                pass

            plt.close()

        # Plot 2: Separate plots for EACH substrate showing Best vs Worst side-by-side with data points
        param = 'log_kcat_Km'
        if param in self.aa_effects['parameter'].unique():
            param_data = self.aa_effects[self.aa_effects['parameter'] == param]

            # Create separate figure for each substrate
            for substrate in self.substrates:
                substrate_data_sig = param_data[
                    (param_data['substrate'] == substrate) &
                    (param_data['significant'])
                ].sort_values('delta', ascending=False)

                substrate_data_insig = param_data[
                    (param_data['substrate'] == substrate) &
                    (~param_data['significant'])
                ].sort_values('delta', ascending=False)

                if len(substrate_data_sig) == 0 and len(substrate_data_insig) == 0:
                    continue

                # Sort by delta (effect size) - significant first, then insignificant
                all_data = pd.concat([substrate_data_sig, substrate_data_insig])

                # Get actual data
                substrate_combined = self.combined_data[self.combined_data['substrate'] == substrate]

                # Adjust figure size based on number of positions
                n_positions = len(all_data)
                fig_width = max(16, n_positions * 1.2)
                fig, ax = plt.subplots(figsize=(fig_width, 10))
                fig.suptitle(f'{substrate} - log(kcat/Km) Distribution: Best vs Worst\n' +
                            f'Green: Significant (p<0.05) | Gray: Not Significant (p≥0.05)',
                            fontsize=13, fontweight='bold')

                # Collect data for violin plots - Best and Worst side by side
                best_violin_data = []
                worst_violin_data = []
                x_positions_best = []
                x_positions_worst = []
                colors = []
                labels = []

                # Store data points for scatter
                best_scatter_data = []
                worst_scatter_data = []

                plot_idx = 0
                for idx, (_, row) in enumerate(all_data.iterrows()):
                    pos_col = f'pos_{row["position"]}'

                    if pos_col not in substrate_combined.columns:
                        continue

                    # Get actual log(kcat/Km) data for best and worst AA
                    best_data = substrate_combined[substrate_combined[pos_col] == row['best_aa']]['log_kcat_Km'].values
                    worst_data = substrate_combined[substrate_combined[pos_col] == row['worst_aa']]['log_kcat_Km'].values

                    if len(best_data) > 1 and len(worst_data) > 1:
                        # Position: Best on left, Worst on right (with gap)
                        x_best = plot_idx * 3
                        x_worst = plot_idx * 3 + 1

                        best_violin_data.append(best_data)
                        worst_violin_data.append(worst_data)
                        x_positions_best.append(x_best)
                        x_positions_worst.append(x_worst)
                        colors.append('green' if row['significant'] else 'lightgray')

                        # Store for scatter plot
                        best_scatter_data.append((x_best, best_data))
                        worst_scatter_data.append((x_worst, worst_data))

                        label = (f"P{row['position']}\n"
                                f"{row['best_aa']}({row['best_mean']:.2f})\n"
                                f"vs\n"
                                f"{row['worst_aa']}({row['worst_mean']:.2f})")
                        labels.append(label)
                        plot_idx += 1

                # Create violin plots for Best and Worst
                if len(best_violin_data) > 0:
                    # Plot Best violins
                    parts_best = ax.violinplot(best_violin_data, positions=x_positions_best, widths=0.8,
                                              showmeans=True, showmedians=True, showextrema=True)

                    # Plot Worst violins
                    parts_worst = ax.violinplot(worst_violin_data, positions=x_positions_worst, widths=0.8,
                                               showmeans=True, showmedians=True, showextrema=True)

                    # Color Best violins (darker shade)
                    for i, pc in enumerate(parts_best['bodies']):
                        if colors[i] == 'green':
                            pc.set_facecolor('darkgreen')
                        else:
                            pc.set_facecolor('gray')
                        pc.set_alpha(0.5)
                        pc.set_edgecolor('black')
                        pc.set_linewidth(1)

                    # Color Worst violins (lighter shade)
                    for i, pc in enumerate(parts_worst['bodies']):
                        if colors[i] == 'green':
                            pc.set_facecolor('lightgreen')
                        else:
                            pc.set_facecolor('lightgray')
                        pc.set_alpha(0.5)
                        pc.set_edgecolor('black')
                        pc.set_linewidth(1)

                    # Style means and medians for both
                    for parts in [parts_best, parts_worst]:
                        parts['cmeans'].set_color('red')
                        parts['cmeans'].set_linewidth(2)
                        parts['cmedians'].set_color('blue')
                        parts['cmedians'].set_linewidth(2)
                        parts['cbars'].set_color('black')
                        parts['cbars'].set_linewidth(1)
                        parts['cmaxes'].set_color('black')
                        parts['cmins'].set_color('black')

                    # Add data points as scatter (with jitter)
                    np.random.seed(42)
                    for x_pos, data in best_scatter_data:
                        # Add jitter to x position
                        jitter = np.random.normal(0, 0.05, size=len(data))
                        ax.scatter(x_pos + jitter, data, alpha=0.3, s=20, color='darkblue', zorder=3)

                    for x_pos, data in worst_scatter_data:
                        jitter = np.random.normal(0, 0.05, size=len(data))
                        ax.scatter(x_pos + jitter, data, alpha=0.3, s=20, color='darkred', zorder=3)

                    # Labels and formatting
                    # Set x-ticks at the midpoint between Best and Worst
                    x_ticks = [(x_positions_best[i] + x_positions_worst[i]) / 2 for i in range(len(x_positions_best))]
                    ax.set_xticks(x_ticks)
                    ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=7)
                    ax.set_ylabel('log(kcat/Km)', fontsize=12)
                    ax.set_xlabel('Position (sorted by Δ)\nDark=Best, Light=Worst', fontsize=12)
                    ax.set_title('Red=mean, Blue=median | Blue dots=Best, Red dots=Worst', fontsize=10)
                    ax.grid(True, alpha=0.3, axis='y')

                    # Add delta and p-values above each pair
                    for i, (_, row) in enumerate(all_data.iterrows()):
                        if i < len(x_ticks):
                            y_max = max(np.max(best_violin_data[i]), np.max(worst_violin_data[i]))

                            if row['significant']:
                                ax.text(x_ticks[i], y_max * 1.02, f"Δ={row['delta']:.3f}\np={row['p_value']:.4f}",
                                       ha='center', va='bottom', fontsize=6, color='darkred',
                                       fontweight='bold')
                            else:
                                ax.text(x_ticks[i], y_max * 1.02, f"Δ={row['delta']:.3f}\np={row['p_value']:.3f}",
                                       ha='center', va='bottom', fontsize=5, color='gray',
                                       style='italic')

                    # Add separation line between significant and insignificant
                    if len(substrate_data_sig) > 0 and len(substrate_data_insig) > 0:
                        sep_idx = len(substrate_data_sig) - 1
                        if sep_idx < len(x_ticks) - 1:
                            sep_x = x_ticks[sep_idx] + 1.5  # Between positions
                            ax.axvline(x=sep_x, color='red', linestyle='--', linewidth=2, alpha=0.5)
                            ax.text(sep_x, ax.get_ylim()[1] * 0.5, 'Significant | Insignificant',
                                   rotation=90, va='center', ha='right', fontsize=9, color='red',
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

                plt.tight_layout()
                filename = f"{substrate}_all_positions_{timestamp}.png"
                filepath = os.path.join(output_dir, filename)
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                print(f"Saved: {filename}")

                # Display image
                try:
                    import subprocess
                    # subprocess.Popen(['start', filepath], shell=True)
                except:
                    pass

                plt.close()

    def plot_sample_sizes_by_position(self, output_dir='./'):
        """
        Plot amino acid sample size (n) distribution for each position
        FOR EACH SUBSTRATE separately with substrate-specific significance
        """
        if self.combined_data is None:
            print("No position data. Run extract_position_data() first.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Get positions from position_info
        positions = sorted(self.position_info.keys())

        # Create SEPARATE figure for EACH substrate
        for substrate in self.substrates:
            n_positions = len(positions)
            n_cols = 5
            n_rows = (n_positions + n_cols - 1) // n_cols

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(25, 5*n_rows))
            fig.suptitle(f'{substrate} - Amino Acid Sample Sizes by Position\n(Unique Sequences for this substrate)',
                        fontsize=16, fontweight='bold')

            if n_rows == 1:
                axes = axes.reshape(1, -1)

            for idx, pos in enumerate(positions):
                row = idx // n_cols
                col = idx % n_cols
                ax = axes[row, col]

                pos_col = f'pos_{pos+1}'  # Column names are 1-indexed
                info = self.position_info[pos]

                # Count UNIQUE sequences FOR THIS SUBSTRATE ONLY
                # Extract AA from aligned_seq at position pos (0-indexed)
                unique_aa_counts = {}
                substrate_df = self.substrate_data[substrate]

                # Get unique sequences for this substrate
                unique_df = substrate_df.drop_duplicates(subset='fasta_id')

                for idx, row in unique_df.iterrows():
                    seq = row['aligned_seq']
                    # pos is 0-indexed PDB position, need to convert to alignment position
                    pdb_residue = pos + 1
                    if pdb_residue in self.pdb_to_aln:
                        aln_pos = self.pdb_to_aln[pdb_residue]
                        if aln_pos < len(seq):
                            aa = seq[aln_pos]
                            if aa != '-':  # Skip gaps
                                if aa not in unique_aa_counts:
                                    unique_aa_counts[aa] = 0
                                unique_aa_counts[aa] += 1

                if len(unique_aa_counts) == 0:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                    ax.set_title(f"Pos {pos+1}\n{info['label']}")
                    continue

                # Sort by amino acid letter
                aa_counts_series = pd.Series(unique_aa_counts).sort_index()

                # Find Best AA and get p-values for THIS SUBSTRATE
                best_aa = None
                worst_aa = None
                aa_pvalues = {}  # Store p-values for each AA vs Best

                if hasattr(self, 'aa_effects') and self.aa_effects is not None:
                    # Get effects for this position and THIS SUBSTRATE
                    pos_effects = self.aa_effects[
                        (self.aa_effects['position'] == pos+1) &
                        (self.aa_effects['substrate'] == substrate) &
                        (self.aa_effects['parameter'] == 'log_kcat_Km')
                    ]

                    # First try significant positions
                    pos_effects_sig = pos_effects[pos_effects['significant']]
                    if len(pos_effects_sig) > 0:
                        pos_effects = pos_effects_sig

                    if len(pos_effects) > 0:
                        # Take the first row
                        best_aa = pos_effects.iloc[0]['best_aa']
                        worst_aa = pos_effects.iloc[0]['worst_aa']

                        # Parse pairwise p-values
                        if 'all_aa_pvalues_vs_best' in pos_effects.columns:
                            import ast
                            try:
                                aa_pvalues = ast.literal_eval(pos_effects.iloc[0]['all_aa_pvalues_vs_best'])
                            except Exception as e:
                                aa_pvalues = {}

                # Create bar plot with colors indicating Best/Worst
                colors = []

                for aa, count in aa_counts_series.items():
                    if aa == best_aa:
                        colors.append('darkgreen')
                    elif aa == worst_aa:
                        colors.append('darkred')
                    elif count == 1:
                        colors.append('red')
                    elif count < 5:
                        colors.append('lightblue')
                    else:
                        colors.append('steelblue')

                bars = ax.bar(range(len(aa_counts_series)), aa_counts_series.values,
                             color=colors, edgecolor='black', alpha=0.7)

                ax.set_xticks(range(len(aa_counts_series)))
                ax.set_xticklabels(aa_counts_series.index, fontsize=8)
                ax.set_ylabel('Unique Seqs', fontsize=9)

                # Add Best/Worst to title if found
                title = f"Pos {pos+1}: {info['label']}\n({info['site']}, w={info['window']})"
                if best_aa and worst_aa:
                    title += f"\nBest:{best_aa} Worst:{worst_aa}"
                ax.set_title(title, fontsize=9)
                ax.grid(True, alpha=0.3, axis='y')

                # Add value labels on bars with p-value indicator
                for i, (aa, count) in enumerate(aa_counts_series.items()):
                    if count > 0:
                        label_text = str(count)
                        fontweight = 'normal'
                        color = 'black'

                        # Add significance indicator if this AA has a p-value vs Best
                        if aa == best_aa:
                            # Best AA - no label needed
                            pass
                        elif aa in aa_pvalues and aa_pvalues[aa] is not None:
                            if aa_pvalues[aa] < 0.05:
                                label_text += f"\n(p={aa_pvalues[aa]:.3f})"
                                fontweight = 'bold'
                                color = 'red'
                            else:
                                label_text += f"\n(ns)"
                                fontweight = 'normal'
                                color = 'gray'

                        ax.text(i, count + max(aa_counts_series.values)*0.02, label_text,
                               ha='center', va='bottom', fontsize=6, fontweight=fontweight, color=color)

                # Rotate x labels if too many
                if len(aa_counts_series) > 10:
                    ax.tick_params(axis='x', rotation=45)

            # Hide empty subplots
            for idx in range(n_positions, n_rows * n_cols):
                row = idx // n_cols
                col = idx % n_cols
                axes[row, col].axis('off')

            # Add legend
            from matplotlib.patches import Patch
            import matplotlib.lines as mlines
            legend_elements = [
                Patch(facecolor='darkgreen', edgecolor='black', alpha=0.7, label='Best AA'),
                Patch(facecolor='darkred', edgecolor='black', alpha=0.7, label='Worst AA'),
                mlines.Line2D([], [], color='red', marker='None', linestyle='None',
                             label='p<0.05 vs Best (red text)'),
                mlines.Line2D([], [], color='gray', marker='None', linestyle='None',
                             label='ns vs Best (gray text)'),
                Patch(facecolor='red', edgecolor='black', alpha=0.7, label='n=1 (excluded)'),
                Patch(facecolor='lightblue', edgecolor='black', alpha=0.7, label='n=2-4 (low)'),
                Patch(facecolor='steelblue', edgecolor='black', alpha=0.7, label='n>=5 (usable)')
            ]
            fig.legend(handles=legend_elements, loc='upper right', fontsize=9)

            plt.tight_layout(rect=[0, 0, 1, 0.98])
            filename = f"{substrate}_sample_sizes_by_position_{timestamp}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filename}")

            # Display image
            try:
                import subprocess
                # subprocess.Popen(['start', filepath], shell=True)
            except:
                pass

            plt.close()

        # Also create a summary table with identity information
        summary_data = []
        for pos in positions:
            info = self.position_info[pos]

            # Calculate unique sequence counts ACROSS ALL SUBSTRATES (no duplication)
            # Extract AA from aligned_seq at position pos (0-indexed)
            unique_aa_counts = {}
            aa_identity_stats = {}
            processed_fasta_ids = set()  # Track which sequences we've already counted

            for substrate, df in self.substrate_data.items():
                unique_df = df.drop_duplicates(subset='fasta_id')

                # Extract amino acid at position pos from aligned_seq
                for idx, row in unique_df.iterrows():
                    fasta_id = row['fasta_id']

                    # Skip if we've already counted this sequence from another substrate
                    if fasta_id in processed_fasta_ids:
                        continue

                    processed_fasta_ids.add(fasta_id)

                    seq = row['aligned_seq']
                    # pos is 0-indexed PDB position, need to convert to alignment position
                    pdb_residue = pos + 1
                    if pdb_residue in self.pdb_to_aln:
                        aln_pos = self.pdb_to_aln[pdb_residue]
                        if aln_pos < len(seq):
                            aa = seq[aln_pos]
                            if aa != '-':  # Skip gaps
                                if aa not in unique_aa_counts:
                                    unique_aa_counts[aa] = 0
                                unique_aa_counts[aa] += 1

            # Calculate identity for each AA group (using all substrates)
            for aa in unique_aa_counts.keys():
                # Get all unique sequences with this AA at this position
                aa_sequences = []
                seen_sequences = set()

                for substrate, df in self.substrate_data.items():
                    unique_df = df.drop_duplicates(subset='fasta_id')
                    for idx, row in unique_df.iterrows():
                        fasta_id = row['fasta_id']
                        if fasta_id not in seen_sequences:
                            seq = row['aligned_seq']
                            # pos is 0-indexed PDB position, convert to alignment position
                            pdb_residue = pos + 1
                            if pdb_residue in self.pdb_to_aln:
                                aln_pos = self.pdb_to_aln[pdb_residue]
                                if aln_pos < len(seq) and seq[aln_pos] == aa:
                                    aa_sequences.append(seq)
                                    seen_sequences.add(fasta_id)

                if len(aa_sequences) >= 2:
                    identity_info = calculate_average_pairwise_identity(aa_sequences)
                    aa_identity_stats[aa] = identity_info['mean_identity']
                elif len(aa_sequences) == 1:
                    aa_identity_stats[aa] = 100.0

            if len(unique_aa_counts) == 0:
                continue

            aa_counts_series = pd.Series(unique_aa_counts)
            n_total_unique = aa_counts_series.sum()
            n_unique_aa = len(aa_counts_series)
            n_singleton = (aa_counts_series == 1).sum()
            n_usable = (aa_counts_series >= 5).sum()
            max_n = aa_counts_series.max() if len(aa_counts_series) > 0 else 0

            # Average identity across all AA groups
            avg_identity = np.mean(list(aa_identity_stats.values())) if aa_identity_stats else 0.0

            summary_data.append({
                'position': pos + 1,
                'label': info['label'],
                'site': info['site'],
                'window': info['window'],
                'total_unique_seqs': n_total_unique,
                'unique_AAs': n_unique_aa,
                'n=1 (excluded)': n_singleton,
                'n>=5 (usable)': n_usable,
                'max_n': max_n,
                'avg_identity_%': round(avg_identity, 1)
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
