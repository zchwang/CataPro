"""
Simple alignment parser without Biopython dependency
Parses CLUSTAL format alignment files
"""

import pandas as pd
import re
import os


def parse_clustal_alignment(aln_file):
    """
    Parse CLUSTAL format alignment file

    Returns:
    --------
    DataFrame with columns: fasta_id, aligned_seq
    """
    print(f"Parsing alignment from {aln_file}...")

    sequences = {}
    current_block = []

    with open(aln_file, 'r') as f:
        lines = f.readlines()

    # Skip header
    start_line = 0
    for i, line in enumerate(lines):
        if line.startswith('CLUSTAL'):
            start_line = i + 1
            break

    # Parse alignment blocks
    for line in lines[start_line:]:
        line = line.strip()

        # Skip empty lines and conservation lines
        if not line or line.startswith('*') or line.startswith(':') or line.startswith('.'):
            continue

        # Parse sequence line
        parts = line.split()
        if len(parts) >= 2:
            seq_id = parts[0]
            seq_part = parts[1]

            if seq_id not in sequences:
                sequences[seq_id] = []
            sequences[seq_id].append(seq_part)

    # Concatenate sequence parts
    final_sequences = {}
    for seq_id, parts in sequences.items():
        final_sequences[seq_id] = ''.join(parts)

    # Create DataFrame
    data = []
    for seq_id, seq in final_sequences.items():
        data.append({
            'fasta_id': seq_id,
            'aligned_seq': seq,
            'seq_length': len(seq.replace('-', ''))
        })

    df = pd.DataFrame(data)

    print(f"  Parsed {len(df)} sequences")
    print(f"  Alignment length: {len(df.iloc[0]['aligned_seq']) if len(df) > 0 else 0}")

    return df


def merge_with_catapro(alignment_df, catapro_file, substrate_name):
    """
    Merge alignment with CataPro predictions

    Parameters:
    -----------
    alignment_df : DataFrame
        DataFrame with fasta_id and aligned_seq
    catapro_file : str
        Path to CataPro prediction CSV
    substrate_name : str
        Name of substrate
    """
    print(f"\nMerging with {substrate_name} predictions...")

    if not os.path.exists(catapro_file):
        print(f"  Warning: {catapro_file} not found")
        return None

    # Load CataPro data
    catapro_df = pd.read_csv(catapro_file)

    if 'Unnamed: 0' in catapro_df.columns:
        catapro_df = catapro_df.drop('Unnamed: 0', axis=1)

    print(f"  CataPro: {len(catapro_df)} sequences")
    print(f"  Alignment: {len(alignment_df)} sequences")

    # Merge on fasta_id
    merged = alignment_df.merge(catapro_df, on='fasta_id', how='inner')

    print(f"  Merged: {len(merged)} sequences")

    if len(merged) == 0:
        print("  Warning: No matching sequences!")
        print(f"  Alignment IDs sample: {alignment_df['fasta_id'].head().tolist()}")
        print(f"  CataPro IDs sample: {catapro_df['fasta_id'].head().tolist()}")

    return merged


def extract_binding_site_info(pdb_file):
    """
    Simple extraction of residue numbers from PDB file

    Returns sequence and residue numbers
    """
    print(f"\nExtracting residue info from {pdb_file}...")

    residues = []
    residue_nums = []
    current_num = None

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Parse PDB ATOM line
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                chain = line[21]

                if res_num != current_num and chain == 'A':
                    residues.append(res_name)
                    residue_nums.append(res_num)
                    current_num = res_num

    print(f"  Found {len(residues)} residues in chain A")
    print(f"  Residue numbers: {min(residue_nums)} - {max(residue_nums)}")

    # Convert 3-letter to 1-letter amino acid code
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    sequence = ''.join([aa_map.get(res, 'X') for res in residues])
    print(f"  Sequence: {sequence[:50]}...")

    return sequence, residue_nums


def main():
    """Main execution"""
    print("="*80)
    print("SIMPLE ALIGNMENT PARSER")
    print("="*80)

    # Parse alignment
    aln_file = 'DWC0THCP211-alignment.aln'
    alignment_df = parse_clustal_alignment(aln_file)

    print(f"\nAlignment DataFrame head:")
    print(alignment_df.head())

    # Extract PDB info (optional)
    pdb_file = 'pdb3rfv_chainA.pdb'
    if os.path.exists(pdb_file):
        pdb_seq, residue_nums = extract_binding_site_info(pdb_file)
        print(f"\nPDB sequence length: {len(pdb_seq)}")
        print(f"Alignment length: {len(alignment_df.iloc[0]['aligned_seq'])}")

    # Define substrate files
    substrate_files = {
        'NADplus': 'catapro_predict_score_NADplus.csv',
        'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
        'alpha_L_Guluronate': 'catapro_predict_score_α-L-Guluronate.csv',
        'beta_D_Galacturonate': 'catapro_predict_score_β-D-Galacturonate.csv',
        'beta_D_Mannuronate': 'catapro_predict_score_β-D-Mannuronate.csv'
    }

    # Merge with CataPro predictions and save
    print("\n" + "="*80)
    print("Creating substrate files with aligned sequences...")
    print("="*80)

    created_files = []

    for substrate, catapro_file in substrate_files.items():
        merged = merge_with_catapro(alignment_df, catapro_file, substrate)

        if merged is not None and len(merged) > 0:
            output_file = f"catapro_with_alignment_{substrate}.csv"
            merged.to_csv(output_file, index=False)
            print(f"  Saved: {output_file}")
            created_files.append(output_file)

    print(f"\n✓ Created {len(created_files)} files")

    # Generate example binding sites
    print("\n" + "="*80)
    print("EXAMPLE: Defining Binding Sites")
    print("="*80)
    print("\nYou can manually define binding sites from clustering analysis:")
    print("\nExample code for sequence_comparison_analyzer.py:")
    print("""
binding_sites = {
    84: 'NAD_ribose',      # Position 85 in 1-indexed
    85: 'NAD_adenine',     # Position 86 in 1-indexed
    122: 'substrate',      # Position 123 in 1-indexed
    155: 'catalytic',      # Position 156 in 1-indexed
    180: 'cofactor'        # Position 181 in 1-indexed
}

analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=3)
    """)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"✓ Parsed alignment: {len(alignment_df)} sequences")
    print(f"✓ Created {len(created_files)} substrate files with aligned sequences")
    print(f"\nNext steps:")
    print(f"  1. Define your binding site positions from clustering")
    print(f"  2. Run sequence_comparison_analyzer.py with the created files")
    print("="*80)


if __name__ == "__main__":
    main()
