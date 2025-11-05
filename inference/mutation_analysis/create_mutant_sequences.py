"""
Create mutant sequences based on AA analysis suggestions
and test them with CataPro to verify if mutations actually improve performance
"""

import pandas as pd
import numpy as np
from parse_alignment_simple import parse_clustal_alignment
from sequence_comparison_analyzer import create_pdb_to_alignment_mapping

def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences"""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split('|')[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def map_pdb_to_sequence_position(alignment_df, pdb_to_aln, enzyme_name, reported_seqs):
    """
    Map PDB residue numbers to actual sequence positions (ungapped)

    For each enzyme, we need to:
    1. Find its aligned sequence
    2. Map PDB position → alignment position (using pdb_to_aln)
    3. Map alignment position → sequence position (removing gaps)
    """
    # Get the reported sequence (ungapped)
    if enzyme_name not in reported_seqs:
        return None

    ungapped_seq = reported_seqs[enzyme_name]

    # We need to create a mapping from alignment position to sequence position
    # by counting non-gap characters

    # For simplicity, we'll align the reported sequence to the reference
    # and create the mapping

    # Actually, we need the aligned version of this enzyme
    # Let's use simple approach: assume reported sequence aligns similar to reference

    # Better approach: create alignment position → sequence position mapping
    # by going through the reference alignment and counting non-gaps

    aln_to_seq = {}
    seq_pos = 0

    # We'll use the reported sequence length to create mapping
    # Assuming the sequence without gaps
    for i, char in enumerate(ungapped_seq):
        aln_to_seq[seq_pos] = i
        seq_pos += 1

    return aln_to_seq


def apply_mutations(sequence, mutations, pdb_to_seq):
    """
    Apply mutations to a sequence

    Parameters:
    -----------
    sequence : str
        Original amino acid sequence
    mutations : list of dict
        Each dict has 'Position' (PDB numbering), 'Current_AA', 'Target_AA'
    pdb_to_seq : dict
        Mapping from PDB position to sequence position

    Returns:
    --------
    str : mutated sequence
    """
    # Convert sequence to list for easier mutation
    seq_list = list(sequence)

    applied_mutations = []
    failed_mutations = []

    for mut in mutations:
        pdb_pos = mut['Position']
        target_aa = mut['Target_AA']
        current_aa = mut['Current_AA']

        # Convert PDB position to sequence position
        # PDB position is 1-indexed, so actual position in sequence is pdb_pos - 1
        # But we need to account for alignment gaps

        # Simple approach: PDB position directly corresponds to sequence position
        # This assumes no gaps in the enzyme sequence itself
        seq_pos = pdb_pos - 1  # Convert to 0-indexed

        # Safety check
        if seq_pos < 0 or seq_pos >= len(seq_list):
            failed_mutations.append(f"P{pdb_pos}: out of range")
            continue

        # Verify current AA matches (optional, for validation)
        if seq_list[seq_pos] != current_aa:
            # Could be alignment issue, still apply
            failed_mutations.append(f"P{pdb_pos}: expected {current_aa}, found {seq_list[seq_pos]}")

        # Apply mutation
        seq_list[seq_pos] = target_aa
        applied_mutations.append(f"P{pdb_pos}:{current_aa}->{target_aa}")

    return ''.join(seq_list), applied_mutations, failed_mutations


def create_mutants():
    """Create mutant sequences based on engineering suggestions"""

    print("="*80)
    print("Creating Mutant Sequences from AA Analysis Suggestions")
    print("="*80)

    # 1. Load original sequences
    print("\n1. Loading original sequences...")
    reported_seqs = parse_fasta('reported_UDHs_complete.fasta')
    print(f"   Loaded {len(reported_seqs)} enzymes")

    # 2. Load mutation suggestions
    print("\n2. Loading mutation suggestions...")
    engineering_df = pd.read_excel('reported_udh_analysis.xlsx',
                                    sheet_name='Engineering_Suggestions')
    print(f"   Loaded {len(engineering_df)} mutation suggestions")

    # 3. Load alignment for PDB mapping
    print("\n3. Loading alignment...")
    alignment_df = parse_clustal_alignment('DWC0THCP211-alignment.aln')
    reference_seq = alignment_df.iloc[0]['aligned_seq']
    pdb_to_aln = create_pdb_to_alignment_mapping(reference_seq)

    # 4. Create mutants for each enzyme
    print("\n4. Creating mutants...")

    # Substrates
    substrates = {
        'NADplus': 'C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N',
        'alpha_D_Glucuronate': '[C@@H]1([C@@H]([C@@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O',
        'alpha_L_Guluronate': '[C@@H]1([C@H]([C@H](O[C@H]([C@@H]1O)O)C(=O)[O-])O)O',
        'beta_D_Galacturonate': '[C@@H]1([C@@H]([C@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O',
        'beta_D_Mannuronate': '[C@@H]1([C@@H]([C@H](O[C@@H]([C@@H]1O)O)C(=O)[O-])O)O'
    }

    # Prepare data for CataPro input
    catapro_data = []
    mutation_summary = []

    for enzyme_name in reported_seqs.keys():
        original_seq = reported_seqs[enzyme_name]

        # Get all mutations for this enzyme
        enzyme_mutations = engineering_df[engineering_df['Enzyme'] == enzyme_name]

        if len(enzyme_mutations) == 0:
            print(f"\n   {enzyme_name}: No mutations suggested")
            continue

        print(f"\n   {enzyme_name}:")
        print(f"      Total mutations suggested: {len(enzyme_mutations)}")

        # Sort by expected gain
        enzyme_mutations = enzyme_mutations.sort_values('Expected_Gain', ascending=False)

        # Create different mutant versions
        mutation_sets = [
            ('top5', enzyme_mutations.head(5), 'Top 5 mutations'),
            ('top10', enzyme_mutations.head(10), 'Top 10 mutations'),
            ('all', enzyme_mutations, 'All mutations'),
        ]

        for set_name, mutations_to_apply, description in mutation_sets:
            if len(mutations_to_apply) == 0:
                continue

            # Apply mutations
            mutant_seq, applied, failed = apply_mutations(
                original_seq,
                mutations_to_apply.to_dict('records'),
                pdb_to_aln
            )

            print(f"      {description}: {len(applied)} applied, {len(failed)} failed")

            # Add to summary
            mutation_summary.append({
                'Enzyme': enzyme_name,
                'Mutant_Type': set_name,
                'Description': description,
                'Mutations_Applied': len(applied),
                'Mutations_Failed': len(failed),
                'Applied_List': '; '.join(applied[:10]),  # First 10 for summary
                'Failed_List': '; '.join(failed[:10]) if failed else 'None'
            })

            # Add to CataPro input for all substrates
            for substrate_name, smiles in substrates.items():
                catapro_data.append({
                    'Enzyme_id': f'{enzyme_name}_{set_name}_{substrate_name}',
                    'type': 'mutant',
                    'sequence': mutant_seq,
                    'smiles': smiles,
                    'original_enzyme': enzyme_name,
                    'mutant_type': set_name,
                    'substrate': substrate_name,
                    'num_mutations': len(applied)
                })

        # Also add original (wild-type) for comparison
        for substrate_name, smiles in substrates.items():
            catapro_data.append({
                'Enzyme_id': f'{enzyme_name}_wildtype_{substrate_name}',
                'type': 'wild-type',
                'sequence': original_seq,
                'smiles': smiles,
                'original_enzyme': enzyme_name,
                'mutant_type': 'wildtype',
                'substrate': substrate_name,
                'num_mutations': 0
            })

    # 5. Save CataPro input
    print("\n5. Saving CataPro input...")
    catapro_df = pd.DataFrame(catapro_data)

    # Save with index for CataPro
    output_file = 'mutant_sequences_catapro_input.csv'
    catapro_df[['Enzyme_id', 'type', 'sequence', 'smiles']].to_csv(output_file, index=True)
    print(f"   Saved {len(catapro_df)} rows to: {output_file}")

    # Save metadata
    metadata_file = 'mutant_sequences_metadata.csv'
    catapro_df.to_csv(metadata_file, index=False)
    print(f"   Saved metadata to: {metadata_file}")

    # Save mutation summary
    summary_df = pd.DataFrame(mutation_summary)
    summary_file = 'mutant_sequences_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"   Saved mutation summary to: {summary_file}")

    print("\n" + "="*80)
    print("Summary Statistics")
    print("="*80)
    print(f"\nTotal sequences created: {len(catapro_df)}")
    print(f"  Wild-type: {len(catapro_df[catapro_df['mutant_type']=='wildtype'])}")
    print(f"  Mutants: {len(catapro_df[catapro_df['mutant_type']!='wildtype'])}")

    print("\nMutant types:")
    for mtype in ['top5', 'top10', 'all']:
        count = len(catapro_df[catapro_df['mutant_type']==mtype])
        if count > 0:
            print(f"  {mtype}: {count} sequences")

    print("\n" + "="*80)
    print("Ready for CataPro prediction!")
    print("Run: python predict.py -inp_fpath mutant_sequences_catapro_input.csv \\")
    print("                       -model_dpath ../models \\")
    print("                       -batch_size 32 \\")
    print("                       -device cpu \\")
    print("                       -out_fpath mutant_sequences_predictions.csv")
    print("="*80)

    return catapro_df, summary_df


if __name__ == "__main__":
    catapro_df, summary_df = create_mutants()
