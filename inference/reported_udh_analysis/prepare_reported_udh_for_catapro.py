"""
Prepare reported UDH enzymes for CataPro prediction
Convert FASTA to CataPro input CSV format
"""

import pandas as pd

def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences"""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)

                # Get enzyme name
                current_name = line[1:].split('|')[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def create_catapro_input(sequences, substrates):
    """
    Create CataPro input CSV

    Format: Enzyme_id, type, sequence, smiles
    """
    rows = []

    for enzyme_name, seq in sequences.items():
        for substrate_name, smiles in substrates.items():
            rows.append({
                'Enzyme_id': f'{enzyme_name}_{substrate_name}',
                'type': 'wild-type',
                'sequence': seq,
                'smiles': smiles
            })

    return pd.DataFrame(rows)


def main():
    # SMILES from previous analysis
    substrates = {
        'NADplus': 'C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N',
        'alpha_D_Glucuronate': '[C@@H]1([C@@H]([C@@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O',
        'alpha_L_Guluronate': '[C@@H]1([C@H]([C@H](O[C@H]([C@@H]1O)O)C(=O)[O-])O)O',
        'beta_D_Galacturonate': '[C@@H]1([C@@H]([C@H](O[C@H]([C@H]1O)O)C(=O)[O-])O)O',
        'beta_D_Mannuronate': '[C@@H]1([C@@H]([C@H](O[C@@H]([C@@H]1O)O)C(=O)[O-])O)O'
    }

    print("Preparing reported UDH enzymes for CataPro prediction...")

    # 1. Load sequences
    print("\n1. Loading sequences from FASTA...")
    sequences = parse_fasta('reported_UDHs_complete.fasta')
    print(f"   Loaded {len(sequences)} enzymes")

    # 2. Create CataPro input
    print("\n2. Creating CataPro input CSV...")
    catapro_input = create_catapro_input(sequences, substrates)
    print(f"   Created {len(catapro_input)} rows ({len(sequences)} enzymes × {len(substrates)} substrates)")

    # 3. Save to CSV with index (required by CataPro)
    output_file = 'reported_udh_catapro_input.csv'
    catapro_input.to_csv(output_file, index=True)
    print(f"\n3. Saved to: {output_file}")

    print("\nReady for CataPro prediction!")
    print("Run: bash run_catapro.sh")


if __name__ == "__main__":
    main()
