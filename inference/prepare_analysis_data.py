"""
Prepare analysis data from alignment file, PDB, and CataPro predictions
"""

import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.PDB import PDBParser, DSSP, NeighborSearch
import os
import warnings
warnings.filterwarnings('ignore')


class DataPreparator:
    """
    Prepare data for sequence comparison analysis
    """

    def __init__(self):
        self.alignment = None
        self.pdb_structure = None
        self.binding_sites = {}

    def load_alignment(self, aln_file, format='clustal'):
        """
        Load alignment file (CLUSTAL, FASTA, etc.)

        Parameters:
        -----------
        aln_file : str
            Path to alignment file
        format : str
            Alignment format ('clustal', 'fasta', 'phylip', etc.)
        """
        print(f"Loading alignment from {aln_file}...")

        try:
            self.alignment = AlignIO.read(aln_file, format)
            print(f"  Loaded {len(self.alignment)} sequences")
            print(f"  Alignment length: {self.alignment.get_alignment_length()}")

            return self.alignment

        except Exception as e:
            print(f"Error loading alignment: {e}")
            print("Trying alternative formats...")

            # Try other formats
            for fmt in ['fasta', 'phylip', 'stockholm']:
                try:
                    self.alignment = AlignIO.read(aln_file, fmt)
                    print(f"  Successfully loaded as {fmt} format")
                    print(f"  Loaded {len(self.alignment)} sequences")
                    return self.alignment
                except:
                    continue

            raise ValueError(f"Could not load alignment file: {aln_file}")

    def alignment_to_dataframe(self):
        """
        Convert alignment to pandas DataFrame

        Returns:
        --------
        DataFrame with columns: fasta_id, aligned_seq
        """
        if self.alignment is None:
            raise ValueError("No alignment loaded. Call load_alignment() first.")

        data = []
        for record in self.alignment:
            data.append({
                'fasta_id': record.id,
                'aligned_seq': str(record.seq),
                'seq_length': len(str(record.seq).replace('-', ''))
            })

        df = pd.DataFrame(data)
        print(f"\nAlignment DataFrame:")
        print(f"  Sequences: {len(df)}")
        print(f"  Alignment length: {len(df.iloc[0]['aligned_seq'])}")

        return df

    def load_pdb_structure(self, pdb_file):
        """
        Load PDB structure

        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        """
        print(f"\nLoading PDB structure from {pdb_file}...")

        parser = PDBParser(QUIET=True)
        structure_id = os.path.basename(pdb_file).replace('.pdb', '')
        self.pdb_structure = parser.get_structure(structure_id, pdb_file)

        # Get chain info
        model = self.pdb_structure[0]
        chains = list(model.get_chains())

        print(f"  Structure: {structure_id}")
        print(f"  Chains: {[c.id for c in chains]}")
        print(f"  Residues: {sum(len(list(c.get_residues())) for c in chains)}")

        return self.pdb_structure

    def extract_binding_sites_from_pdb(self, ligand_name='NAD', distance_cutoff=5.0, chain_id='A'):
        """
        Extract binding site residues from PDB based on proximity to ligand

        Parameters:
        -----------
        ligand_name : str
            Name of ligand (e.g., 'NAD', 'NAP', 'GLC')
        distance_cutoff : float
            Distance cutoff in Angstroms
        chain_id : str
            Chain ID to analyze
        """
        if self.pdb_structure is None:
            raise ValueError("No PDB structure loaded")

        print(f"\nExtracting binding sites for {ligand_name}...")

        model = self.pdb_structure[0]
        chain = model[chain_id]

        # Find ligand
        ligand_atoms = []
        for residue in chain:
            if residue.get_resname() == ligand_name:
                ligand_atoms.extend(residue.get_atoms())

        if not ligand_atoms:
            print(f"  Warning: Ligand {ligand_name} not found in chain {chain_id}")
            print(f"  Available hetero residues:")
            hetero_res = set()
            for residue in chain:
                if residue.id[0] != ' ':  # Hetero residue
                    hetero_res.add(residue.get_resname())
            print(f"    {sorted(hetero_res)}")
            return {}

        print(f"  Found ligand {ligand_name} with {len(ligand_atoms)} atoms")

        # Find nearby protein residues
        protein_atoms = []
        for residue in chain:
            if residue.id[0] == ' ':  # Protein residue
                protein_atoms.extend(residue.get_atoms())

        ns = NeighborSearch(protein_atoms)

        binding_residues = set()
        for ligand_atom in ligand_atoms:
            nearby = ns.search(ligand_atom.coord, distance_cutoff, level='R')
            binding_residues.update(nearby)

        # Create binding site dictionary
        binding_sites = {}
        for residue in sorted(binding_residues, key=lambda r: r.id[1]):
            res_num = residue.id[1]
            res_name = residue.get_resname()
            binding_sites[res_num - 1] = f"{ligand_name}_{res_num}"  # 0-indexed

        print(f"  Found {len(binding_sites)} binding site residues")
        print(f"  Positions (0-indexed): {sorted(binding_sites.keys())[:20]}...")

        self.binding_sites[ligand_name] = binding_sites
        return binding_sites

    def manually_define_binding_sites(self, positions, site_name='binding_site'):
        """
        Manually define binding site positions

        Parameters:
        -----------
        positions : list
            List of positions (0-indexed or 1-indexed, specify below)
        site_name : str
            Name for this binding site
        one_indexed : bool
            If True, convert from 1-indexed to 0-indexed
        """
        binding_sites = {}
        for i, pos in enumerate(positions):
            binding_sites[pos] = f"{site_name}_{i+1}"

        self.binding_sites[site_name] = binding_sites
        print(f"\nManually defined {len(binding_sites)} positions for {site_name}")
        print(f"  Positions: {sorted(binding_sites.keys())}")

        return binding_sites

    def merge_with_catapro_predictions(self, alignment_df, catapro_file,
                                      on='fasta_id', how='inner'):
        """
        Merge alignment data with CataPro predictions

        Parameters:
        -----------
        alignment_df : DataFrame
            DataFrame with fasta_id and aligned_seq
        catapro_file : str
            Path to CataPro prediction file
        on : str
            Column to merge on
        how : str
            Merge type ('inner', 'left', 'outer')
        """
        print(f"\nMerging with CataPro predictions from {catapro_file}...")

        catapro_df = pd.read_csv(catapro_file)

        if 'Unnamed: 0' in catapro_df.columns:
            catapro_df = catapro_df.drop('Unnamed: 0', axis=1)

        print(f"  CataPro data: {len(catapro_df)} sequences")

        merged = alignment_df.merge(catapro_df, on=on, how=how)

        print(f"  Merged data: {len(merged)} sequences")
        print(f"  Columns: {merged.columns.tolist()}")

        return merged

    def create_substrate_files(self, alignment_df, catapro_files_dict, output_dir='./'):
        """
        Create separate files for each substrate with aligned sequences

        Parameters:
        -----------
        alignment_df : DataFrame
            DataFrame with fasta_id and aligned_seq
        catapro_files_dict : dict
            {substrate_name: catapro_file_path}
        output_dir : str
            Output directory
        """
        print("\nCreating substrate-specific files with aligned sequences...")

        created_files = {}

        for substrate, catapro_file in catapro_files_dict.items():
            print(f"\n  Processing {substrate}...")

            if not os.path.exists(catapro_file):
                print(f"    Warning: {catapro_file} not found, skipping")
                continue

            # Merge alignment with CataPro predictions
            merged = self.merge_with_catapro_predictions(
                alignment_df, catapro_file, on='fasta_id', how='inner'
            )

            # Save
            output_file = os.path.join(output_dir, f"catapro_with_alignment_{substrate}.csv")
            merged.to_csv(output_file, index=False)

            print(f"    Saved: {output_file} ({len(merged)} sequences)")
            created_files[substrate] = output_file

        print(f"\nCreated {len(created_files)} substrate files")
        return created_files

    def generate_summary_report(self, output_file='data_preparation_summary.txt'):
        """Generate summary report"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("DATA PREPARATION SUMMARY\n")
            f.write("="*80 + "\n\n")

            # Alignment info
            if self.alignment is not None:
                f.write("ALIGNMENT\n")
                f.write("-"*80 + "\n")
                f.write(f"Sequences: {len(self.alignment)}\n")
                f.write(f"Alignment length: {self.alignment.get_alignment_length()}\n")
                f.write(f"First sequence ID: {self.alignment[0].id}\n\n")

            # PDB info
            if self.pdb_structure is not None:
                f.write("PDB STRUCTURE\n")
                f.write("-"*80 + "\n")
                model = self.pdb_structure[0]
                chains = list(model.get_chains())
                f.write(f"Structure ID: {self.pdb_structure.id}\n")
                f.write(f"Chains: {[c.id for c in chains]}\n")
                f.write(f"Total residues: {sum(len(list(c.get_residues())) for c in chains)}\n\n")

            # Binding sites
            if self.binding_sites:
                f.write("BINDING SITES\n")
                f.write("-"*80 + "\n")
                for site_name, positions in self.binding_sites.items():
                    f.write(f"\n{site_name}:\n")
                    f.write(f"  Positions (0-indexed): {sorted(positions.keys())}\n")
                    f.write(f"  Total: {len(positions)} residues\n")

            f.write("\n" + "="*80 + "\n")

        print(f"\nSummary saved to: {output_file}")


def main():
    """Example usage"""
    print("="*80)
    print("DATA PREPARATION FOR SEQUENCE COMPARISON")
    print("="*80)

    # Initialize
    prep = DataPreparator()

    # Load alignment
    aln_file = 'DWC0THCP211-alignment.aln'
    prep.load_alignment(aln_file, format='clustal')

    # Convert to DataFrame
    alignment_df = prep.alignment_to_dataframe()

    # Load PDB (optional - for automatic binding site detection)
    pdb_file = 'pdb3rfv_chainA.pdb'
    if os.path.exists(pdb_file):
        prep.load_pdb_structure(pdb_file)

        # Try to extract binding sites from PDB
        # Try common ligand names
        for ligand in ['NAD', 'NAP', 'NDP', 'NAI']:
            sites = prep.extract_binding_sites_from_pdb(ligand, distance_cutoff=5.0)
            if sites:
                break

    # Or manually define binding sites
    # (Replace these with your actual binding site positions from clustering)
    print("\n" + "="*80)
    print("OPTION: Manually define binding sites")
    print("="*80)
    print("If PDB extraction didn't work or you want to use clustering positions:")
    print("Edit the script and uncomment the following lines:")
    print()
    print("# Example binding sites (0-indexed)")
    print("# prep.manually_define_binding_sites([84, 85, 122, 155, 180], 'NAD_binding')")

    # Create substrate files with aligned sequences
    catapro_files = {
        'NADplus': 'catapro_predict_score_NADplus.csv',
        'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
        'alpha_L_Guluronate': 'catapro_predict_score_α-L-Guluronate.csv',
        'beta_D_Galacturonate': 'catapro_predict_score_β-D-Galacturonate.csv',
        'beta_D_Mannuronate': 'catapro_predict_score_β-D-Mannuronate.csv'
    }

    # Check which files exist
    existing_files = {k: v for k, v in catapro_files.items() if os.path.exists(v)}

    if existing_files:
        print("\n" + "="*80)
        print("Creating substrate files with aligned sequences...")
        print("="*80)

        created_files = prep.create_substrate_files(alignment_df, existing_files)

        print("\n✓ Preparation complete!")
        print(f"Created {len(created_files)} files ready for analysis:")
        for substrate, filepath in created_files.items():
            print(f"  - {filepath}")

    else:
        print("\nWarning: No CataPro prediction files found")
        print("Expected files:", list(catapro_files.values()))

    # Generate summary
    prep.generate_summary_report()

    print("\n" + "="*80)
    print("Next steps:")
    print("  1. Review data_preparation_summary.txt")
    print("  2. If needed, manually define binding sites in the script")
    print("  3. Run sequence comparison analysis with prepared files")
    print("="*80)


if __name__ == "__main__":
    main()
