"""
Analyze reported UDH enzymes to find which ones are closest to optimal
and suggest engineering strategies
"""

import pandas as pd
import numpy as np
from parse_alignment_simple import parse_clustal_alignment
import os


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

                # Get enzyme name (e.g., "AtUdh", "PsUdh")
                current_name = line[1:].split('|')[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def load_best_amino_acids(excel_file, substrates=['NADplus', 'alpha_D_Glucuronate',
                                                   'alpha_L_Guluronate', 'beta_D_Galacturonate',
                                                   'beta_D_Mannuronate']):
    """
    Load best amino acids for each position from Excel file

    Returns:
    --------
    dict : {substrate: {position: {'best_aa': X, 'worst_aa': Y, 'delta': Z}}}
    """
    # Read Position_Information sheet
    position_info = pd.read_excel(excel_file, sheet_name='Position_Information')

    best_aas = {}

    for substrate in substrates:
        best_aas[substrate] = {}

        for _, row in position_info.iterrows():
            pos = row['position']

            # Get Best/Worst for this substrate
            best_col = f'{substrate}_best_AA'
            worst_col = f'{substrate}_worst_AA'
            delta_col = f'{substrate}_delta'
            sig_col = f'{substrate}_significant'

            if best_col in row and pd.notna(row[best_col]):
                best_aas[substrate][pos] = {
                    'best_aa': row[best_col],
                    'worst_aa': row[worst_col] if worst_col in row else None,
                    'delta': row[delta_col] if delta_col in row else None,
                    'significant': row[sig_col] if sig_col in row else False
                }

    return best_aas


def map_reported_udh_to_alignment(reported_seq, alignment_df, reference_id='lcl|Query_10001'):
    """
    Map reported UDH sequence to alignment by finding best match

    Returns:
    --------
    str : aligned sequence with gaps
    """
    # For simplicity, we'll align using the reference sequence
    # In practice, this should use proper alignment tools

    # Find the reference sequence in alignment
    ref_aligned = alignment_df[alignment_df['fasta_id'] == reference_id]['aligned_seq'].iloc[0]

    # For now, assume sequences are similar enough to map directly
    # A proper implementation would use pairwise alignment

    # Simple approach: insert gaps based on reference
    aligned_reported = []
    reported_idx = 0

    for aln_char in ref_aligned:
        if aln_char == '-':
            aligned_reported.append('-')
        else:
            if reported_idx < len(reported_seq):
                aligned_reported.append(reported_seq[reported_idx])
                reported_idx += 1
            else:
                aligned_reported.append('-')

    return ''.join(aligned_reported)


def analyze_reported_udh(enzyme_name, aligned_seq, pdb_to_aln, best_aas, substrate):
    """
    Analyze a reported UDH enzyme

    Returns:
    --------
    dict : Analysis results
    """
    results = {
        'enzyme': enzyme_name,
        'substrate': substrate,
        'total_positions': 0,
        'best_matches': 0,
        'worst_matches': 0,
        'neutral_matches': 0,
        'score': 0.0,
        'details': [],
        'mutations_needed': []
    }

    substrate_best = best_aas[substrate]

    for pdb_pos, info in substrate_best.items():
        if pdb_pos not in pdb_to_aln:
            continue

        aln_pos = pdb_to_aln[pdb_pos]
        if aln_pos >= len(aligned_seq):
            continue

        current_aa = aligned_seq[aln_pos]
        if current_aa == '-':
            continue

        results['total_positions'] += 1
        best_aa = info['best_aa']
        worst_aa = info['worst_aa']
        delta = info['delta']
        significant = info['significant']

        if current_aa == best_aa:
            results['best_matches'] += 1
            if significant and delta is not None:
                results['score'] += delta
        elif current_aa == worst_aa:
            results['worst_matches'] += 1
            if significant and delta is not None:
                results['score'] -= delta
        else:
            results['neutral_matches'] += 1

        results['details'].append({
            'position': pdb_pos,
            'current_aa': current_aa,
            'best_aa': best_aa,
            'worst_aa': worst_aa,
            'delta': delta,
            'significant': significant,
            'status': 'BEST' if current_aa == best_aa else ('WORST' if current_aa == worst_aa else 'NEUTRAL')
        })

        # If not best AA and significant, suggest mutation
        if current_aa != best_aa and significant and delta is not None and delta > 0.3:
            results['mutations_needed'].append({
                'position': pdb_pos,
                'current': current_aa,
                'target': best_aa,
                'expected_gain': delta
            })

    # Normalize score
    if results['total_positions'] > 0:
        results['normalized_score'] = results['score'] / results['total_positions']
        results['best_match_ratio'] = results['best_matches'] / results['total_positions']

    return results


def main():
    """Main analysis"""
    print("="*80)
    print("Reported UDH 효소 분석")
    print("="*80)

    # 1. Load reported UDH sequences
    print("\n1. Reported UDH 서열 로딩...")
    fasta_file = 'reported_UDHs_complete.fasta'
    reported_udhs = parse_fasta(fasta_file)
    print(f"  로딩된 효소: {len(reported_udhs)}개")
    for name in reported_udhs.keys():
        print(f"    - {name}")

    # 2. Load alignment to get PDB→Alignment mapping
    print("\n2. Alignment 로딩...")
    alignment_df = parse_clustal_alignment('DWC0THCP211-alignment.aln')
    reference_seq = alignment_df.iloc[0]['aligned_seq']

    # Create PDB to alignment mapping
    from sequence_comparison_analyzer import create_pdb_to_alignment_mapping
    pdb_to_aln = create_pdb_to_alignment_mapping(reference_seq)
    print(f"  PDB → Alignment 매핑 생성 완료")

    # 3. Load best amino acids from Excel
    print("\n3. Excel에서 최적 아미노산 정보 로딩...")
    excel_file = 'analysis_output_20251105_171401/active_site_residues_20251105_171401.xlsx'
    best_aas = load_best_amino_acids(excel_file)
    print(f"  로딩 완료: {len(best_aas)} 기질")

    # 4. Analyze each reported UDH for each substrate
    print("\n4. 각 효소 분석 중...")

    all_results = []

    for enzyme_name, seq in reported_udhs.items():
        # Map to alignment (simple version - assumes similar to reference)
        aligned_seq = map_reported_udh_to_alignment(seq, alignment_df)

        for substrate in best_aas.keys():
            result = analyze_reported_udh(enzyme_name, aligned_seq, pdb_to_aln, best_aas, substrate)
            all_results.append(result)

    # 5. Summarize results
    print("\n" + "="*80)
    print("분석 결과 요약")
    print("="*80)

    # Group by substrate
    substrates = ['NADplus', 'alpha_D_Glucuronate', 'alpha_L_Guluronate',
                  'beta_D_Galacturonate', 'beta_D_Mannuronate']

    for substrate in substrates:
        print(f"\n{'='*80}")
        print(f"{substrate} - 상위 5개 효소")
        print("-"*80)

        substrate_results = [r for r in all_results if r['substrate'] == substrate]
        substrate_results.sort(key=lambda x: x['score'], reverse=True)

        for i, result in enumerate(substrate_results[:5], 1):
            print(f"\n{i}. {result['enzyme']}")
            print(f"   Score: {result['score']:.3f} (normalized: {result['normalized_score']:.3f})")
            print(f"   Best AA 일치: {result['best_matches']}/{result['total_positions']} ({result['best_match_ratio']*100:.1f}%)")
            print(f"   Worst AA: {result['worst_matches']}")
            print(f"   Neutral: {result['neutral_matches']}")

            if result['mutations_needed']:
                print(f"\n   권장 Mutations ({len(result['mutations_needed'])}개):")
                for mut in result['mutations_needed'][:5]:  # Show top 5
                    print(f"     P{mut['position']}: {mut['current']} → {mut['target']} (gain={mut['expected_gain']:.3f})")

    # 6. Create detailed Excel report
    print("\n" + "="*80)
    print("상세 Excel 리포트 생성...")
    print("="*80)

    # Create summary DataFrame
    summary_data = []
    for result in all_results:
        summary_data.append({
            'Enzyme': result['enzyme'],
            'Substrate': result['substrate'],
            'Score': result['score'],
            'Normalized_Score': result['normalized_score'],
            'Total_Positions': result['total_positions'],
            'Best_Matches': result['best_matches'],
            'Best_Match_Ratio': result['best_match_ratio'],
            'Worst_Matches': result['worst_matches'],
            'Neutral_Matches': result['neutral_matches'],
            'Mutations_Needed': len(result['mutations_needed'])
        })

    summary_df = pd.DataFrame(summary_data)

    # Create detailed mutations DataFrame
    mutations_data = []
    for result in all_results:
        for mut in result['mutations_needed']:
            mutations_data.append({
                'Enzyme': result['enzyme'],
                'Substrate': result['substrate'],
                'Position': mut['position'],
                'Current_AA': mut['current'],
                'Target_AA': mut['target'],
                'Expected_Gain': mut['expected_gain']
            })

    mutations_df = pd.DataFrame(mutations_data)

    # Save to Excel
    output_file = 'reported_udh_analysis.xlsx'
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        mutations_df.to_excel(writer, sheet_name='Engineering_Suggestions', index=False)

        # Detailed results for each enzyme
        for enzyme_name in reported_udhs.keys():
            enzyme_results = [r for r in all_results if r['enzyme'] == enzyme_name]

            details_data = []
            for result in enzyme_results:
                for detail in result['details']:
                    details_data.append({
                        'Substrate': result['substrate'],
                        'Position': detail['position'],
                        'Current_AA': detail['current_aa'],
                        'Best_AA': detail['best_aa'],
                        'Worst_AA': detail['worst_aa'],
                        'Delta': detail['delta'],
                        'Significant': detail['significant'],
                        'Status': detail['status']
                    })

            if details_data:
                details_df = pd.DataFrame(details_data)
                sheet_name = enzyme_name[:31]  # Excel sheet name limit
                details_df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\nExcel 리포트 저장 완료: {output_file}")
    print("\n분석 완료!")


if __name__ == "__main__":
    main()
