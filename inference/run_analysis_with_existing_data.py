"""
Run sequence comparison analysis with existing data files
Uses: DWC0THCP211-alignment.aln + catapro_predict_score.csv
"""

import pandas as pd
import numpy as np
from parse_alignment_simple import parse_clustal_alignment
from sequence_comparison_analyzer import SequenceComparisonAnalyzer
import os


def identify_substrates_by_smiles(df):
    """
    Identify substrates from SMILES patterns
    """
    smiles_map = {}

    for idx, (smiles, count) in enumerate(df['smiles'].value_counts().items()):
        # NAD+ has much longer SMILES (contains adenine, ribose, phosphate)
        if 'N+' in smiles and len(smiles) > 100:
            smiles_map['NADplus'] = smiles
        elif idx == 0:  # First uronate
            smiles_map['uronate_1'] = smiles
        elif idx == 1:
            smiles_map['uronate_2'] = smiles
        elif idx == 2:
            smiles_map['uronate_3'] = smiles
        elif idx == 3:
            smiles_map['uronate_4'] = smiles

    # Based on common patterns in uronate SMILES
    # α/β determined by stereochemistry, D/L by configuration
    uronate_names = [
        'alpha_D_Glucuronate',
        'alpha_L_Guluronate',
        'beta_D_Galacturonate',
        'beta_D_Mannuronate'
    ]

    substrate_smiles = {}
    uronate_idx = 0

    for name, smiles in smiles_map.items():
        if name == 'NADplus':
            substrate_smiles['NADplus'] = smiles
        else:
            if uronate_idx < len(uronate_names):
                substrate_smiles[uronate_names[uronate_idx]] = smiles
                uronate_idx += 1

    return substrate_smiles


def match_alignment_to_catapro(alignment_df, catapro_df):
    """
    Match alignment sequences to CataPro by position
    Assumes both are in same order
    """
    print("\n서열 매칭 중...")

    # Both should have 1181 sequences
    n_unique_catapro = catapro_df['fasta_id'].nunique()
    if len(alignment_df) != n_unique_catapro:
        print(f"  경고: 서열 수 불일치!")
        print(f"    Alignment: {len(alignment_df)}")
        print(f"    CataPro unique IDs: {n_unique_catapro}")

    # Try matching by position (assume same order)
    # Each sequence in alignment should match 5 rows in CataPro (5 substrates)

    # Get first sequence per fasta_id in CataPro
    catapro_unique = catapro_df.drop_duplicates(subset=['fasta_id'], keep='first')

    if len(alignment_df) == len(catapro_unique):
        print(f"  [OK] 서열 수 일치: {len(alignment_df)}")

        # Add position index for matching
        alignment_df['seq_idx'] = range(len(alignment_df))
        catapro_unique['seq_idx'] = range(len(catapro_unique))

        # Map fasta_id
        id_map = dict(zip(catapro_unique['seq_idx'], catapro_unique['fasta_id']))
        alignment_df['catapro_id'] = alignment_df['seq_idx'].map(id_map)

        print(f"  매칭 완료!")
        print(f"    Alignment ID 예: {alignment_df['fasta_id'].iloc[0]}")
        print(f"    CataPro ID 예: {alignment_df['catapro_id'].iloc[0]}")

        return alignment_df
    else:
        print(f"  매칭 실패 - 수동 매칭 필요")
        return None


def create_substrate_dataframes(alignment_df, catapro_df, substrate_smiles):
    """
    Create separate DataFrames for each substrate with aligned sequences
    """
    print("\n기질별 데이터 생성 중...")

    substrate_dfs = {}

    for substrate, smiles in substrate_smiles.items():
        print(f"\n  {substrate}...")

        # Get CataPro data for this substrate
        sub_catapro = catapro_df[catapro_df['smiles'] == smiles].copy()

        # Merge with alignment using catapro_id
        merged = alignment_df.merge(
            sub_catapro,
            left_on='catapro_id',
            right_on='fasta_id',
            how='inner'
        )

        # Rename for consistency
        merged['fasta_id'] = merged['catapro_id']

        # Add log columns for analysis
        if 'pred_log10[kcat(s^-1)]' in merged.columns:
            merged['log_kcat'] = merged['pred_log10[kcat(s^-1)]']
        if 'pred_log10[Km(mM)]' in merged.columns:
            merged['log_Km'] = merged['pred_log10[Km(mM)]']
        if 'pred_log10[kcat/Km(s^-1mM^-1)]' in merged.columns:
            merged['log_kcat_Km'] = merged['pred_log10[kcat/Km(s^-1mM^-1)]']

        substrate_dfs[substrate] = merged
        print(f"    {len(merged)} 서열")

    return substrate_dfs


def main():
    """Main execution"""
    print("="*80)
    print("기존 데이터로 Sequence Comparison 분석")
    print("="*80)

    # 1. Load alignment
    print("\n1. Alignment 로딩...")
    alignment_df = parse_clustal_alignment('DWC0THCP211-alignment.aln')

    # 2. Load CataPro data
    print("\n2. CataPro 데이터 로딩...")
    catapro_df = pd.read_csv('catapro_predict_score.csv')
    if 'Unnamed: 0' in catapro_df.columns:
        catapro_df = catapro_df.drop('Unnamed: 0', axis=1)

    print(f"  총 {len(catapro_df)} 행")

    # 3. Identify substrates
    print("\n3. 기질 식별 (SMILES 기반)...")
    substrate_smiles = identify_substrates_by_smiles(catapro_df)

    for substrate, smiles in substrate_smiles.items():
        count = len(catapro_df[catapro_df['smiles'] == smiles])
        print(f"  {substrate}: {count} 서열")
        print(f"    SMILES: {smiles[:60]}...")

    # 4. Match alignment to CataPro
    print("\n4. Alignment와 CataPro 매칭...")
    alignment_matched = match_alignment_to_catapro(alignment_df, catapro_df)

    if alignment_matched is None:
        print("\n오류: 매칭 실패")
        return

    # 5. Create substrate DataFrames
    print("\n5. 기질별 데이터 준비...")
    substrate_dfs = create_substrate_dataframes(
        alignment_matched, catapro_df, substrate_smiles
    )

    # 6. Define binding sites
    print("\n6. Binding sites 정의...")
    print("\nPDB 3RFV (Uronate dehydrogenase)의 알려진 binding sites 사용:")

    # Known binding sites for UDH family (0-indexed for Python)
    # Based on conserved positions in alignment
    binding_sites = {
        11: 'NAD_binding_1',   # Gly-rich loop
        12: 'NAD_binding_2',
        13: 'NAD_binding_3',
        86: 'catalytic_1',     # Active site
        109: 'catalytic_2',
        135: 'substrate_1',    # Substrate binding
        158: 'substrate_2',
        180: 'cofactor',       # Cofactor interaction
    }

    print(f"  정의된 core sites: {len(binding_sites)}")
    for pos, name in sorted(binding_sites.items()):
        print(f"    Position {pos+1} (1-indexed): {name}")

    window_size = 3
    print(f"  Window size: ±{window_size}")

    # 7. Initialize analyzer
    print("\n7. Analyzer 초기화...")
    analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=window_size)

    total_positions = len(analyzer.positions_with_window)
    print(f"  총 분석 positions: {total_positions}")
    print(f"  Core sites: {len(binding_sites)}")

    # 8. Prepare data in memory (don't write files)
    print("\n8. 데이터 통합...")

    # Manually set substrate data
    analyzer.substrate_data = substrate_dfs
    analyzer.substrates = list(substrate_dfs.keys())

    # Extract position data
    print("\n9. Position 데이터 추출...")
    analyzer.extract_position_data(aligned_seq_col='aligned_seq')

    # 10. Analyze amino acid effects
    print("\n10. 아미노산 효과 분석 중...")
    print("    (이 과정은 몇 분 소요될 수 있습니다...)")

    aa_effects = analyzer.analyze_amino_acid_effects(
        kinetic_params=['log_kcat', 'log_Km', 'log_kcat_Km'],
        min_sequences=5
    )

    print(f"\n분석 완료!")
    print(f"  총 분석: {len(aa_effects)} position-substrate-parameter 조합")
    print(f"  유의한 결과 (p<0.05): {aa_effects['significant'].sum()}")

    # 11. Find substrate-specific residues
    print("\n11. 기질 특이적 잔기 탐색...")
    substrate_specific = analyzer.find_substrate_specific_residues(
        parameter='log_kcat_Km',
        min_delta=0.1
    )

    # 12. Show results
    print("\n" + "="*80)
    print("주요 결과")
    print("="*80)

    # Top effects for NADplus
    if 'NADplus' in aa_effects['substrate'].values:
        print("\nNADplus에 대한 log(kcat/Km) 상위 10 positions:")
        print("-"*80)

        nadplus_kcat_km = aa_effects[
            (aa_effects['substrate'] == 'NADplus') &
            (aa_effects['parameter'] == 'log_kcat_Km') &
            (aa_effects['significant'])
        ].nlargest(10, 'delta')

        if len(nadplus_kcat_km) > 0:
            for i, (_, row) in enumerate(nadplus_kcat_km.iterrows(), 1):
                print(f"\n{i}. Position {row['position']} ({row['label']}):")
                print(f"   최고: {row['best_aa']} (평균={row['best_mean']:.3f}, n={row['best_count']}, unique_seqs={row['best_unique_seqs']})")
                print(f"        서열 다양성: {row['best_seq_identity_mean']:.1f}% identity (range: {row['best_seq_identity_min']:.1f}-{row['best_seq_identity_max']:.1f}%)")
                print(f"   최저: {row['worst_aa']} (평균={row['worst_mean']:.3f}, n={row['worst_count']}, unique_seqs={row['worst_unique_seqs']})")
                print(f"        서열 다양성: {row['worst_seq_identity_mean']:.1f}% identity (range: {row['worst_seq_identity_min']:.1f}-{row['worst_seq_identity_max']:.1f}%)")
                print(f"   차이: {row['delta']:.3f} (p={row['p_value']:.4f})")
        else:
            print("  유의한 결과 없음")

    # Substrate-specific positions
    if len(substrate_specific) > 0:
        print("\n" + "="*80)
        print(f"기질 특이적 잔기: {len(substrate_specific)}개 발견")
        print("="*80)

        for i, (_, row) in enumerate(substrate_specific.nlargest(5, 'max_delta').iterrows(), 1):
            print(f"\n{i}. Position {row['position']} ({row['label']}):")
            print(f"   서로 다른 선호: {row['n_different_preferences']}개")
            print(f"   최대 효과 차이: {row['max_delta']:.3f}")

    # 13. Generate outputs
    print("\n" + "="*80)
    print("결과 저장 중...")
    print("="*80)

    report_file = analyzer.generate_report()
    analyzer.plot_amino_acid_effects()

    # 14. Plot sample sizes by position
    print("\n" + "="*80)
    print("Position별 아미노산 분포 시각화...")
    print("="*80)
    sample_size_summary = analyzer.plot_sample_sizes_by_position()
    print("\nSample size summary:")
    print(sample_size_summary.to_string())

    # Summary by substrate
    print("\n" + "="*80)
    print("기질별 요약")
    print("="*80)

    for substrate in analyzer.substrates:
        sub_data = aa_effects[
            (aa_effects['substrate'] == substrate) &
            (aa_effects['parameter'] == 'log_kcat_Km')
        ]
        sig_data = sub_data[sub_data['significant']]

        print(f"\n{substrate}:")
        print(f"  분석된 positions: {len(sub_data['position'].unique())}")
        print(f"  유의한 효과: {len(sig_data)}")

        if len(sig_data) > 0:
            max_effect = sig_data.nlargest(1, 'delta').iloc[0]
            print(f"  최대 효과: Position {max_effect['position']} "
                  f"({max_effect['best_aa']} vs {max_effect['worst_aa']}, "
                  f"delta={max_effect['delta']:.3f})")

    print("\n" + "="*80)
    print("분석 완료!")
    print("="*80)
    print(f"\n생성된 파일:")
    print(f"  - {report_file}")
    print(f"  - aa_effects_*.png (parameter별)")
    print(f"  - substrate_comparison_top_positions_*.png")
    print(f"  - sample_sizes_by_position_*.png (각 position별 n수 분포)")
    print(f"  - sample_size_summary_*.csv (n수 요약 테이블)")


if __name__ == "__main__":
    main()
