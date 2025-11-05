"""
Demo script for running sequence comparison analysis
This shows what the analysis will tell you
"""

from sequence_comparison_analyzer import SequenceComparisonAnalyzer
import pandas as pd
import os


def main():
    print("="*80)
    print("SEQUENCE COMPARISON ANALYSIS - DEMO")
    print("="*80)
    print("\n이 분석은 다음을 알려줍니다:")
    print("1. 각 position에서 어떤 amino acid가 kinetic parameter를 높이는지")
    print("2. NAD vs 4개 uronate에서 선호하는 amino acid가 다른지")
    print("3. CataPro가 어떤 sequence를 높은 활성이라고 예측하는지")
    print("="*80)

    # Define binding sites (example - adjust based on your actual binding sites)
    binding_sites = {
        84: 'NAD_ribose',
        85: 'NAD_adenine',
        122: 'substrate_binding',
        155: 'catalytic_site',
        180: 'cofactor_binding'
    }

    window_size = 3

    # Initialize analyzer
    analyzer = SequenceComparisonAnalyzer(binding_sites, window_size=window_size)

    print(f"\n분석 설정:")
    print(f"  Core binding sites: {len(binding_sites)}")
    print(f"  Window size: ±{window_size}")
    print(f"  Total positions to analyze: {len(analyzer.positions_with_window)}")

    # Check for required files
    substrate_files = {
        'NADplus': 'catapro_predict_score_NADplus.csv',
        'alpha_D_Glucuronate': 'catapro_predict_score_α-D-Glucuronate.csv',
        'alpha_L_Guluronate': 'catapro_predict_score_α-L-Guluronate.csv',
        'beta_D_Galacturonate': 'catapro_predict_score_β-D-Galacturonate.csv',
        'beta_D_Mannuronate': 'catapro_predict_score_β-D-Mannuronate.csv'
    }

    print(f"\n필요한 파일 확인:")
    missing_files = []
    for substrate, filepath in substrate_files.items():
        exists = os.path.exists(filepath)
        status = "O" if exists else "X"
        print(f"  [{status}] {filepath}")
        if not exists:
            missing_files.append(filepath)

    if missing_files:
        print(f"\n경고: {len(missing_files)}개 파일이 없습니다.")
        print("\n다음 파일들이 필요합니다:")
        for f in missing_files:
            print(f"  - {f}")
        print("\n이 파일들은 다음을 포함해야 합니다:")
        print("  1. fasta_id 컬럼")
        print("  2. aligned_seq 컬럼 (MSA 정렬된 서열)")
        print("  3. pred_log10[kcat(s^-1)] 컬럼")
        print("  4. pred_log10[Km(mM)] 컬럼")
        print("  5. pred_log10[kcat/Km(s^-1mM^-1)] 컬럼")
        print("\n파일을 준비한 후 다시 실행해주세요.")
        return

    # Load and analyze
    print("\n" + "="*80)
    print("데이터 로딩 중...")
    print("="*80)

    try:
        # Load substrate data
        analyzer.load_substrate_data(substrate_files)

        # Extract position data
        print("\nBinding site position 데이터 추출 중...")
        analyzer.extract_position_data(aligned_seq_col='aligned_seq')

        # Analyze amino acid effects
        print("\nAmino acid 효과 분석 중...")
        aa_effects = analyzer.analyze_amino_acid_effects(min_sequences=5)

        print(f"\n분석 결과:")
        print(f"  총 분석: {len(aa_effects)} position-substrate-parameter 조합")
        print(f"  유의한 결과 (p<0.05): {aa_effects['significant'].sum()}")

        # Show examples
        print("\n" + "="*80)
        print("예제 결과: NADplus에 대한 log(kcat/Km) 상위 5 positions")
        print("="*80)

        nadplus_kcat_km = aa_effects[
            (aa_effects['substrate'] == 'NADplus') &
            (aa_effects['parameter'] == 'log_kcat_Km') &
            (aa_effects['significant'])
        ].nlargest(5, 'delta')

        if len(nadplus_kcat_km) > 0:
            for _, row in nadplus_kcat_km.iterrows():
                print(f"\nPosition {row['position']} ({row['label']}):")
                print(f"  최고 amino acid: {row['best_aa']} (평균 = {row['best_mean']:.3f})")
                print(f"  최저 amino acid: {row['worst_aa']} (평균 = {row['worst_mean']:.3f})")
                print(f"  차이 (delta): {row['delta']:.3f}")
                print(f"  P-value: {row['p_value']:.4f}")
                print(f"  해석: {row['best_aa']}를 가진 서열들이 {row['worst_aa']}보다")
                print(f"        log(kcat/Km)이 {row['delta']:.3f} 높음")

        # Find substrate-specific residues
        print("\n" + "="*80)
        print("기질 특이적 잔기 탐색 중...")
        print("="*80)

        substrate_specific = analyzer.find_substrate_specific_residues(
            parameter='log_kcat_Km',
            min_delta=0.1
        )

        if len(substrate_specific) > 0:
            print(f"\n발견: {len(substrate_specific)}개 positions이 기질에 따라 다른 amino acid 선호")
            print("\n상위 5개:")

            for _, row in substrate_specific.nlargest(5, 'max_delta').iterrows():
                print(f"\nPosition {row['position']} ({row['label']}):")
                print(f"  서로 다른 선호: {row['n_different_preferences']}개")
                print(f"  최대 효과 차이: {row['max_delta']:.3f}")
                print(f"  해석: 이 위치는 기질 특이성에 중요함")

        # Generate outputs
        print("\n" + "="*80)
        print("결과 파일 생성 중...")
        print("="*80)

        report_file = analyzer.generate_report()
        analyzer.plot_amino_acid_effects()

        print("\n✓ 분석 완료!")
        print(f"  - 리포트: {report_file}")
        print(f"  - 시각화: aa_effects_*.png, substrate_comparison_*.png")

        # Summary
        print("\n" + "="*80)
        print("결과 요약")
        print("="*80)

        for substrate in analyzer.substrates:
            sub_data = aa_effects[aa_effects['substrate'] == substrate]
            sig_data = sub_data[sub_data['significant']]

            print(f"\n{substrate}:")
            print(f"  분석된 positions: {len(sub_data['position'].unique())}")
            print(f"  유의한 효과: {len(sig_data)}")

            if len(sig_data) > 0:
                max_effect = sig_data.nlargest(1, 'delta').iloc[0]
                print(f"  최대 효과: Position {max_effect['position']} "
                      f"({max_effect['best_aa']} vs {max_effect['worst_aa']}, "
                      f"delta={max_effect['delta']:.3f})")

    except FileNotFoundError as e:
        print(f"\n오류: 파일을 찾을 수 없습니다: {e}")
    except KeyError as e:
        print(f"\n오류: 필요한 컬럼이 없습니다: {e}")
        print("파일에 'aligned_seq' 컬럼이 있는지 확인해주세요.")
    except Exception as e:
        print(f"\n예상치 못한 오류: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "="*80)


if __name__ == "__main__":
    main()
