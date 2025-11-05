"""
아미노산 분석에서 제안된 mutation과 CataPro 예측 성능 비교
Analyze correlation between mutation suggestions and CataPro predictions
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_data():
    """Load amino acid analysis and CataPro predictions"""

    # Load amino acid analysis results
    summary_df = pd.read_excel('reported_udh_analysis.xlsx', sheet_name='Summary')
    engineering_df = pd.read_excel('reported_udh_analysis.xlsx', sheet_name='Engineering_Suggestions')

    # Load CataPro predictions summary
    catapro_df = pd.read_csv('reported_udh_rankings_summary.csv')

    return summary_df, engineering_df, catapro_df

def analyze_mutation_correlation():
    """Analyze correlation between mutations needed and CataPro performance"""

    print("="*80)
    print("Mutation Suggestions from AA Analysis vs CataPro Predictions")
    print("="*80)

    summary_df, engineering_df, catapro_df = load_data()

    # Substrate list
    substrates = ['NADplus', 'alpha_D_Glucuronate', 'alpha_L_Guluronate',
                  'beta_D_Galacturonate', 'beta_D_Mannuronate']

    substrate_names_kr = {
        'NADplus': 'NAD+',
        'alpha_D_Glucuronate': 'α-D-글루쿠론산',
        'alpha_L_Guluronate': 'α-L-굴루론산',
        'beta_D_Galacturonate': 'β-D-갈락투론산',
        'beta_D_Mannuronate': 'β-D-만누론산'
    }

    # Get unique enzymes
    enzymes = summary_df['Enzyme'].unique()

    print(f"\nAnalysis target: {len(enzymes)}enzymes")
    print(f"substrates: {len(substrates)}개")

    # === 1. Enzyme-level comparison ===
    print("\n" + "="*80)
    print("1. Enzyme Comparison: AA Analysis Score vs CataPro Predictions")
    print("="*80)

    comparison_data = []

    for enzyme in enzymes:
        # Get amino acid analysis scores for each substrate
        enzyme_summary = summary_df[summary_df['Enzyme'] == enzyme]

        # Get CataPro performance
        catapro_row = catapro_df[catapro_df['Enzyme'] == enzyme]

        if len(catapro_row) == 0:
            continue

        catapro_score = catapro_row['Mean_kcat_over_Km'].iloc[0]

        for substrate in substrates:
            substrate_row = enzyme_summary[enzyme_summary['Substrate'] == substrate]

            if len(substrate_row) == 0:
                continue

            aa_score = substrate_row['Score'].iloc[0]
            normalized_score = substrate_row['Normalized_Score'].iloc[0]
            best_matches = substrate_row['Best_Matches'].iloc[0]
            mutations_needed = substrate_row['Mutations_Needed'].iloc[0]

            comparison_data.append({
                'Enzyme': enzyme,
                'Substrate': substrate,
                'AA_Score': aa_score,
                'AA_Normalized_Score': normalized_score,
                'Best_Match_Count': best_matches,
                'Mutations_Needed': mutations_needed,
                'CataPro_kcat_over_Km': catapro_score
            })

    comparison_df = pd.DataFrame(comparison_data)

    # Calculate average scores per enzyme
    enzyme_avg = comparison_df.groupby('Enzyme').agg({
        'AA_Score': 'mean',
        'AA_Normalized_Score': 'mean',
        'Best_Match_Count': 'mean',
        'Mutations_Needed': 'mean',
        'CataPro_kcat_over_Km': 'first'
    }).reset_index()

    # Sort by CataPro score
    enzyme_avg = enzyme_avg.sort_values('CataPro_kcat_over_Km', ascending=False)
    enzyme_avg['CataPro_Rank'] = range(1, len(enzyme_avg) + 1)

    # Sort by AA score
    enzyme_avg_by_aa = enzyme_avg.sort_values('AA_Score', ascending=False).copy()
    enzyme_avg_by_aa['AA_Rank'] = range(1, len(enzyme_avg_by_aa) + 1)

    # Merge ranks
    enzyme_avg = enzyme_avg.merge(
        enzyme_avg_by_aa[['Enzyme', 'AA_Rank']],
        on='Enzyme'
    )

    print("\nAverage scores and rankings per enzyme:")
    print("-"*120)
    print(f"{'Enzyme':<10} {'AA_Score':>8} {'AA_Rank':>6} {'Optimal_AA':>7} {'Mutations':>8} {'CataPro':>10} {'CataPro_Rank':>11} {'Rank_Diff':>8}")
    print("-"*120)

    for _, row in enzyme_avg.iterrows():
        rank_diff = int(row['AA_Rank']) - int(row['CataPro_Rank'])
        print(f"{row['Enzyme']:<10} {row['AA_Score']:>8.2f} {int(row['AA_Rank']):>6} "
              f"{row['Best_Match_Count']:>7.1f} {row['Mutations_Needed']:>8.1f} "
              f"{row['CataPro_kcat_over_Km']:>10.3f} {int(row['CataPro_Rank']):>11} "
              f"{rank_diff:>+8}")

    # === 2. Correlation analysis ===
    print("\n" + "="*80)
    print("2. Correlation Analysis")
    print("="*80)

    from scipy.stats import pearsonr, spearmanr

    # AA Score vs CataPro
    corr_aa, p_aa = pearsonr(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'])
    scorr_aa, sp_aa = spearmanr(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'])

    print(f"\nAA Score vs CataPro kcat/Km:")
    print(f"  Pearson correlation: {corr_aa:.3f} (p={p_aa:.4f})")
    print(f"  Spearman rank correlation: {scorr_aa:.3f} (p={sp_aa:.4f})")

    # Best matches vs CataPro
    corr_best, p_best = pearsonr(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'])
    scorr_best, sp_best = spearmanr(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'])

    print(f"\nOptimal AA matches vs CataPro kcat/Km:")
    print(f"  Pearson correlation: {corr_best:.3f} (p={p_best:.4f})")
    print(f"  Spearman rank correlation: {scorr_best:.3f} (p={sp_best:.4f})")

    # Mutations needed vs CataPro (negative correlation expected)
    corr_mut, p_mut = pearsonr(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'])
    scorr_mut, sp_mut = spearmanr(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'])

    print(f"\nMutations needed vs CataPro kcat/Km:")
    print(f"  Pearson correlation: {corr_mut:.3f} (p={p_mut:.4f})")
    print(f"  Spearman rank correlation: {scorr_mut:.3f} (p={sp_mut:.4f})")

    # === 3. Top performers analysis ===
    print("\n" + "="*80)
    print("3. 상위 Enzyme 비교 분석")
    print("="*80)

    print("\n[CataPro 상위 5enzymes]")
    print("-"*80)
    for i, (_, row) in enumerate(enzyme_avg.head(5).iterrows(), 1):
        print(f"\n{i}. {row['Enzyme']}")
        print(f"   CataPro kcat/Km: {row['CataPro_kcat_over_Km']:.3f} (순위 #{int(row['CataPro_Rank'])})")
        print(f"   AA 분석 점수: {row['AA_Score']:.2f} (순위 #{int(row['AA_Rank'])})")
        print(f"   평균 최적 AA: {row['Best_Match_Count']:.1f}/44개")
        print(f"   평균 필요 변이: {row['Mutations_Needed']:.1f}개")

        # Get mutation details
        enzyme_mutations = engineering_df[engineering_df['Enzyme'] == row['Enzyme']]
        if len(enzyme_mutations) > 0:
            print(f"   총 제안 변이: {len(enzyme_mutations)}개")
            top_muts = enzyme_mutations.nlargest(3, 'Expected_Gain')
            if len(top_muts) > 0:
                print(f"   상위 3개 변이:")
                for _, mut in top_muts.iterrows():
                    print(f"     - P{mut['Position']}: {mut['Current_AA']} → {mut['Target_AA']} "
                          f"(substrates: {substrate_names_kr[mut['Substrate']]}, 기대효과: {mut['Expected_Gain']:.3f})")

    print("\n[AA 분석 상위 5enzymes]")
    print("-"*80)
    for i, (_, row) in enumerate(enzyme_avg_by_aa.head(5).iterrows(), 1):
        print(f"\n{i}. {row['Enzyme']}")
        print(f"   AA 분석 점수: {row['AA_Score']:.2f} (순위 #{int(row['AA_Rank'])})")
        print(f"   CataPro kcat/Km: {row['CataPro_kcat_over_Km']:.3f} (순위 #{int(row['CataPro_Rank'])})")
        print(f"   평균 최적 AA: {row['Best_Match_Count']:.1f}/44개")
        print(f"   평균 필요 변이: {row['Mutations_Needed']:.1f}개")

    # === 4. Visualizations ===
    print("\n" + "="*80)
    print("4. Generating visualizations...")
    print("="*80)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: AA Score vs CataPro
    ax = axes[0, 0]
    scatter = ax.scatter(enzyme_avg['AA_Score'], enzyme_avg['CataPro_kcat_over_Km'],
                        c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                        s=100, alpha=0.7, edgecolors='black')
    for _, row in enzyme_avg.iterrows():
        ax.annotate(row['Enzyme'], (row['AA_Score'], row['CataPro_kcat_over_Km']),
                   xytext=(5, 5), textcoords='offset points', fontsize=8)
    ax.set_xlabel('AA Analysis Score', fontsize=11, fontweight='bold')
    ax.set_ylabel('CataPro Predicted kcat/Km', fontsize=11, fontweight='bold')
    ax.set_title(f'AA Score vs CataPro Performance\n(correlation: {corr_aa:.3f})',
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='CataPro Rank')

    # Plot 2: Best matches vs CataPro
    ax = axes[0, 1]
    scatter = ax.scatter(enzyme_avg['Best_Match_Count'], enzyme_avg['CataPro_kcat_over_Km'],
                        c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                        s=100, alpha=0.7, edgecolors='black')
    for _, row in enzyme_avg.iterrows():
        ax.annotate(row['Enzyme'], (row['Best_Match_Count'], row['CataPro_kcat_over_Km']),
                   xytext=(5, 5), textcoords='offset points', fontsize=8)
    ax.set_xlabel('최적 아미노산 매치 수', fontsize=11, fontweight='bold')
    ax.set_ylabel('CataPro 예측 kcat/Km', fontsize=11, fontweight='bold')
    ax.set_title(f'최적 AA vs CataPro 성능\n(상관계수: {corr_best:.3f})',
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='CataPro 순위')

    # Plot 3: Mutations needed vs CataPro
    ax = axes[1, 0]
    scatter = ax.scatter(enzyme_avg['Mutations_Needed'], enzyme_avg['CataPro_kcat_over_Km'],
                        c=enzyme_avg['CataPro_Rank'], cmap='RdYlGn_r',
                        s=100, alpha=0.7, edgecolors='black')
    for _, row in enzyme_avg.iterrows():
        ax.annotate(row['Enzyme'], (row['Mutations_Needed'], row['CataPro_kcat_over_Km']),
                   xytext=(5, 5), textcoords='offset points', fontsize=8)
    ax.set_xlabel('필요 변이 수', fontsize=11, fontweight='bold')
    ax.set_ylabel('CataPro 예측 kcat/Km', fontsize=11, fontweight='bold')
    ax.set_title(f'필요 변이 vs CataPro 성능\n(상관계수: {corr_mut:.3f})',
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='CataPro 순위')

    # Plot 4: Rank comparison
    ax = axes[1, 1]
    x = np.arange(len(enzyme_avg))
    width = 0.35

    bars1 = ax.bar(x - width/2, enzyme_avg['AA_Rank'], width, label='AA 분석 순위', alpha=0.7)
    bars2 = ax.bar(x + width/2, enzyme_avg['CataPro_Rank'], width, label='CataPro 순위', alpha=0.7)

    ax.set_xlabel('Enzyme', fontsize=11, fontweight='bold')
    ax.set_ylabel('순위 (낮을수록 좋음)', fontsize=11, fontweight='bold')
    ax.set_title('Enzyme별 순위 비교', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(enzyme_avg['Enzyme'], rotation=45, ha='right', fontsize=9)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig('mutation_catapro_correlation.png', dpi=300, bbox_inches='tight')
    print("Saved: mutation_catapro_correlation.png")
    plt.close()

    # === 5. Save detailed comparison ===
    enzyme_avg.to_csv('enzyme_comparison_aa_vs_catapro.csv', index=False)
    print("Saved: enzyme_comparison_aa_vs_catapro.csv")

    # === 6. Key insights ===
    print("\n" + "="*80)
    print("5. Key Findings")
    print("="*80)

    print(f"\n✓ 상관관계:")
    if abs(corr_aa) < 0.3:
        strength = "매우 약한"
    elif abs(corr_aa) < 0.5:
        strength = "약한"
    elif abs(corr_aa) < 0.7:
        strength = "중간"
    else:
        strength = "강한"

    print(f"  - AA 점수와 CataPro 성능: {strength} 상관관계 (r={corr_aa:.3f})")

    if p_aa < 0.05:
        print(f"  - 통계적으로 유의미함 (p={p_aa:.4f})")
    else:
        print(f"  - 통계적으로 유의미하지 않음 (p={p_aa:.4f})")

    print(f"\n✓ 순위 불일치:")
    rank_differences = abs(enzyme_avg['AA_Rank'] - enzyme_avg['CataPro_Rank'])
    avg_rank_diff = rank_differences.mean()
    max_rank_diff = rank_differences.max()

    print(f"  - 평균 순위 차이: {avg_rank_diff:.1f}")
    print(f"  - 최대 순위 차이: {int(max_rank_diff)}")

    worst_mismatch = enzyme_avg.loc[rank_differences.idxmax()]
    print(f"  - 가장 큰 불일치: {worst_mismatch['Enzyme']} "
          f"(AA #{int(worst_mismatch['AA_Rank'])} vs CataPro #{int(worst_mismatch['CataPro_Rank'])})")

    print(f"\n✓ 변이 제안의 효과 예측:")
    if abs(corr_mut) < 0.3:
        print(f"  - 필요 변이 수와 CataPro 성능의 상관관계가 약함 (r={corr_mut:.3f})")
        print(f"  - 즉, active site의 최적 아미노산만으로는 전체 Enzyme 성능을 설명하기 어려움")
    else:
        print(f"  - 필요 변이 수가 CataPro 성능과 관련있음 (r={corr_mut:.3f})")

    print("\n" + "="*80)
    print("분석 완료!")
    print("="*80)

    return enzyme_avg, comparison_df

if __name__ == "__main__":
    enzyme_avg, comparison_df = analyze_mutation_correlation()
