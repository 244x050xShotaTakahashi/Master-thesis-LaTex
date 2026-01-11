#!/usr/bin/env python3
"""
統計処理スクリプト: 複数シードの結果から平均・標準偏差を計算

使用方法:
    python scripts/summarize_charge_dist_aor_statistics.py

入力:
    results/charge_distribution_study/summary/aor_vs_cutoff.csv

出力:
    results/charge_distribution_study/summary/aor_statistics.csv
    results/charge_distribution_study/summary/aor_mean_matrix.csv
    results/charge_distribution_study/summary/aor_std_matrix.csv
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse


def compute_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """分布×カットオフごとに統計を計算"""
    
    # 完了したケースのみを使用
    df_completed = df[df['status'] == 'completed'].copy()
    
    if df_completed.empty:
        print("Error: No completed cases found.")
        return pd.DataFrame()
    
    # 分布×カットオフでグループ化して統計を計算
    stats = df_completed.groupby(['distribution', 'cutoff_m']).agg({
        'angle_deg': ['mean', 'std', 'min', 'max', 'count'],
        'r_squared': ['mean', 'std'],
        'elapsed_seconds': ['mean', 'sum'],
    }).reset_index()
    
    # カラム名をフラット化
    stats.columns = [
        'distribution', 'cutoff_m',
        'angle_mean', 'angle_std', 'angle_min', 'angle_max', 'n_samples',
        'r_squared_mean', 'r_squared_std',
        'elapsed_mean', 'elapsed_total',
    ]
    
    # 標準誤差を計算
    stats['angle_se'] = stats['angle_std'] / np.sqrt(stats['n_samples'])
    
    # 95%信頼区間を計算 (t分布の近似として1.96を使用)
    stats['angle_ci95_lower'] = stats['angle_mean'] - 1.96 * stats['angle_se']
    stats['angle_ci95_upper'] = stats['angle_mean'] + 1.96 * stats['angle_se']
    
    return stats


def create_mean_matrix(stats: pd.DataFrame) -> pd.DataFrame:
    """平均値のマトリックス形式"""
    if stats.empty:
        return pd.DataFrame()
    
    pivot = stats.pivot_table(
        index='cutoff_m',
        columns='distribution',
        values='angle_mean',
        aggfunc='first'
    )
    pivot = pivot.sort_index(ascending=False)
    return pivot


def create_std_matrix(stats: pd.DataFrame) -> pd.DataFrame:
    """標準偏差のマトリックス形式"""
    if stats.empty:
        return pd.DataFrame()
    
    pivot = stats.pivot_table(
        index='cutoff_m',
        columns='distribution',
        values='angle_std',
        aggfunc='first'
    )
    pivot = pivot.sort_index(ascending=False)
    return pivot


def create_formatted_matrix(stats: pd.DataFrame) -> pd.DataFrame:
    """平均±標準偏差の形式でマトリックスを作成"""
    if stats.empty:
        return pd.DataFrame()
    
    # 平均と標準偏差を組み合わせた文字列を作成
    stats = stats.copy()
    stats['formatted'] = stats.apply(
        lambda row: f"{row['angle_mean']:.2f}±{row['angle_std']:.2f}" 
                    if pd.notna(row['angle_std']) else f"{row['angle_mean']:.2f}",
        axis=1
    )
    
    pivot = stats.pivot_table(
        index='cutoff_m',
        columns='distribution',
        values='formatted',
        aggfunc='first'
    )
    pivot = pivot.sort_index(ascending=False)
    return pivot


def main():
    parser = argparse.ArgumentParser(
        description='複数シードの安息角結果から統計を計算'
    )
    parser.add_argument(
        '--input-csv',
        type=str,
        default='results/charge_distribution_study/summary/aor_vs_cutoff.csv',
        help='入力CSVファイル（aor_vs_cutoff.csv）'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='results/charge_distribution_study/summary',
        help='出力ディレクトリ'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("複数シード安息角統計処理")
    print("=" * 60)
    print(f"入力ファイル: {args.input_csv}")
    print(f"出力ディレクトリ: {args.output_dir}")
    print()
    
    # 入力CSVを読み込み
    if not os.path.exists(args.input_csv):
        print(f"Error: Input file not found: {args.input_csv}")
        print("Please run summarize_charge_dist_aor.py first.")
        sys.exit(1)
    
    df = pd.read_csv(args.input_csv)
    
    if df.empty:
        print("Error: Input CSV is empty.")
        sys.exit(1)
    
    # シード数を確認
    seed_ids = df['seed_id'].unique()
    n_seeds = len(seed_ids)
    
    print(f"検出されたシード数: {n_seeds}")
    print(f"シードID: {', '.join(sorted(seed_ids))}")
    print()
    
    # 統計を計算
    stats = compute_statistics(df)
    
    if stats.empty:
        print("Error: Could not compute statistics.")
        sys.exit(1)
    
    # 出力ディレクトリを作成
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 統計結果を保存
    stats_csv = os.path.join(args.output_dir, 'aor_statistics.csv')
    stats.to_csv(stats_csv, index=False)
    print(f"統計結果を保存: {stats_csv}")
    
    # 平均マトリックスを保存
    mean_matrix = create_mean_matrix(stats)
    if not mean_matrix.empty:
        mean_csv = os.path.join(args.output_dir, 'aor_mean_matrix.csv')
        mean_matrix.to_csv(mean_csv)
        print(f"平均マトリックスを保存: {mean_csv}")
    
    # 標準偏差マトリックスを保存
    std_matrix = create_std_matrix(stats)
    if not std_matrix.empty:
        std_csv = os.path.join(args.output_dir, 'aor_std_matrix.csv')
        std_matrix.to_csv(std_csv)
        print(f"標準偏差マトリックスを保存: {std_csv}")
    
    # フォーマット済みマトリックス（平均±標準偏差）を保存
    formatted_matrix = create_formatted_matrix(stats)
    if not formatted_matrix.empty:
        formatted_csv = os.path.join(args.output_dir, 'aor_formatted_matrix.csv')
        formatted_matrix.to_csv(formatted_csv)
        print(f"フォーマット済みマトリックスを保存: {formatted_csv}")
    
    print()
    print("-" * 60)
    print("統計サマリー")
    print("-" * 60)
    print()
    
    # サマリーを表示
    print(f"サンプル数: {stats['n_samples'].iloc[0]} シード/条件")
    print()
    
    if not mean_matrix.empty:
        print("安息角 平均 [deg]:")
        print(mean_matrix.round(2).to_string())
        print()
    
    if not std_matrix.empty:
        print("安息角 標準偏差 [deg]:")
        print(std_matrix.round(3).to_string())
        print()
    
    if not formatted_matrix.empty:
        print("安息角 平均±標準偏差 [deg]:")
        print(formatted_matrix.to_string())
        print()
    
    # 分布間の比較
    print("-" * 60)
    print("分布間比較（全カットオフ平均）")
    print("-" * 60)
    
    dist_summary = stats.groupby('distribution').agg({
        'angle_mean': 'mean',
        'angle_std': 'mean',
        'n_samples': 'sum',
    }).round(3)
    print(dist_summary.to_string())
    print()
    
    # カットオフ間の比較
    print("-" * 60)
    print("カットオフ間比較（全分布平均）")
    print("-" * 60)
    
    cutoff_summary = stats.groupby('cutoff_m').agg({
        'angle_mean': 'mean',
        'angle_std': 'mean',
        'n_samples': 'sum',
    }).sort_index(ascending=False).round(3)
    print(cutoff_summary.to_string())
    print()
    
    print("=" * 60)
    print("統計処理完了")
    print("=" * 60)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

