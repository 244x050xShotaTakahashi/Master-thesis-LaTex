#!/usr/bin/env python3
"""
結果集約スクリプト: 帯電分布×カットオフ距離の安息角データを集約

使用方法:
    python scripts/summarize_charge_dist_aor.py

出力:
    results/charge_distribution_study/summary/aor_vs_cutoff.csv
    results/charge_distribution_study/summary/aor_matrix.csv

対応ディレクトリ構造:
    - シード1（既存）: bimodal_rc_0.02/, normal_rc_0.02/, ...
    - シード2-5（新規）: bimodal_seed2_rc_0.02/, normal_seed3_rc_0.015/, ...
"""

import os
import sys
import glob
import re
import pandas as pd
import argparse
from pathlib import Path


def parse_repose_angle_csv(filepath: str) -> dict:
    """repose_angle_results.csv から安息角データを読み取る
    
    ファイル形式: parameter,value の縦形式
    """
    try:
        df = pd.read_csv(filepath)
        if df.empty:
            return None
        
        # parameter,value 形式を辞書に変換
        data = dict(zip(df['parameter'], df['value']))
        
        # 安息角データを抽出
        angle_deg = None
        r_squared = None
        
        # right_angle_deg または average_angle_deg を使用
        if 'right_angle_deg' in data:
            try:
                angle_deg = float(data['right_angle_deg'])
            except (ValueError, TypeError):
                pass
        elif 'average_angle_deg' in data:
            try:
                angle_deg = float(data['average_angle_deg'])
            except (ValueError, TypeError):
                pass
        
        if 'right_r_squared' in data:
            try:
                r_squared = float(data['right_r_squared'])
            except (ValueError, TypeError):
                pass
        
        return {
            'angle_deg': angle_deg,
            'r_squared': r_squared,
            'num_particles': data.get('num_particles', None),
            'simulation_time': data.get('simulation_time', None),
        }
    except Exception as e:
        print(f"Warning: Failed to parse {filepath}: {e}")
        return None


def parse_timing_csv(filepath: str) -> dict:
    """timing.csv から計算時間データを読み取る"""
    try:
        df = pd.read_csv(filepath)
        if df.empty:
            return None
        
        row = df.iloc[0]
        return {
            'elapsed_seconds': row.get('elapsed_seconds', None),
        }
    except Exception as e:
        print(f"Warning: Failed to parse {filepath}: {e}")
        return None


def parse_case_name(case_name: str) -> dict:
    """ケース名から分布、シードID、カットオフを抽出
    
    対応パターン:
    - bimodal_rc_0.02 → distribution=bimodal, seed_id=seed1, cutoff=0.02
    - bimodal_seed2_rc_0.02 → distribution=bimodal, seed_id=seed2, cutoff=0.02
    """
    # パターン1: {dist}_seed{N}_rc_{cutoff}
    match = re.match(r'^(bimodal|normal|uniform)_seed(\d+)_rc_(.+)$', case_name)
    if match:
        return {
            'distribution': match.group(1),
            'seed_id': f'seed{match.group(2)}',
            'cutoff_m': float(match.group(3)),
        }
    
    # パターン2: {dist}_rc_{cutoff} (シード1、既存形式)
    match = re.match(r'^(bimodal|normal|uniform)_rc_(.+)$', case_name)
    if match:
        return {
            'distribution': match.group(1),
            'seed_id': 'seed1',
            'cutoff_m': float(match.group(2)),
        }
    
    return None


def collect_phase2_results(base_dir: str) -> pd.DataFrame:
    """Phase 2 の全結果を収集（複数シード対応）"""
    results = []
    
    # ディレクトリを走査
    if not os.path.isdir(base_dir):
        print(f"Error: Directory not found: {base_dir}")
        return pd.DataFrame()
    
    for case_name in os.listdir(base_dir):
        case_dir = os.path.join(base_dir, case_name)
        
        if not os.path.isdir(case_dir):
            continue
        
        # ケース名をパース
        parsed = parse_case_name(case_name)
        if parsed is None:
            print(f"Warning: Could not parse case name: {case_name}")
            continue
        
        result = {
            'distribution': parsed['distribution'],
            'seed_id': parsed['seed_id'],
            'cutoff_m': parsed['cutoff_m'],
            'case_name': case_name,
            'status': 'not_found',
        }
        
        # repose_angle_results.csv を探す
        repose_file = os.path.join(case_dir, 'repose_angle_results.csv')
        if os.path.exists(repose_file):
            repose_data = parse_repose_angle_csv(repose_file)
            if repose_data:
                result.update(repose_data)
                result['status'] = 'completed'
        
        # timing.csv を探す
        timing_file = os.path.join(case_dir, 'timing.csv')
        if os.path.exists(timing_file):
            timing_data = parse_timing_csv(timing_file)
            if timing_data:
                result.update(timing_data)
        
        results.append(result)
    
    # DataFrameを作成してソート
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(['distribution', 'seed_id', 'cutoff_m']).reset_index(drop=True)
    
    return df


def create_matrix_view(df: pd.DataFrame, seed_id: str = None) -> pd.DataFrame:
    """分布×カットオフのマトリックス形式に変換
    
    seed_id が指定された場合はそのシードのみ、None の場合は全シードの平均
    """
    if df.empty:
        return pd.DataFrame()
    
    # angle_deg が存在し、有効な値があるか確認
    if 'angle_deg' not in df.columns:
        return pd.DataFrame()
    
    # シードでフィルタリング
    if seed_id is not None:
        df_filtered = df[df['seed_id'] == seed_id]
    else:
        df_filtered = df
    
    # angle_deg が全て None/NaN の場合は空のDataFrameを返す
    valid_angles = df_filtered['angle_deg'].dropna()
    if valid_angles.empty:
        return pd.DataFrame()
    
    # angle_deg をピボットテーブルに
    if seed_id is not None:
        # 特定シードの場合は first
        pivot = df_filtered.pivot_table(
            index='cutoff_m',
            columns='distribution',
            values='angle_deg',
            aggfunc='first'
        )
    else:
        # 全シードの場合は mean
        pivot = df_filtered.pivot_table(
            index='cutoff_m',
            columns='distribution',
            values='angle_deg',
            aggfunc='mean'
        )
    
    # カットオフ距離でソート（降順: 大きい方が上）
    pivot = pivot.sort_index(ascending=False)
    
    return pivot


def main():
    parser = argparse.ArgumentParser(
        description='帯電分布×カットオフ距離の安息角結果を集約（複数シード対応）'
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default='results/charge_distribution_study/phase2_aor',
        help='Phase 2 結果のベースディレクトリ'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='results/charge_distribution_study/summary',
        help='出力ディレクトリ'
    )
    
    args = parser.parse_args()
    
    # 出力ディレクトリを作成
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("帯電分布×カットオフ距離 安息角結果集約（複数シード対応）")
    print("=" * 60)
    print(f"入力ディレクトリ: {args.base_dir}")
    print(f"出力ディレクトリ: {args.output_dir}")
    print()
    
    # 結果を収集
    df = collect_phase2_results(args.base_dir)
    
    if df.empty:
        print("Error: No results found.")
        sys.exit(1)
    
    # 完了ケースの数を表示
    completed = df[df['status'] == 'completed']
    total = len(df)
    
    # シードごとの統計
    seed_counts = df.groupby('seed_id').size()
    seed_completed = df[df['status'] == 'completed'].groupby('seed_id').size()
    
    print(f"収集結果: {len(completed)}/{total} ケース完了")
    print()
    print("シード別:")
    for seed_id in sorted(df['seed_id'].unique()):
        count = seed_counts.get(seed_id, 0)
        comp = seed_completed.get(seed_id, 0)
        print(f"  {seed_id}: {comp}/{count} ケース完了")
    print()
    
    # 詳細結果を保存
    output_csv = os.path.join(args.output_dir, 'aor_vs_cutoff.csv')
    df.to_csv(output_csv, index=False)
    print(f"詳細結果を保存: {output_csv}")
    
    # 全シード平均のマトリックス形式で保存
    matrix = create_matrix_view(df, seed_id=None)
    if not matrix.empty:
        matrix_csv = os.path.join(args.output_dir, 'aor_matrix.csv')
        matrix.to_csv(matrix_csv)
        print(f"マトリックス形式（全シード平均）を保存: {matrix_csv}")
        
        print()
        print("安息角 [deg] マトリックス（全シード平均）:")
        print("-" * 60)
        print(matrix.to_string())
    
    # 各シードのマトリックスも保存
    seed_ids = sorted(df['seed_id'].unique())
    if len(seed_ids) > 1:
        print()
        print("シード別マトリックス:")
        for seed_id in seed_ids:
            matrix_seed = create_matrix_view(df, seed_id=seed_id)
            if not matrix_seed.empty:
                matrix_seed_csv = os.path.join(args.output_dir, f'aor_matrix_{seed_id}.csv')
                matrix_seed.to_csv(matrix_seed_csv)
                print(f"  {seed_id}: {matrix_seed_csv}")
    
    print()
    print("=" * 60)
    print("集約完了")
    print("=" * 60)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
