#!/usr/bin/env python3
"""
particles.csvから最終ステップの電荷・半径データを抽出して
charges_step*.datとradii_step*.datファイルを生成するスクリプト

使用方法:
    python scripts/extract_distribution_from_csv.py --file particles.csv --output-dir .
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def extract_distribution_from_csv(csv_file: str, output_dir: str = "."):
    """
    particles.csvから最終ステップの電荷・半径データを抽出
    
    Parameters
    ----------
    csv_file : str
        particles.csvファイルのパス
    output_dir : str
        出力ディレクトリ
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # CSVファイルを読み込み
    df = pd.read_csv(csv_file, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    
    # 必要なカラムの確認
    required_cols = ['step', 'time', 'radius']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSVファイルに必要なカラムが含まれていません: {required_cols}")
    
    # 最終ステップを取得
    final_step = df['step'].max()
    final_time = df[df['step'] == final_step]['time'].iloc[0]
    
    # 最終ステップのデータを抽出
    final_df = df[df['step'] == final_step]
    
    # 電荷データの抽出と保存
    if 'charge' in final_df.columns:
        charges = final_df['charge'].values
        charges_file = output_path / f"charges_step{int(final_step)}.dat"
        with open(charges_file, 'w') as f:
            f.write(f"# 分布: particles.csv charge (final step, step={int(final_step)}, time={final_time:.6e})\n")
            for charge in charges:
                f.write(f"{charge:.16e}\n")
        print(f"電荷データを保存しました: {charges_file} ({len(charges)} 個)")
    else:
        print("警告: particles.csvに'charge'列がありません。電荷データは生成されません。")
    
    # 半径データの抽出と保存（直径として保存）
    radii = final_df['radius'].values
    diameters = radii * 2.0  # 直径に変換
    radii_file = output_path / f"radii_step{int(final_step)}.dat"
    with open(radii_file, 'w') as f:
        f.write(f"# 分布: particles.csv diameter (=2*radius) (final step, step={int(final_step)}, time={final_time:.6e})\n")
        for diameter in diameters:
            f.write(f"{diameter:.16e}\n")
    print(f"半径データを保存しました: {radii_file} ({len(diameters)} 個)")
    
    return int(final_step), final_time


def main():
    parser = argparse.ArgumentParser(
        description='particles.csvから最終ステップの電荷・半径データを抽出',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--file', '-f', type=str, required=True,
                        help='particles.csvファイルのパス')
    parser.add_argument('--output-dir', '-o', type=str, default='.',
                        help='出力ディレクトリ (デフォルト: .)')
    
    args = parser.parse_args()
    
    try:
        step, time = extract_distribution_from_csv(args.file, args.output_dir)
        print(f"\n完了: ステップ {step}, 時間 {time:.6e} s")
        return 0
    except Exception as e:
        print(f"エラー: {e}", file=__import__('sys').stderr)
        return 1


if __name__ == '__main__':
    exit(main())



