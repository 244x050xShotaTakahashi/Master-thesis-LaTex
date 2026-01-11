#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
摩擦係数と安息角の関係をプロットするスクリプト

複数の摩擦係数での安息角測定結果を統合して、
摩擦係数 vs 安息角のグラフを作成します。
理論予測（tan(θ) ≈ μ）との比較も行います。
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys
from typing import List, Dict
import argparse


def read_results_from_directory(results_dir: str) -> pd.DataFrame:
    """
    結果ディレクトリから全ケースの測定結果を読み込む
    
    Args:
        results_dir: 結果ディレクトリのパス
        
    Returns:
        測定結果を含むDataFrame
    """
    results_dir = Path(results_dir)
    
    # 各ケースのディレクトリを検索
    case_dirs = sorted([d for d in results_dir.iterdir() if d.is_dir() and d.name.startswith('friction_')])
    
    if not case_dirs:
        raise ValueError(f"結果ディレクトリが見つかりません: {results_dir}")
    
    data_list = []
    
    for case_dir in case_dirs:
        # ケース名から摩擦係数を抽出
        case_name = case_dir.name
        try:
            friction_coeff = float(case_name.replace('friction_', ''))
        except ValueError:
            print(f"警告: ディレクトリ名から摩擦係数を抽出できません: {case_name}")
            continue
        
        # 結果CSVファイルを読み込み
        csv_file = case_dir / f"{case_name}_results.csv"
        
        if not csv_file.exists():
            print(f"警告: 結果ファイルが見つかりません: {csv_file}")
            continue
        
        # CSVファイルを読み込み（parameter, value形式）
        df = pd.read_csv(csv_file)
        
        # データを辞書に変換
        result_dict = {'friction_coeff': friction_coeff}
        for _, row in df.iterrows():
            param = row['parameter']
            value = row['value']
            
            # 数値データを変換
            if param in ['left_angle_deg', 'right_angle_deg', 'average_angle_deg', 
                        'std_angle_deg', 'left_r_squared', 'right_r_squared']:
                try:
                    result_dict[param] = float(value)
                except (ValueError, TypeError):
                    result_dict[param] = np.nan
        
        data_list.append(result_dict)
    
    if not data_list:
        raise ValueError("有効な結果データが見つかりませんでした")
    
    # DataFrameに変換
    results_df = pd.DataFrame(data_list)
    results_df = results_df.sort_values('friction_coeff')
    
    return results_df


def plot_friction_vs_angle(results_df: pd.DataFrame, output_file: str = None):
    """
    摩擦係数 vs 安息角のグラフを作成
    
    Args:
        results_df: 測定結果のDataFrame
        output_file: 出力ファイル名
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # データの準備
    mu = results_df['friction_coeff'].values
    avg_angle = results_df['average_angle_deg'].values
    std_angle = results_df.get('std_angle_deg', pd.Series([0]*len(mu))).values
    left_angle = results_df.get('left_angle_deg', pd.Series([np.nan]*len(mu))).values
    right_angle = results_df.get('right_angle_deg', pd.Series([np.nan]*len(mu))).values
    
    # 理論曲線: θ = arctan(μ)
    mu_theory = np.linspace(0, max(mu) * 1.1, 100)
    angle_theory = np.degrees(np.arctan(mu_theory))
    
    # === グラフ1: 平均安息角 vs 摩擦係数 ===
    ax1.errorbar(mu, avg_angle, yerr=std_angle, fmt='o-', 
                 color='darkblue', markersize=10, linewidth=2,
                 capsize=5, capthick=2, label='測定値（平均）', zorder=3)
    
    # 左右の斜面も表示
    ax1.scatter(mu, left_angle, marker='^', s=100, color='red', 
               alpha=0.6, label='左斜面', zorder=2)
    ax1.scatter(mu, right_angle, marker='v', s=100, color='blue', 
               alpha=0.6, label='右斜面', zorder=2)
    
    # 理論曲線
    ax1.plot(mu_theory, angle_theory, '--', color='gray', linewidth=2,
            label='理論: θ = arctan(μ)', zorder=1)
    
    ax1.set_xlabel('摩擦係数 μ [-]', fontsize=14, weight='bold')
    ax1.set_ylabel('安息角 θ [度]', fontsize=14, weight='bold')
    ax1.set_title('安息角と摩擦係数の関係', fontsize=16, weight='bold', pad=20)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(fontsize=11, loc='upper left', framealpha=0.9)
    ax1.set_xlim(0, max(mu) * 1.05)
    ax1.set_ylim(0, max(avg_angle) * 1.15)
    
    # === グラフ2: tan(θ) vs μ（理論的には線形関係）===
    tan_avg = np.tan(np.radians(avg_angle))
    tan_left = np.tan(np.radians(left_angle))
    tan_right = np.tan(np.radians(right_angle))
    
    ax2.scatter(mu, tan_avg, marker='o', s=120, color='darkblue', 
               label='測定値（平均）', zorder=3)
    ax2.scatter(mu, tan_left, marker='^', s=100, color='red', 
               alpha=0.6, label='左斜面', zorder=2)
    ax2.scatter(mu, tan_right, marker='v', s=100, color='blue', 
               alpha=0.6, label='右斜面', zorder=2)
    
    # 理論直線: tan(θ) = μ
    mu_line = np.linspace(0, max(mu) * 1.1, 100)
    ax2.plot(mu_line, mu_line, '--', color='gray', linewidth=2,
            label='理論: tan(θ) = μ', zorder=1)
    
    # 線形フィッティング
    valid_mask = ~np.isnan(tan_avg)
    if valid_mask.sum() >= 2:
        coeffs = np.polyfit(mu[valid_mask], tan_avg[valid_mask], 1)
        fit_line = np.poly1d(coeffs)
        ax2.plot(mu_line, fit_line(mu_line), '-', color='green', linewidth=2,
                label=f'フィット: tan(θ) = {coeffs[0]:.3f}μ + {coeffs[1]:.3f}', zorder=2)
    
    ax2.set_xlabel('摩擦係数 μ [-]', fontsize=14, weight='bold')
    ax2.set_ylabel('tan(θ) [-]', fontsize=14, weight='bold')
    ax2.set_title('tan(θ) と摩擦係数の関係（線形性チェック）', fontsize=16, weight='bold', pad=20)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(fontsize=11, loc='upper left', framealpha=0.9)
    ax2.set_xlim(0, max(mu) * 1.05)
    ax2.set_ylim(0, max(tan_avg[valid_mask]) * 1.15 if valid_mask.sum() > 0 else 1)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"グラフを保存しました: {output_file}")
    
    plt.close()


def create_summary_table(results_df: pd.DataFrame, output_file: str = None):
    """
    結果のサマリーテーブルを作成
    
    Args:
        results_df: 測定結果のDataFrame
        output_file: 出力ファイル名（テキストファイル）
    """
    summary_lines = []
    summary_lines.append("=" * 80)
    summary_lines.append("安息角測定実験 - 結果サマリー")
    summary_lines.append("=" * 80)
    summary_lines.append("")
    summary_lines.append(f"{'摩擦係数':>10} | {'左斜面[°]':>12} | {'右斜面[°]':>12} | {'平均[°]':>10} | {'標準偏差':>10} | {'tan(θ)':>10}")
    summary_lines.append("-" * 80)
    
    for _, row in results_df.iterrows():
        mu = row['friction_coeff']
        left = row.get('left_angle_deg', np.nan)
        right = row.get('right_angle_deg', np.nan)
        avg = row.get('average_angle_deg', np.nan)
        std = row.get('std_angle_deg', 0)
        tan_theta = np.tan(np.radians(avg)) if not np.isnan(avg) else np.nan
        
        summary_lines.append(f"{mu:10.2f} | {left:12.2f} | {right:12.2f} | {avg:10.2f} | {std:10.2f} | {tan_theta:10.3f}")
    
    summary_lines.append("=" * 80)
    summary_lines.append("")
    summary_lines.append("理論予測: θ = arctan(μ)")
    summary_lines.append("")
    
    # 理論値との比較
    summary_lines.append(f"{'摩擦係数':>10} | {'測定[°]':>10} | {'理論[°]':>10} | {'差[°]':>10} | {'相対誤差[%]':>12}")
    summary_lines.append("-" * 80)
    
    for _, row in results_df.iterrows():
        mu = row['friction_coeff']
        avg = row.get('average_angle_deg', np.nan)
        theory = np.degrees(np.arctan(mu))
        diff = avg - theory if not np.isnan(avg) else np.nan
        rel_error = (diff / theory * 100) if not np.isnan(diff) and theory != 0 else np.nan
        
        summary_lines.append(f"{mu:10.2f} | {avg:10.2f} | {theory:10.2f} | {diff:10.2f} | {rel_error:12.1f}")
    
    summary_lines.append("=" * 80)
    
    summary_text = "\n".join(summary_lines)
    
    # コンソールに出力
    print("\n" + summary_text)
    
    # ファイルに保存
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(summary_text)
        print(f"\nサマリーを保存しました: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='摩擦係数と安息角の関係をプロット')
    parser.add_argument('--input', '-i', default='results/friction_study',
                       help='結果ディレクトリ (デフォルト: results/friction_study)')
    parser.add_argument('--output', '-o', default='results/friction_study/friction_study_summary.png',
                       help='出力ファイル名 (デフォルト: results/friction_study/friction_study_summary.png)')
    parser.add_argument('--summary', '-s', default=None,
                       help='サマリーテキストファイル名（オプション）')
    
    args = parser.parse_args()
    
    try:
        print(f"結果を読み込んでいます: {args.input}")
        results_df = read_results_from_directory(args.input)
        print(f"読み込んだケース数: {len(results_df)}")
        print(f"摩擦係数の範囲: {results_df['friction_coeff'].min():.1f} - {results_df['friction_coeff'].max():.1f}")
        
        # グラフ作成
        print("\nグラフを作成中...")
        plot_friction_vs_angle(results_df, args.output)
        
        # サマリーテーブル作成
        summary_file = args.summary
        if summary_file is None:
            # デフォルトはグラフと同じディレクトリ
            output_path = Path(args.output)
            summary_file = str(output_path.parent / "summary.txt")
        
        create_summary_table(results_df, summary_file)
        
        print("\n処理が完了しました。")
        
    except Exception as e:
        print(f"エラーが発生しました: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()




