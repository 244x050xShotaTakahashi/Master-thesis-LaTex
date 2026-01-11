#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2次元斜面検証用DEMシミュレーション - 接触力プロットスクリプト

使用方法:
    python plots/plot_slope_forces.py

出力:
    plots/slope_forces.png - 接触力の時系列（法線力・接線力）
    plots/slope_force_ratio.png - 接線力/法線力の比（摩擦係数との比較）
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# 日本語フォント設定
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

def load_data(filename='data/slope_trace.csv'):
    """CSVファイルからデータを読み込む"""
    if not os.path.exists(filename):
        print(f"エラー: {filename} が見つかりません")
        sys.exit(1)
    
    df = pd.read_csv(filename)
    print(f"データ読み込み完了: {len(df)} 行")
    
    # 力の列が存在するか確認
    if 'force_normal' not in df.columns or 'force_tangent' not in df.columns:
        print("エラー: CSVファイルに 'force_normal' または 'force_tangent' 列が存在しません")
        print(f"利用可能な列: {df.columns.tolist()}")
        sys.exit(1)
    
    return df

def read_friction_coeff_from_input(input_path='input/slope_input.dat'):
    """入力ファイルから摩擦係数 FRICTION_COEFF を読み取る。
    見つからない場合はデフォルト値を返す。
    """
    fallback = 0.2
    
    try:
        if not os.path.exists(input_path):
            print(f"警告: {input_path} が見つかりません。フォールバック値 {fallback} を使用します。")
            return fallback
        
        with open(input_path, 'r', encoding='utf-8') as f:
            for raw in f:
                line = raw.strip()
                # コメント行と空行をスキップ
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0] == 'FRICTION_COEFF':
                    try:
                        mu = float(parts[1])
                        print(f"FRICTION_COEFF 読み込み: {mu}")
                        return mu
                    except ValueError:
                        print(f"警告: FRICTION_COEFF の値 '{parts[1]}' を数値に変換できません")
                        continue
        
        # FRICTION_COEFFが見つからなかった場合
        print(f"警告: {input_path} に FRICTION_COEFF が見つかりません。フォールバック値 {fallback} を使用します。")
        return fallback
        
    except Exception as e:
        print(f"警告: {input_path} の読み込み中にエラーが発生しました: {e}")
        print(f"フォールバック値 {fallback} を使用します。")
        return fallback

def plot_forces(df, output_file='plots/slope_forces.png'):
    """接触力（法線力・接線力）の時系列をプロット"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # 上段: 法線力
    ax1.plot(df['time'], df['force_normal'], 'b-', linewidth=1.5, label='Normal Force')
    ax1.set_xlabel('Time [s]', fontsize=12)
    ax1.set_ylabel('Force [N]', fontsize=12)
    ax1.set_title('Contact Force: Normal Component', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # 下段: 接線力（絶対値）
    ax2.plot(df['time'], np.abs(df['force_tangent']), 'r-', linewidth=1.5, label='Tangential Force (abs)')
    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Force [N]', fontsize=12)
    ax2.set_title('Contact Force: Tangential Component (Absolute Value)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"接触力プロット保存: {output_file}")
    plt.close()

def plot_force_ratio(df, mu, output_file='plots/slope_force_ratio.png'):
    """接線力/法線力の比と摩擦係数を比較プロット"""
    # 接触時のデータのみ抽出（法線力が0でない）
    contact_df = df[df['force_normal'] > 1e-12].copy()
    
    if len(contact_df) == 0:
        print("警告: 接触データがありません。force_ratio プロットをスキップします。")
        return
    
    # 力の比を計算
    contact_df['force_ratio'] = np.abs(contact_df['force_tangent']) / contact_df['force_normal']
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 力の比
    ax.plot(contact_df['time'], contact_df['force_ratio'], 'b-', linewidth=1.5, 
            label='|F_tangent| / F_normal', alpha=0.8)
    
    # 摩擦係数（理論的な上限）
    ax.axhline(y=mu, color='r', linestyle='--', linewidth=2, 
               label=f'Friction Coefficient μ = {mu:.3f}', alpha=0.7)
    
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Force Ratio [-]', fontsize=12)
    ax.set_title('Ratio of Tangential to Normal Force', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, mu * 1.5)  # 摩擦係数の1.5倍までの範囲
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"力比プロット保存: {output_file}")
    plt.close()

def plot_forces_combined(df, output_file='plots/slope_forces_combined.png'):
    """法線力と接線力を1つのグラフに重ねてプロット"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 両方の力をプロット
    ax.plot(df['time'], df['force_normal'], 'b-', linewidth=1.5, label='Normal Force', alpha=0.8)
    ax.plot(df['time'], np.abs(df['force_tangent']), 'r-', linewidth=1.5, 
            label='Tangential Force (abs)', alpha=0.8)
    
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Force [N]', fontsize=12)
    ax.set_title('Contact Forces: Normal and Tangential Components', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"複合力プロット保存: {output_file}")
    plt.close()

def print_force_statistics(df, mu):
    """力に関する統計情報を表示"""
    # 接触時のデータのみ
    contact_df = df[df['contact'] == 1]
    
    if len(contact_df) == 0:
        print("\n警告: 接触データがありません")
        return
    
    print("\n" + "="*60)
    print("接触力の統計情報")
    print("="*60)
    
    # 法線力の統計
    print("\n法線力 (Normal Force):")
    print(f"  平均値: {contact_df['force_normal'].mean():.6e} N")
    print(f"  最大値: {contact_df['force_normal'].max():.6e} N")
    print(f"  最小値: {contact_df['force_normal'].min():.6e} N")
    print(f"  標準偏差: {contact_df['force_normal'].std():.6e} N")
    
    # 接線力の統計
    print("\n接線力 (Tangential Force):")
    print(f"  平均値: {contact_df['force_tangent'].mean():.6e} N")
    print(f"  最大値 (絶対値): {np.abs(contact_df['force_tangent']).max():.6e} N")
    print(f"  最小値 (絶対値): {np.abs(contact_df['force_tangent']).min():.6e} N")
    print(f"  標準偏差: {contact_df['force_tangent'].std():.6e} N")
    
    # 力の比の統計
    valid_contact = contact_df[contact_df['force_normal'] > 1e-12]
    if len(valid_contact) > 0:
        force_ratios = np.abs(valid_contact['force_tangent']) / valid_contact['force_normal']
        print("\n力の比 |F_tangent| / F_normal:")
        print(f"  平均値: {force_ratios.mean():.6f}")
        print(f"  最大値: {force_ratios.max():.6f}")
        print(f"  最小値: {force_ratios.min():.6f}")
        print(f"  標準偏差: {force_ratios.std():.6f}")
        print(f"  摩擦係数 μ: {mu:.6f}")
        print(f"  平均比/μ: {force_ratios.mean()/mu:.6f} (1.0に近いほど滑り条件)")
    
    print("="*60)

def main():
    """メイン処理"""
    print("="*60)
    print("2次元斜面検証用DEM - 接触力プロット")
    print("="*60)
    
    # データ読み込み
    df = load_data()
    
    # plotsディレクトリの作成
    os.makedirs('plots', exist_ok=True)
    
    # 摩擦係数を入力ファイルから取得
    mu = read_friction_coeff_from_input('input/slope_input.dat')
    
    # 各種プロット作成
    plot_forces(df)
    plot_force_ratio(df, mu)
    plot_forces_combined(df)
    
    # 統計情報表示
    print_force_statistics(df, mu)
    
    print("\nプロット完了!")
    print("="*60)

if __name__ == '__main__':
    main()



