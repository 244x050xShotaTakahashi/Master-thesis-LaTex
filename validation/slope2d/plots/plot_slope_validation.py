#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2次元斜面検証用DEMシミュレーション結果のプロットスクリプト

使用方法:
    python plots/plot_slope_validation.py

出力:
    plots/slope_velocity.png - 速度の時系列比較
    plots/slope_trajectory.png - 粒子の軌跡
    plots/slope_error.png - 相対誤差の時系列
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
    return df

def read_slope_angle_from_input(input_path='input/slope_input.dat'):
    """入力ファイルから斜面角度 SLOPE_ANGLE [度] を読み取り、ラジアンに変換して返す。
    コメント行(#/!)と空行を無視し、見つからない場合は安全なフォールバック値を返す。
    """
    PI_VAL = np.pi
    fallback = 0.785398163  # π/4 rad (45度)
    
    try:
        if not os.path.exists(input_path):
            print(f"警告: {input_path} が見つかりません。フォールバック値 {fallback} rad を使用します。")
            return fallback
        
        with open(input_path, 'r', encoding='utf-8') as f:
            for raw in f:
                line = raw.strip()
                # コメント行と空行をスキップ
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0] == 'SLOPE_ANGLE':
                    try:
                        angle_deg = float(parts[1])
                        angle_rad = angle_deg / 180.0 * PI_VAL
                        print(f"SLOPE_ANGLE 読み込み: {angle_deg}° = {angle_rad:.6f} rad")
                        return angle_rad
                    except ValueError:
                        print(f"警告: SLOPE_ANGLE の値 '{parts[1]}' を数値に変換できません")
                        continue
        
        # SLOPE_ANGLEが見つからなかった場合
        print(f"警告: {input_path} に SLOPE_ANGLE が見つかりません。フォールバック値 {fallback} rad を使用します。")
        return fallback
        
    except Exception as e:
        print(f"警告: {input_path} の読み込み中にエラーが発生しました: {e}")
        print(f"フォールバック値 {fallback} rad を使用します。")
        return fallback

def add_velocity_components(df, slope_angle):
    """データフレームに斜面平行/法線方向の速度成分を追加する。
    v_parallel = vx*cosθ + vz*sinθ
    v_normal   = -vx*sinθ + vz*cosθ
    """
    if 'vx' not in df.columns or 'vz' not in df.columns:
        raise KeyError("DataFrameに 'vx' または 'vz' 列が存在しません")
    c = np.cos(slope_angle)
    s = np.sin(slope_angle)
    df['v_parallel'] = df['vx'] * c + df['vz'] * s
    df['v_normal'] = -df['vx'] * s + df['vz'] * c

def plot_velocity_comparison(df, output_file='plots/slope_velocity.png'):
    """速度の数値解と理論解を比較プロット"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # 上段: 速度の時系列
    ax1.plot(df['time'], df['v_magnitude'], 'b-', linewidth=2, label='Numerical (DEM)')
    ax1.plot(df['time'], df['omega'] * 0.01, 'g--', linewidth=2, label='Velocity')
    ax1.plot(df['time'], df['v_theory'], 'r--', linewidth=2, label='Theoretical')
    ax1.set_xlabel('Time [s]', fontsize=12)
    ax1.set_ylabel('Velocity Magnitude [m/s]', fontsize=12)
    ax1.set_title('Particle Velocity on Slope: Numerical vs Theoretical', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    # ax1.set_xlim(0, 0.0002)
    # ax1.set_ylim(0, 0.005)
    
    # 下段: 相対誤差
    ax2.plot(df['time'], df['error_percent'], 'g-', linewidth=1.5)
    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Relative Error [%]', fontsize=12)
    ax2.set_title('Relative Error', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='1% error')
    ax2.legend(fontsize=10)
    # ax2.set_xlim(0, 0.002)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"速度比較プロット保存: {output_file}")
    plt.close()

def plot_trajectory(df, output_file='plots/slope_trajectory.png'):
    """粒子の軌跡をプロット"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 粒子の軌跡
    scatter = ax.scatter(df['x'], df['z'], c=df['time'], cmap='viridis', 
                        s=20, alpha=0.6, edgecolors='none')
    
    # カラーバー
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Time [s]', fontsize=11)
    
    # 斜面の描画（最初の接触点から推定）
    contact_data = df[df['contact'] == 1]
    if len(contact_data) > 0:
        x_contact = contact_data['x'].iloc[0]
        z_contact = contact_data['z'].iloc[0]
        
        # 斜面の傾きを推定
        x_slope = np.linspace(0, df['x'].max() * 1.2, 100)
        # 斜面方程式: z = tan(θ) * x （原点を通る）
        # 接触点から傾きを逆算
        if x_contact > 1e-10:
            slope_tan = z_contact / x_contact
            z_slope = slope_tan * x_slope
            ax.plot(x_slope, z_slope, 'k--', linewidth=2, alpha=0.7, label='Slope')
    
    ax.set_xlabel('X Position [m]', fontsize=12)
    ax.set_ylabel('Z Position [m]', fontsize=12)
    ax.set_title('Particle Trajectory on Slope', fontsize=14, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"軌跡プロット保存: {output_file}")
    plt.close()

def plot_angular_velocity(df, output_file='plots/slope_angular_velocity.png'):
    """角速度の時系列をプロット"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(df['time'], df['omega'], 'b-', linewidth=2, label='Angular Velocity')
    
    # 転がり条件の理論角速度: ω = v/r
    # 理論角速度を計算（半径が必要なので近似）
    # 入力ファイルから r = 0.01 m と仮定
    r_approx = 0.01
    omega_theory = df['v_theory'] / r_approx
    ax.plot(df['time'], omega_theory, 'r--', linewidth=2, label='Theoretical (rolling)')
    
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Angular Velocity [rad/s]', fontsize=12)
    ax.set_title('Angular Velocity: Numerical vs Theoretical', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    # ax.set_ylim(0, 0.1)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"角速度プロット保存: {output_file}")
    plt.close()

def plot_velocity_components(df, output_file='plots/slope_velocity_components.png'):
    """斜面平行/法線方向の速度成分の時系列をプロット"""
    if 'v_parallel' not in df.columns or 'v_normal' not in df.columns:
        raise KeyError("DataFrameに 'v_parallel' または 'v_normal' 列が存在しません")

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(df['time'], df['v_parallel'], 'b-', linewidth=2, label='Parallel to slope')
    ax.plot(df['time'], df['v_normal'], 'r-', linewidth=2, label='Normal to slope')

    ax.axhline(0.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Velocity [m/s]', fontsize=12)
    ax.set_title('Velocity Components relative to Slope', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.02)
    ax.set_ylim(-0.025, 0.025)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"成分速度プロット保存: {output_file}")
    plt.close()

def print_statistics(df):
    """統計情報を表示"""
    print("\n" + "="*60)
    print("統計情報")
    print("="*60)
    
    # 最終時刻のデータ
    final_data = df.iloc[-1]
    print(f"最終時刻: {final_data['time']:.4f} s")
    print(f"最終速度（数値解）: {final_data['v_magnitude']:.6f} m/s")
    print(f"最終速度（理論解）: {final_data['v_theory']:.6f} m/s")
    print(f"最終相対誤差: {final_data['error_percent']:.8f} %")
    
    # 統計
    print(f"\n平均相対誤差: {df['error_percent'].mean():.4f} %")
    print(f"最大相対誤差: {df['error_percent'].max():.4f} %")
    print(f"最小相対誤差: {df['error_percent'].min():.4f} %")
    
    # 接触状態の統計
    contact_ratio = df['contact'].sum() / len(df) * 100
    print(f"\n接触時間の割合: {contact_ratio:.2f} %")
    
    print("="*60)

def main():
    """メイン処理"""
    print("="*60)
    print("2次元斜面検証用DEMシミュレーション結果プロット")
    print("="*60)
    
    # データ読み込み
    df = load_data()
    
    # plotsディレクトリの作成
    os.makedirs('plots', exist_ok=True)
    
    # 斜面角度を入力から取得し、速度成分を追加
    angle = read_slope_angle_from_input('input/slope_input.dat')
    add_velocity_components(df, angle)

    # 各種プロット作成
    plot_velocity_comparison(df)
    plot_trajectory(df)
    plot_angular_velocity(df)
    plot_velocity_components(df)
    
    # 統計情報表示
    print_statistics(df)
    
    print("\nプロット完了!")
    print("="*60)

if __name__ == '__main__':
    main()



