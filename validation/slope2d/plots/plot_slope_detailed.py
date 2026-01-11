#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PEM-SLOPE用詳細プロットスクリプト

pem2d_slope_kv.f90の出力データから、詳細な速度成分・角速度・力のプロットを作成します。

使用方法:
    python plots/plot_slope_detailed.py [--csv data/pem_slope_trace.csv] [--out plots/]

出力:
    plots/slope_velocity.png - 速度の時系列比較（数値解 vs 理論解）
    plots/slope_velocity_components.png - 速度成分（vx, vz, 斜面平行・法線方向）
    plots/slope_forces.png - 接触力（法線力Fn, 接線力Ft）
    plots/slope_angular_velocity.png - 角速度（omega）
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# 日本語フォント設定
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False


def load_data(filename='data/pem_slope_trace.csv'):
    """CSVファイルからデータを読み込む"""
    if not os.path.exists(filename):
        print(f"エラー: {filename} が見つかりません")
        sys.exit(1)
    
    # 型推定の分割読み込みを無効化し、混在型警告を回避
    df = pd.read_csv(filename, low_memory=False)
    print(f"データ読み込み完了: {len(df)} 行")
    
    # 必要な列の存在確認
    required_cols = ['time', 'vx', 'vz', 'omega', 'v_tangent', 'v_theory', 'error_percent', 'Fn', 'Ft']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"エラー: 必要な列が不足しています: {missing_cols}")
        print(f"利用可能な列: {df.columns.tolist()}")
        sys.exit(1)
    
    # 数値列を明示的に float に変換（混在型をNaNへ強制）
    numeric_candidates = required_cols + [
        'u_generalized', 'u_generalized_theory', 'vg_theory', 'Ft_theory'
    ]
    for col in numeric_candidates:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    return df


def read_slope_angle_from_input(input_path='input/slope_input.dat'):
    """入力ファイルから斜面角度 SLOPE_ANGLE [度] を読み取り、ラジアンに変換して返す。
    コメント行(#/!)と空行を無視し、見つからない場合は安全なフォールバック値を返す。
    """
    PI_VAL = np.pi
    fallback = 0.523599  # π/6 rad (30度)
    
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


def read_particle_radius_from_input(input_path='input/slope_input.dat'):
    """入力ファイルから粒子半径 PARTICLE_RADIUS [m] を読み取る。
    見つからない場合はデフォルト値を返す。
    """
    fallback = 0.01  # 1cm
    
    try:
        if not os.path.exists(input_path):
            print(f"警告: {input_path} が見つかりません。フォールバック値 {fallback} m を使用します。")
            return fallback
        
        with open(input_path, 'r', encoding='utf-8') as f:
            for raw in f:
                line = raw.strip()
                # コメント行と空行をスキップ
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0] == 'PARTICLE_RADIUS':
                    try:
                        r = float(parts[1])
                        print(f"PARTICLE_RADIUS 読み込み: {r} m")
                        return r
                    except ValueError:
                        print(f"警告: PARTICLE_RADIUS の値 '{parts[1]}' を数値に変換できません")
                        continue
        
        # PARTICLE_RADIUSが見つからなかった場合
        print(f"警告: {input_path} に PARTICLE_RADIUS が見つかりません。フォールバック値 {fallback} m を使用します。")
        return fallback
        
    except Exception as e:
        print(f"警告: {input_path} の読み込み中にエラーが発生しました: {e}")
        print(f"フォールバック値 {fallback} m を使用します。")
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
    ax1.plot(df['time'], df['v_tangent'], 'b-', linewidth=2, label='Numerical (v_tangent)')
    ax1.plot(df['time'], df['v_theory'], 'r--', linewidth=2, label='Theoretical')
    ax1.set_xlabel('Time [s]', fontsize=12)
    ax1.set_ylabel('Velocity [m/s]', fontsize=12)
    ax1.set_title('Particle Velocity on Slope: Numerical vs Theoretical', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    # ax1.set_xlim(0.0, 0.002)  
    # ax1.set_ylim(-0.01, 0.001)
    
    # 下段: 相対誤差
    ax2.plot(df['time'], df['error_percent'], 'g-', linewidth=1.5)
    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Relative Error [%]', fontsize=12)
    ax2.set_title('Relative Error', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5, label='1% error')
    ax2.legend(fontsize=10)
    # ax2.set_xlim(0.0, 0.002)  
    # ax2.set_ylim(-0.001, 0.001)
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"速度比較プロット保存: {output_file}")
    plt.close()


def plot_velocity_components(df, output_file='plots/slope_velocity_components.png'):
    """斜面平行/法線方向の速度成分の時系列をプロット"""
    if 'v_parallel' not in df.columns or 'v_normal' not in df.columns:
        raise KeyError("DataFrameに 'v_parallel' または 'v_normal' 列が存在しません")

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(df['time'], df['vx'], 'b-', linewidth=1.5, label='vx (horizontal)', alpha=0.7)
    ax.plot(df['time'], df['vz'], 'r-', linewidth=1.5, label='vz (vertical)', alpha=0.7)
    ax.plot(df['time'], df['v_parallel'], 'g-', linewidth=2, label='Parallel to slope')
    ax.plot(df['time'], df['v_normal'], 'm-', linewidth=2, label='Normal to slope')

    ax.axhline(0.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Velocity [m/s]', fontsize=12)
    ax.set_title('Velocity Components relative to Slope', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.0, 0.002)  
    ax.set_ylim(-0.01, 0.001)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"成分速度プロット保存: {output_file}")
    plt.close()


def plot_forces(df, output_file='plots/slope_forces.png'):
    """接触力（法線力・接線力）の時系列をプロット"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # 上段: 法線力
    ax1.plot(df['time'], df['Fn'], 'b-', linewidth=1.5, label='Normal Force (Fn)')
    ax1.set_xlabel('Time [s]', fontsize=12)
    ax1.set_ylabel('Force [N]', fontsize=12)
    ax1.set_title('Contact Force: Normal Component', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0.0, 0.002)  
    ax1.set_ylim(-0.001, 0.2)
    
    # 下段: 接線力（絶対値）
    ax2.plot(df['time'], np.abs(df['Ft']), 'r-', linewidth=1.5, label='Tangential Force |Ft|')
    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Force [N]', fontsize=12)
    ax2.set_title('Contact Force: Tangential Component (Absolute Value)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0.0, 0.002)  
    ax2.set_ylim(-0.001, 0.06)
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"接触力プロット保存: {output_file}")
    plt.close()


def plot_tangential_force_compare(df, output_file='plots/slope_forces_compare.png'):
    """接線力の数値 vs 理論を比較プロット（符号付き）"""
    if 'Ft' not in df.columns or 'Ft_theory' not in df.columns:
        print("スキップ: Ft_theory がCSVに無いため、接線力比較プロットを作成しません。")
        return

    fig, ax = plt.subplots(figsize=(12, 5))
    sub = df[['time','Ft','Ft_theory']].dropna()
    ax.plot(sub['time'], sub['Ft'], 'b-', linewidth=1.8, label='Ft (numerical)')
    ax.plot(sub['time'], sub['Ft_theory'], 'r--', linewidth=1.8, label='Ft_theory')
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Tangential Force [N]', fontsize=12)
    ax.set_title('Tangential Force: Numerical vs Theory', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.0, 0.002)  

    # RMSE を表示
    try:
        rmse = float(np.sqrt(np.mean((sub['Ft'] - sub['Ft_theory'])**2)))
        ax.text(0.02, 0.95, f"RMSE = {rmse:.3e} N", transform=ax.transAxes,
                fontsize=10, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    except Exception:
        pass

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"接線力比較プロット保存: {output_file}")
    plt.close()


def plot_generalized_velocity_compare(df, output_file='plots/slope_generalized_velocity.png'):
    """一般化速度 u' の数値 vs 理論を比較プロット
    数値: u_generalized = v_parallel + r*omega（CSV出力列）
    理論: vg_theory = u'(t)
    """
    if 'u_generalized' not in df.columns or 'vg_theory' not in df.columns:
        print("スキップ: u_generalized または vg_theory がCSVに無いため、一般化速度比較プロットを作成しません。")
        return

    fig, ax = plt.subplots(figsize=(12, 5))
    sub = df[['time','u_generalized','vg_theory']].dropna()
    ax.plot(sub['time'], sub['u_generalized'], 'g-', linewidth=1.8, label="u' (numerical)")
    ax.plot(sub['time'], -1.0 * sub['vg_theory'], 'm--', linewidth=1.8, label="u' (theory)")
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel("Generalized Velocity u' [m/s]", fontsize=12)
    ax.set_title("Generalized Velocity u': Numerical vs Theory", fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.0, 0.002)  
    ax.set_ylim(-0.0002, 0.0002)

    # RMSE を表示
    try:
        rmse = float(np.sqrt(np.mean((sub['u_generalized'] - sub['vg_theory'])**2)))
        ax.text(0.02, 0.95, f"RMSE = {rmse:.3e} m/s", transform=ax.transAxes,
                fontsize=10, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    except Exception:
        pass

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"一般化速度比較プロット保存: {output_file}")
    plt.close()

def plot_angular_velocity(df, particle_radius, output_file='plots/slope_angular_velocity.png'):
    """角速度の時系列をプロット"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(df['time'], df['omega'], 'b-', linewidth=2, label='Angular Velocity (numerical)')
    
    # 転がり条件の理論角速度: ω = v_theory / r
    omega_theory = df['v_theory'] / particle_radius
    ax.plot(df['time'], omega_theory, 'r--', linewidth=2, label='Theoretical (rolling: v/r)')
    
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Angular Velocity [rad/s]', fontsize=12)
    ax.set_title('Angular Velocity: Numerical vs Theoretical', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"角速度プロット保存: {output_file}")
    plt.close()


def print_statistics(df, mu):
    """統計情報を表示"""
    print("\n" + "="*60)
    print("統計情報")
    print("="*60)
    
    # 最終時刻のデータ
    final_data = df.iloc[-1]
    print(f"最終時刻: {final_data['time']:.6f} s")
    print(f"最終速度（数値解）: {final_data['v_tangent']:.6e} m/s")
    print(f"最終速度（理論解）: {final_data['v_theory']:.6e} m/s")
    print(f"最終相対誤差: {final_data['error_percent']:.8f} %")
    
    # 統計
    print(f"\n平均相対誤差: {df['error_percent'].mean():.4f} %")
    print(f"最大相対誤差: {df['error_percent'].max():.4f} %")
    print(f"最小相対誤差: {df['error_percent'].min():.4f} %")
    
    # 接触状態の統計
    if 'contact' in df.columns:
        contact_ratio = df['contact'].sum() / len(df) * 100
        print(f"\n接触時間の割合: {contact_ratio:.2f} %")
    
    # 力の統計
    contact_df = df[df['Fn'] > 1e-12]
    if len(contact_df) > 0:
        print(f"\n接触力の統計 (接触時のみ):")
        print(f"  法線力 Fn:")
        print(f"    平均値: {contact_df['Fn'].mean():.6e} N")
        print(f"    最大値: {contact_df['Fn'].max():.6e} N")
        print(f"  接線力 Ft:")
        print(f"    平均値 (絶対値): {np.abs(contact_df['Ft']).mean():.6e} N")
        print(f"    最大値 (絶対値): {np.abs(contact_df['Ft']).max():.6e} N")
        
        # 力の比の統計
        force_ratios = np.abs(contact_df['Ft']) / contact_df['Fn']
        print(f"  力の比 |Ft| / Fn:")
        print(f"    平均値: {force_ratios.mean():.6f}")
        print(f"    最大値: {force_ratios.max():.6f}")
        print(f"  摩擦係数 μ: {mu:.6f}")
    
    print("="*60)


def main():
    """メイン処理"""
    parser = argparse.ArgumentParser(description='PEM-SLOPE用詳細プロットスクリプト')
    parser.add_argument('--csv', default='data/pem_slope_trace.csv', 
                        help='入力CSVファイルのパス (デフォルト: data/pem_slope_trace.csv)')
    parser.add_argument('--out', default='plots/', 
                        help='出力ディレクトリ (デフォルト: plots/)')
    parser.add_argument('--input', default='input/slope_input.dat',
                        help='入力パラメータファイルのパス (デフォルト: input/slope_input.dat)')
    args = parser.parse_args()
    
    print("="*60)
    print("PEM-SLOPE用詳細プロットスクリプト")
    print("="*60)
    
    # データ読み込み
    df = load_data(args.csv)
    
    # 出力ディレクトリの作成
    os.makedirs(args.out, exist_ok=True)
    
    # 入力パラメータを読み込み
    slope_angle = read_slope_angle_from_input(args.input)
    friction_coeff = read_friction_coeff_from_input(args.input)
    particle_radius = read_particle_radius_from_input(args.input)
    
    # 速度成分を追加
    add_velocity_components(df, slope_angle)
    
    # 各種プロット作成（既存）
    plot_velocity_comparison(df, os.path.join(args.out, 'slope_velocity.png'))
    plot_velocity_components(df, os.path.join(args.out, 'slope_velocity_components.png'))
    plot_forces(df, os.path.join(args.out, 'slope_forces.png'))
    plot_angular_velocity(df, particle_radius, os.path.join(args.out, 'slope_angular_velocity.png'))

    # 追加: 理論 vs 数値の比較プロット（列がある場合のみ）
    if 'Ft_theory' in df.columns:
        plot_tangential_force_compare(df, os.path.join(args.out, 'slope_forces_compare.png'))
    else:
        print('注意: Ft_theory 列が無いため、接線力比較プロットはスキップします。')

    if 'u_generalized' in df.columns and 'vg_theory' in df.columns:
        plot_generalized_velocity_compare(df, os.path.join(args.out, 'slope_generalized_velocity.png'))
    else:
        print("注意: u_generalized または vg_theory 列が無いため、一般化速度比較プロットはスキップします。")
    
    # 統計情報表示
    print_statistics(df, friction_coeff)
    
    print("\nプロット完了!")
    print("="*60)


if __name__ == '__main__':
    main()



