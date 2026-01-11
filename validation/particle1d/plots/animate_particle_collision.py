#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
粒子間衝突シミュレーションのアニメーション
particle_collision_validation.csv から粒子の動きを可視化
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.font_manager as fm
from matplotlib.patches import Circle


def setup_japanese_font():
    """日本語フォントを自動検出して設定"""
    candidates = [
        'Noto Sans CJK JP', 'Noto Sans JP', 'Noto Serif CJK JP',
        'IPAexGothic', 'IPAPGothic', 'VL Gothic', 'TakaoGothic', 'TakaoPGothic',
        'Yu Gothic', 'Hiragino Sans', 'MS Gothic'
    ]

    found_name = None
    for name in candidates:
        try:
            path = fm.findfont(fm.FontProperties(family=name), fallback_to_default=False)
            if os.path.exists(path):
                found_name = name
                break
        except Exception:
            pass

    if found_name is None:
        for path in fm.findSystemFonts(fontpaths=None, fontext='ttf'):
            base = os.path.basename(path).lower()
            if ('noto' in base and ('cjk' in base or 'jp' in base)) or \
               ('ipa' in base) or ('takao' in base) or ('vl-gothic' in base):
                try:
                    prop = fm.FontProperties(fname=path)
                    found_name = prop.get_name()
                    break
                except Exception:
                    continue

    if found_name:
        plt.rcParams['font.family'] = found_name
        plt.rcParams['font.sans-serif'] = [found_name, 'Noto Sans CJK JP', 'DejaVu Sans', 'sans-serif']
    else:
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Noto Sans CJK JP', 'DejaVu Sans', 'sans-serif']

    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.size'] = 10
    plt.rcParams['animation.embed_limit'] = 100  # HTML埋め込みサイズ上限を100MBに設定


def load_csv(path='data/particle_collision_validation.csv'):
    """CSVファイルを読み込む"""
    if not os.path.exists(path):
        print(f"エラー: {path} が見つかりません")
        return None
    try:
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception as e:
        print(f"読み込みエラー: {e}")
        return None


def calculate_particle_positions(df, mass1=1.0, mass2=1.0, radius1=1.0, radius2=1.0):
    """
    CSVデータから粒子の位置を計算
    
    Parameters:
    - df: CSVデータフレーム
    - mass1, mass2: 粒子の質量
    - radius1, radius2: 粒子の半径
    
    Returns:
    - x1_arr, x2_arr: 各時刻における粒子1, 2のx座標
    """
    # 初期条件
    x1_init = 3.0
    x2_init = x1_init + radius1 + radius2 + 0.001
    v1_init = 1.0
    v2_init = 0.0
    
    # 時刻配列
    t_arr = df['time'].values
    dt_arr = np.diff(t_arr, prepend=0.0)
    
    # オーバーラップと相対速度
    overlap_arr = df['overlap_numeric'].values
    v_rel_arr = df['v_rel_numeric'].values
    
    # 位置配列の初期化
    x1_arr = np.zeros(len(t_arr))
    x2_arr = np.zeros(len(t_arr))
    x1_arr[0] = x1_init
    x2_arr[0] = x2_init
    
    # 速度配列（蛙飛び法で半ステップずれている）
    v1 = v1_init
    v2 = v2_init
    
    # 接触力から速度と位置を更新
    for i in range(1, len(t_arr)):
        dt = dt_arr[i]
        
        # オーバーラップから接触力を計算（簡易版）
        # F = k * δ + c * v_rel
        k = 1.0  # stiffness_k
        e = 0.25  # restitution_coeff_e
        m_eff = (mass1 * mass2) / (mass1 + mass2)
        c = -2.0 * np.log(e) * np.sqrt(m_eff * k / (np.log(e)**2 + np.pi**2))
        
        overlap = overlap_arr[i-1]
        v_rel = v_rel_arr[i-1]
        
        if overlap > 0:
            force = k * overlap + c * v_rel
        else:
            force = 0.0
        
        # 加速度
        a1 = -force / mass1
        a2 = force / mass2
        
        # 速度更新（蛙飛び法の近似）
        v1 += a1 * dt
        v2 += a2 * dt
        
        # 位置更新
        x1_arr[i] = x1_arr[i-1] + v1 * dt
        x2_arr[i] = x2_arr[i-1] + v2 * dt
    
    return x1_arr, x2_arr


def create_animation(df, output_file='data/particle_collision_animation.mp4', fps=30, skip_frames=10):
    """
    アニメーションを作成
    
    Parameters:
    - df: データフレーム
    - output_file: 出力ファイル名
    - fps: フレームレート
    - skip_frames: フレームをスキップする数（データ量削減）
    """
    # パラメータ
    mass1 = mass2 = 1.0
    radius1 = radius2 = 1.0
    
    # 粒子位置を計算
    x1_arr, x2_arr = calculate_particle_positions(df, mass1, mass2, radius1, radius2)
    
    # データをスキップ（アニメーション高速化）
    indices = np.arange(0, len(df), skip_frames)
    t_arr = df['time'].values[indices]
    x1_arr = x1_arr[indices]
    x2_arr = x2_arr[indices]
    overlap_arr = df['overlap_numeric'].values[indices]
    v_rel_arr = df['v_rel_numeric'].values[indices]
    
    # プロット設定
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[2, 1, 1], hspace=0.3, wspace=0.3)
    
    # 上段: 粒子の動き
    ax_particles = fig.add_subplot(gs[0, :])
    ax_particles.set_xlim(0, 10)
    ax_particles.set_ylim(-2, 2)
    ax_particles.set_aspect('equal')
    ax_particles.set_xlabel('位置 (m)')
    ax_particles.set_ylabel('y (m)')
    ax_particles.set_title('粒子間衝突シミュレーション')
    ax_particles.grid(True, alpha=0.3)
    
    # 中段左: オーバーラップ
    ax_overlap = fig.add_subplot(gs[1, 0])
    ax_overlap.set_xlabel('時刻 (ms)')
    ax_overlap.set_ylabel('オーバーラップ (mm)')
    ax_overlap.set_title('オーバーラップの時間変化')
    ax_overlap.grid(True, alpha=0.3)
    
    # 中段右: 相対速度
    ax_velocity = fig.add_subplot(gs[1, 1])
    ax_velocity.set_xlabel('時刻 (ms)')
    ax_velocity.set_ylabel('相対速度 (m/s)')
    ax_velocity.set_title('相対速度の時間変化')
    ax_velocity.grid(True, alpha=0.3)
    
    # 下段: エネルギー
    ax_energy = fig.add_subplot(gs[2, :])
    ax_energy.set_xlabel('時刻 (ms)')
    ax_energy.set_ylabel('運動エネルギー (J)')
    ax_energy.set_title('運動エネルギーの時間変化')
    ax_energy.grid(True, alpha=0.3)
    
    # 粒子の円
    circle1 = Circle((x1_arr[0], 0), radius1, color='C0', alpha=0.7, label='粒子1')
    circle2 = Circle((x2_arr[0], 0), radius2, color='C1', alpha=0.7, label='粒子2')
    ax_particles.add_patch(circle1)
    ax_particles.add_patch(circle2)
    ax_particles.legend()
    
    # 時刻テキスト
    time_text = ax_particles.text(0.02, 0.95, '', transform=ax_particles.transAxes, 
                                    fontsize=12, verticalalignment='top',
                                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 軌跡プロット用の線
    line_overlap, = ax_overlap.plot([], [], 'C0-', lw=2)
    line_velocity, = ax_velocity.plot([], [], 'C0-', lw=2)
    
    # エネルギー計算（速度を数値微分で計算）
    dt_arr = np.diff(t_arr, prepend=t_arr[1] - t_arr[0])  # ゼロ除算回避
    dt_arr[0] = dt_arr[1]  # 最初の時間刻みを設定
    v1_arr = np.diff(x1_arr, prepend=x1_arr[0]) / dt_arr
    v2_arr = np.diff(x2_arr, prepend=x2_arr[0]) / dt_arr
    ke1_arr = 0.5 * mass1 * v1_arr**2
    ke2_arr = 0.5 * mass2 * v2_arr**2
    ke_total = ke1_arr + ke2_arr
    
    line_energy, = ax_energy.plot([], [], 'C2-', lw=2, label='総運動エネルギー')
    ax_energy.legend()
    
    def init():
        """初期化"""
        circle1.center = (x1_arr[0], 0)
        circle2.center = (x2_arr[0], 0)
        time_text.set_text('')
        line_overlap.set_data([], [])
        line_velocity.set_data([], [])
        line_energy.set_data([], [])
        return circle1, circle2, time_text, line_overlap, line_velocity, line_energy
    
    def animate(frame):
        """各フレームの更新"""
        # 粒子位置の更新
        circle1.center = (x1_arr[frame], 0)
        circle2.center = (x2_arr[frame], 0)
        
        # 接触を色で表現
        if overlap_arr[frame] > 0:
            circle1.set_color('red')
            circle2.set_color('red')
        else:
            circle1.set_color('C0')
            circle2.set_color('C1')
        
        # 時刻表示
        time_text.set_text(f'時刻: {t_arr[frame]*1e3:.2f} ms\n'
                           f'オーバーラップ: {overlap_arr[frame]*1e3:.4f} mm\n'
                           f'相対速度: {v_rel_arr[frame]:.4f} m/s')
        
        # オーバーラップのプロット
        line_overlap.set_data(t_arr[:frame+1] * 1e3, overlap_arr[:frame+1] * 1e3)
        ax_overlap.set_xlim(0, t_arr[-1] * 1e3)
        ax_overlap.set_ylim(np.min(overlap_arr) * 1e3 - 0.1, np.max(overlap_arr) * 1e3 + 0.1)
        
        # 速度のプロット
        line_velocity.set_data(t_arr[:frame+1] * 1e3, v_rel_arr[:frame+1])
        ax_velocity.set_xlim(0, t_arr[-1] * 1e3)
        ax_velocity.set_ylim(np.min(v_rel_arr) - 0.1, np.max(v_rel_arr) + 0.1)
        
        # エネルギーのプロット
        line_energy.set_data(t_arr[:frame+1] * 1e3, ke_total[:frame+1])
        ax_energy.set_xlim(0, t_arr[-1] * 1e3)
        ax_energy.set_ylim(0, np.max(ke_total) * 1.1)
        
        return circle1, circle2, time_text, line_overlap, line_velocity, line_energy
    
    # アニメーション作成
    print(f"アニメーション作成中... (フレーム数: {len(indices)})")
    anim = animation.FuncAnimation(fig, animate, init_func=init, 
                                    frames=len(indices), interval=1000/fps, 
                                    blit=True, repeat=True)
    
    # ファイル保存
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # HTML5として保存（ブラウザで再生可能）
    html_file = output_file.replace('.mp4', '.html')
    try:
        print(f"HTML5形式で保存中... {html_file}")
        with open(html_file, 'w') as f:
            f.write(anim.to_jshtml())
        print(f"HTML5保存完了: {html_file}")
        print(f"ブラウザで {html_file} を開いて再生してください")
    except Exception as e:
        print(f"HTML5保存エラー: {e}")
    
    # MP4として保存（ffmpegが必要）
    if os.system('which ffmpeg > /dev/null 2>&1') == 0:
        try:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, bitrate=1800)
            anim.save(output_file, writer=writer)
            print(f"MP4保存完了: {output_file}")
        except Exception as e:
            print(f"MP4保存エラー: {e}")
    else:
        print("ffmpegが見つかりません。MP4保存をスキップします。")
    
    return anim


def main():
    setup_japanese_font()
    df = load_csv()
    
    if df is None:
        return
    
    print("データ読み込み完了")
    print(f"データ点数: {len(df)}")
    print(f"時間範囲: {df['time'].min():.4e} ~ {df['time'].max():.4e} 秒")
    
    # アニメーション作成（データ量が多い場合はskip_framesを増やす）
    # 100フレームごとにスキップして高速化
    skip = max(1, len(df) // 500)  # 約500フレームに調整
    print(f"フレームスキップ数: {skip}")
    
    anim = create_animation(df, skip_frames=skip, fps=30)
    
    print("完了")


if __name__ == '__main__':
    main()

