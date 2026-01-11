#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
粒子間衝突検証プログラム
particle_collision_validation.csv から理論解と数値解のオーバーラップを可視化
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def setup_japanese_font():
    """日本語フォントを自動検出して設定（なければ可能な範囲でフォールバック）"""
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

    # 2段階目: ファイル名から推測
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
        # 最低限のフォールバック
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Noto Sans CJK JP', 'DejaVu Sans', 'sans-serif']

    # マイナス記号の文字化け対策
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.size'] = 12


def load_csv(path='data/particle_collision_validation.csv'):
    if not os.path.exists(path):
        print(f"エラー: {path} が見つかりません")
        return None
    try:
        df = pd.read_csv(path)
        # 列名の正規化（余分な空白を除去）
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception as e:
        print(f"読み込みエラー: {e}")
        return None


def calculate_restitution_coefficient(df):
    """
    接触開始時と終了時の速度から跳ね返り係数を計算
    Returns: (e_numeric, v_before, v_after) or (None, None, None)
    """
    if df is None or df.empty:
        return None, None, None
    
    # v_rel_numeric を使用（CSVに含まれる相対速度）
    if 'v_rel_numeric' not in df.columns:
        print("警告: v_rel_numeric列が見つかりません")
        return None, None, None
    
    v_rel = df['v_rel_numeric'].values
    
    # 速度の符号変化を検出（正→負への変化点 = 最大オーバーラップ点）
    sign_changes = []
    for i in range(len(v_rel) - 1):
        if v_rel[i] > 0 and v_rel[i+1] <= 0:
            sign_changes.append(i)
    
    if len(sign_changes) == 0:
        print("警告: 速度の符号変化が見つかりませんでした")
        # フォールバック: 最初と最後の絶対値を使用
        v_before = abs(v_rel[0])
        v_after = abs(v_rel[-1])
    else:
        # 接触開始時の速度（最初のデータ点）
        v_before = abs(v_rel[0])
        
        # 接触終了時の速度（最後のデータ点）
        v_after = abs(v_rel[-1])
    
    # 跳ね返り係数 e = |v_after| / |v_before|
    if v_before > 0.0:
        e_numeric = v_after / v_before
    else:
        e_numeric = float('nan')
    
    return e_numeric, v_before, v_after


def plot_overlap(df, output_dir='data'):
    if df is None or df.empty:
        return

    # time, overlap_numeric, overlap_theory を使用
    t = df['time'].values
    num = df['overlap_numeric'].values
    theo = df['overlap_theory'].values
    err = num - theo
    # 速度（数値・理論）をCSVから取得
    v_num = df['v_rel_numeric'].values
    v_theo = df['v_rel_theory'].values

    os.makedirs(output_dir, exist_ok=True)

    # 上段: オーバーラップ比較、下段: 速度比較
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    ax1.plot(t * 1e3, num * 1e3, label='数値解', color='C0', lw=2)
    ax1.plot(t * 1e3, theo * 1e3, label='理論解', color='C3', lw=2, ls='--')
    ax1.set_ylabel('オーバーラップ (mm)')
    ax1.set_title('粒子間衝突：理論解と数値解のオーバーラップ比較')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(t * 1e3, v_num, label='数値解 相対速度', color='C0', lw=2)
    ax2.plot(t * 1e3, v_theo, label='理論解 相対速度', color='C3', lw=2, ls='--')
    ax2.set_xlabel('接触開始からの時間 (ms)')
    ax2.set_ylabel('相対速度 (m/s)')
    ax2.set_title('粒子間衝突：理論解と数値解の相対速度比較')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    out_path = os.path.join(output_dir, 'particle_collision_validation.png')
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"保存: {out_path}")
    plt.show()

    # 簡易統計（オーバーラップと速度）
    rmse = float(np.sqrt(np.mean(err**2)))
    max_abs = float(np.max(np.abs(err)))
    v_err = v_num - v_theo
    v_rmse = float(np.sqrt(np.mean(v_err**2)))
    v_max_abs = float(np.max(np.abs(v_err)))
    
    print("")
    print("=================================")
    print("オーバーラップと速度の統計")
    print("=================================")
    print(f"Overlap RMSE: {rmse:.6e} m  ( {rmse*1e3:.6f} mm )")
    print(f"Overlap MaxAbsErr: {max_abs:.6e} m  ( {max_abs*1e3:.6f} mm )")
    print(f"Velocity RMSE: {v_rmse:.6e} m/s")
    print(f"Velocity MaxAbsErr: {v_max_abs:.6e} m/s")
    print("=================================")
    print("")


def main():
    setup_japanese_font()
    df = load_csv()
    if df is None:
        return
    
    # 理論値（particle_1dcollision.f90で定義されている値）
    restitution_coeff_e = 0.25
    
    # データの基本情報を表示
    print("=================================")
    print("データ情報")
    print("=================================")
    print(f"サンプル数: {len(df)}")
    print(f"接触開始時刻: {df['time'].iloc[0]:.6e} s")
    print(f"接触終了時刻: {df['time'].iloc[-1]:.6e} s")
    print(f"接触継続時間: {(df['time'].iloc[-1] - df['time'].iloc[0]):.6e} s")
    print("=================================")
    print("")
    
    # 数値解から跳ね返り係数を算出
    e_numeric, v_before, v_after = calculate_restitution_coefficient(df)
    if e_numeric is not None:
        print("=================================")
        print("跳ね返り係数の評価")
        print("=================================")
        print(f"接触前相対速度: v_before = {v_before:.6e} m/s")
        print(f"接触後相対速度: v_after = {v_after:.6e} m/s")
        print(f"数値解 跳ね返り係数: e_numeric = {e_numeric:.6f}")
        print(f"理論値 跳ね返り係数: e_theory = {restitution_coeff_e:.6f}")
        
        # 理論値との相対誤差を計算
        if abs(restitution_coeff_e) > 1e-10:
            relative_error = abs(e_numeric - restitution_coeff_e) / abs(restitution_coeff_e) * 100.0
            print(f"相対誤差: {relative_error:.12f} %")
        else:
            print("相対誤差: 計算不可（理論値がゼロに近い）")
        print("=================================")
        print("")
    else:
        print("跳ね返り係数の計算をスキップしました")
        print("")
    
    # グラフのプロット
    plot_overlap(df)


if __name__ == '__main__':
    main()

