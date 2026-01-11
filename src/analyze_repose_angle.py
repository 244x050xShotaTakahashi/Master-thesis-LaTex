#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
安息角測定スクリプト

DEMシミュレーションの出力データから粒子堆積物の安息角を測定します。
最終フレームの粒子座標を読み込み、表面粒子を検出して線形回帰により斜面角度を計算します。
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from typing import Union, Tuple, List, Dict, Optional
import sys
from scipy.stats import linregress


def read_particle_data(filename: Union[str, Path] = "data/graph11.d", container_width: Optional[float] = None) -> Dict:
    """
    粒子データを読み込む (graph11.d または CSV)
    
    Args:
        filename: 入力ファイルパス
        container_width: 容器の幅 (CSVの場合に必要、指定がない場合は0.5またはデータから推定)
        
    Returns:
        粒子データを含む辞書 (x, z, r, num_particles, container_width)
    """
    filename = Path(filename)
    
    if not filename.exists():
        raise FileNotFoundError(f"データファイルが見つかりません: {filename}")
    
    if filename.suffix == '.csv':
        return read_csv_data(filename, container_width)
    else:
        return read_graph11_data(filename)


def read_csv_data(filename: Path, container_width: Optional[float] = None) -> Dict:
    """
    CSVファイルから最終フレームの粒子データを読み込む
    """
    try:
        # pandasで読み込み (ヘッダーの空白除去も考慮)
        df = pd.read_csv(filename, skipinitialspace=True)
        
        # 必要なカラムの確認
        required_cols = ['step', 'time', 'x', 'z', 'radius']
        if not all(col in df.columns for col in required_cols):
             raise ValueError(f"CSVファイルに必要なカラムが含まれていません: {required_cols}")

        # 最終ステップのデータを取得
        last_step = df['step'].max()
        last_frame_df = df[df['step'] == last_step]
        
        num_particles = len(last_frame_df)
        if num_particles == 0:
            raise ValueError("有効な粒子データがありません")
            
        time_val = last_frame_df['time'].iloc[0]
        
        # データの抽出
        x_coords = last_frame_df['x'].values
        z_coords = last_frame_df['z'].values
        radii = last_frame_df['radius'].values
        
        # コンテナ幅の決定
        if container_width is None:
            # x座標の最大値から推定 (マージンを少し考慮するか、単に最大値とするか)
            # ここでは安全のためデフォルト0.5を採用しつつ警告を出すか、
            # あるいは引数で指定させる設計にする。
            # 既存コードとの互換性のため、一旦0.5とするが、最大値と比較してチェックする
            max_x = np.max(x_coords) + np.max(radii)
            est_width = 0.5
            if max_x > 0.55: # 明らかに0.5より大きい場合
                est_width = max_x
                print(f"警告: コンテナ幅が指定されていないため、データから推定します: {est_width:.4f} m")
            else:
                print(f"警告: コンテナ幅が指定されていないため、デフォルト値を使用します: {est_width} m")
            container_width = est_width
            
        return {
            'x': x_coords,
            'z': z_coords,
            'r': radii,
            'num_particles': num_particles,
            'time': time_val,
            'container_width': container_width,
            'container_height': 0.0 # CSVには含まれていないため仮定
        }
        
    except Exception as e:
        raise ValueError(f"CSVファイルの読み込みに失敗しました: {e}")


def read_graph11_data(filename: Path) -> Dict:
    """
    graph11.dファイルから最終フレームの粒子データを読み込む
    フォーマット(推測):
      NumParticles Time Width
      Val1 Val2
      x1 z1 r1 x2 z2 r2 ...
      (追加データがある可能性あり)
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()
            
        tokens = content.split()
        if not tokens:
            raise ValueError("ファイルが空です")
            
        # 粒子数はファイルの先頭にあると仮定
        try:
            num_particles_global = int(tokens[0])
        except ValueError:
            raise ValueError("粒子数を読み取れませんでした")

        # スナップショットの開始位置（ヘッダー）を探索
        # ヘッダー構造: [N, Time, Width, Val1, Val2] と仮定 (5トークン)
        # データ構造: [x, z, r] * N
        
        header_size = 5
        data_size = num_particles_global * 3
        block_min_size = header_size + data_size
        
        candidates = []
        i = 0
        
        while i < len(tokens):
            if tokens[i] == str(num_particles_global):
                # 残りトークン数のチェック
                if i + block_min_size <= len(tokens):
                    # 簡易チェック: TimeとWidthが数値か
                    try:
                        _ = float(tokens[i+1])
                        _ = float(tokens[i+2])
                        candidates.append(i)
                        # 次のブロックへスキップ（重複検出を避けるためデータ分進める）
                        i += block_min_size 
                        continue 
                    except ValueError:
                        pass
            i += 1
            
        if not candidates:
            raise ValueError("有効なデータスナップショットが見つかりませんでした")
            
        # 最後のスナップショットを使用
        last_idx = candidates[-1]
        
        time_val = float(tokens[last_idx+1])
        container_width = float(tokens[last_idx+2])
        
        data_start = last_idx + header_size
        data_end = data_start + data_size
        
        # データを抽出して数値変換
        raw_data = [float(x) for x in tokens[data_start : data_end]]
        data_arr = np.array(raw_data)
        
        x = data_arr[0::3]
        z = data_arr[1::3]
        r = data_arr[2::3]
        
        return {
            'x': x,
            'z': z,
            'r': r,
            'num_particles': num_particles_global,
            'time': time_val,
            'container_width': container_width,
            'container_height': 0.0
        }

    except Exception as e:
        print(f"エラー: graph11.dの解析に失敗しました: {e}", file=sys.stderr)
        raise


def detect_surface_particles(x: np.ndarray, z: np.ndarray, r: np.ndarray, 
                             container_width: float,
                             left_margin: float = 0.1,
                             right_margin: float = 0.1,
                             center_margin: float = 0.1,
                             side: str = 'both') -> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]:
    """
    堆積物の表面粒子を検出
    
    Args:
        x: 粒子のx座標配列
        z: 粒子のz座標配列
        r: 粒子の半径配列
        container_width: 容器の幅
        left_margin: 左側のマージン比率 (0.0-1.0)
        right_margin: 右側のマージン比率 (0.0-1.0)
        center_margin: 中央のマージン比率 (0.0-1.0)
        side: 測定対象の斜面 ('left', 'right', 'both')
        
    Returns:
        左斜面の表面粒子(x, z)、右斜面の表面粒子(x, z)
    """
    # 容器の中央
    center_x = container_width / 2.0
    
    # 各x位置での最も高い粒子を検出
    def get_surface_points(region_mask):
        region_x = x[region_mask]
        region_z = z[region_mask]
        region_r = r[region_mask]
        
        if len(region_x) == 0:
            return np.array([]), np.array([])
        
        # x座標でビンを作成
        num_bins = 20
        x_min, x_max = region_x.min(), region_x.max()
        if x_min == x_max:
             return np.array([x_min]), np.array([region_z[0]])

        bins = np.linspace(x_min, x_max, num_bins)
        
        surface_x = []
        surface_z = []
        
        for i in range(len(bins) - 1):
            bin_mask = (region_x >= bins[i]) & (region_x < bins[i + 1])
            if bin_mask.sum() > 0:
                # このビン内で最も高い粒子（粒子上端）
                z_top = region_z[bin_mask] + region_r[bin_mask]
                max_idx = np.argmax(z_top)
                indices = np.where(bin_mask)[0]
                particle_idx = indices[max_idx]
                
                surface_x.append(region_x[particle_idx])
                surface_z.append(region_z[particle_idx])
        
        return np.array(surface_x), np.array(surface_z)
    
    left_surface_x, left_surface_z = np.array([]), np.array([])
    right_surface_x, right_surface_z = np.array([]), np.array([])
    
    if side in ['left', 'both']:
        # 左斜面: 左壁側から中央に向かって上昇する斜面
        left_region = (x < center_x - center_margin * container_width) & (x > left_margin * container_width)
        left_surface_x, left_surface_z = get_surface_points(left_region)
    
    if side in ['right', 'both']:
        # 右斜面: 右壁側から中央に向かって上昇する斜面（壁引き抜き後の主斜面）
        if side == 'right':
            # 片側モード: 中央マージンなし、右壁側から全体を対象
            right_region = (x > left_margin * container_width) & (x < (1 - right_margin) * container_width)
        else:
            # 両側モード: 中央マージンあり
            right_region = (x > center_x + center_margin * container_width) & (x < (1 - right_margin) * container_width)
        right_surface_x, right_surface_z = get_surface_points(right_region)
    
    # 片側モード（left）でも同様に調整
    if side == 'left':
        left_region = (x > left_margin * container_width) & (x < (1 - right_margin) * container_width)
        left_surface_x, left_surface_z = get_surface_points(left_region)
    
    return (left_surface_x, left_surface_z), (right_surface_x, right_surface_z)


def calculate_repose_angle(surface_x: np.ndarray, surface_z: np.ndarray, 
                           side: str = 'left') -> Tuple[float, float, float, float]:
    """
    表面粒子から安息角を計算
    
    Args:
        surface_x: 表面粒子のx座標
        surface_z: 表面粒子のz座標
        side: 'left' または 'right'
        
    Returns:
        angle_deg: 安息角（度）
        slope: 斜面の傾き
        intercept: 切片
        r_squared: 決定係数
    """
    if len(surface_x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    
    # 線形回帰
    slope, intercept, r_value, p_value, std_err = linregress(surface_x, surface_z)
    
    # 傾きから角度を計算
    angle_rad = np.arctan(abs(slope))
    angle_deg = np.degrees(angle_rad)
    
    # 左斜面の場合は正の傾き、右斜面の場合は負の傾きになるはず
    if side == 'right':
        angle_rad = np.arctan(abs(slope))
        angle_deg = np.degrees(angle_rad)
    
    r_squared = r_value ** 2
    
    return angle_deg, slope, intercept, r_squared


def plot_repose_angle(data: Dict, left_surface: Tuple, right_surface: Tuple,
                     left_fit: Tuple, right_fit: Tuple, output_file: str = None,
                     side: str = 'both'):
    """
    安息角の測定結果を可視化
    
    Args:
        data: 粒子データ
        left_surface: 左斜面の表面粒子座標
        right_surface: 右斜面の表面粒子座標
        left_fit: 左斜面のフィッティング結果 (angle, slope, intercept, r2)
        right_fit: 右斜面のフィッティング結果 (angle, slope, intercept, r2)
        output_file: 出力ファイル名
        side: 測定対象の斜面 ('left', 'right', 'both')
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 全粒子を描画
    x, z, r = data['x'], data['z'], data['r']
    for i in range(len(x)):
        circle = plt.Circle((x[i], z[i]), r[i], color='lightgray', ec='gray', linewidth=0.5)
        ax.add_patch(circle)
    
    # 表面粒子を強調
    left_x, left_z = left_surface
    right_x, right_z = right_surface
    
    if len(left_x) > 0 and side in ['left', 'both']:
        label = '表面粒子' if side == 'left' else '左斜面表面粒子'
        ax.scatter(left_x, left_z, c='red', s=50, marker='o', label=label, zorder=5)
    
    if len(right_x) > 0 and side in ['right', 'both']:
        label = '表面粒子' if side == 'right' else '右斜面表面粒子'
        ax.scatter(right_x, right_z, c='blue', s=50, marker='o', label=label, zorder=5)
    
    # フィッティング直線を描画
    left_angle, left_slope, left_intercept, left_r2 = left_fit
    right_angle, right_slope, right_intercept, right_r2 = right_fit
    
    if not np.isnan(left_angle) and len(left_x) > 0 and side in ['left', 'both']:
        x_fit = np.array([left_x.min(), left_x.max()])
        z_fit = left_slope * x_fit + left_intercept
        label = f'安息角: {left_angle:.1f}° (R²={left_r2:.3f})' if side == 'left' else f'左斜面: {left_angle:.1f}° (R²={left_r2:.3f})'
        ax.plot(x_fit, z_fit, 'r--', linewidth=2, label=label)
    
    if not np.isnan(right_angle) and len(right_x) > 0 and side in ['right', 'both']:
        x_fit = np.array([right_x.min(), right_x.max()])
        z_fit = right_slope * x_fit + right_intercept
        label = f'安息角: {right_angle:.1f}° (R²={right_r2:.3f})' if side == 'right' else f'右斜面: {right_angle:.1f}° (R²={right_r2:.3f})'
        ax.plot(x_fit, z_fit, 'b--', linewidth=2, label=label)
    
    # 安息角を表示
    if side == 'both':
        angles = [a for a in [left_angle, right_angle] if not np.isnan(a)]
        if angles:
            avg_angle = np.mean(angles)
            ax.text(0.5, 0.98, f'平均安息角: {avg_angle:.2f}°', 
                    transform=ax.transAxes, fontsize=14, weight='bold',
                    ha='center', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    elif side == 'left' and not np.isnan(left_angle):
        ax.text(0.5, 0.98, f'安息角: {left_angle:.2f}°', 
                transform=ax.transAxes, fontsize=14, weight='bold',
                ha='center', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    elif side == 'right' and not np.isnan(right_angle):
        ax.text(0.5, 0.98, f'安息角: {right_angle:.2f}°', 
                transform=ax.transAxes, fontsize=14, weight='bold',
                ha='center', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel('x [m]', fontsize=12)
    ax.set_ylabel('z [m]', fontsize=12)
    
    # タイトルを測定モードに応じて変更
    if side == 'left':
        ax.set_title('安息角測定結果 (左斜面)', fontsize=14, weight='bold')
    elif side == 'right':
        ax.set_title('安息角測定結果 (右斜面 - 壁引き抜き後)', fontsize=14, weight='bold')
    else:
        ax.set_title('安息角測定結果', fontsize=14, weight='bold')
    
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    # 容器の底を描画
    ax.axhline(y=0, color='black', linewidth=2)
    ax.axvline(x=0, color='black', linewidth=2)
    # 壁引き抜きモードでは右壁を点線で表示
    if side == 'right':
        ax.axvline(x=data['container_width'], color='gray', linewidth=1, linestyle='--', alpha=0.5)
    else:
        ax.axvline(x=data['container_width'], color='black', linewidth=2)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"プロットを保存しました: {output_file}")
    
    plt.close()


def analyze_repose_angle(data_file: str = "data/graph11.d", 
                        output_dir: str = "results",
                        case_name: str = "repose_angle",
                        x_range: Optional[Tuple[float, float]] = None,
                        z_range: Optional[Tuple[float, float]] = None,
                        container_width: Optional[float] = None,
                        margins: Tuple[float, float, float] = (0.1, 0.1, 0.1),
                        side: str = 'both') -> Dict:
    """
    安息角を測定してCSVに保存
    
    Args:
        data_file: 入力データファイル
        output_dir: 出力ディレクトリ
        case_name: ケース名
        x_range: x座標の範囲 (min, max)
        z_range: z座標の範囲 (min, max)
        container_width: 容器の幅 (CSV入力時用)
        margins: 解析除外マージン (left, right, center)
        side: 測定対象の斜面 ('left', 'right', 'both')
        
    Returns:
        測定結果を含む辞書
    """
    # データ読み込み
    print(f"データを読み込んでいます: {data_file}")
    data = read_particle_data(data_file, container_width)
    print(f"粒子数: {data['num_particles']}, 時刻: {data['time']:.4f} s")
    print(f"測定モード: {side} (斜面)")
    
    # 表面粒子を検出
    print("表面粒子を検出中...")
    
    mask = np.ones(len(data['x']), dtype=bool)
    
    if x_range:
        x_min_str = x_range[0] if x_range[0] is not None else '-inf'
        x_max_str = x_range[1] if x_range[1] is not None else 'inf'
        print(f"x座標範囲: {x_min_str} <= x <= {x_max_str}")
        
        if x_range[0] is not None:
            mask &= (data['x'] >= x_range[0])
        if x_range[1] is not None:
            mask &= (data['x'] <= x_range[1])
        
    if z_range:
        z_min_str = z_range[0] if z_range[0] is not None else '-inf'
        z_max_str = z_range[1] if z_range[1] is not None else 'inf'
        print(f"z座標範囲: {z_min_str} <= z <= {z_max_str}")
        
        if z_range[0] is not None:
            mask &= (data['z'] >= z_range[0])
        if z_range[1] is not None:
            mask &= (data['z'] <= z_range[1])

    filtered_x = data['x'][mask]
    filtered_z = data['z'][mask]
    filtered_r = data['r'][mask]
    
    print(f"解析対象粒子数: {len(filtered_x)} (全粒子数: {len(data['x'])})")
    
    # マージン設定
    left_m, right_m, center_m = margins
    
    # 範囲指定がある場合はマージンを強制的に0にする
    if x_range is not None:
        print("x座標範囲が指定されているため、自動除外マージンを無効化します (0.0)")
        left_m, right_m, center_m = 0.0, 0.0, 0.0
    else:
        print(f"解析除外マージン: Left={left_m}, Right={right_m}, Center={center_m}")
        
    left_surface, right_surface = detect_surface_particles(
        filtered_x, filtered_z, filtered_r, data['container_width'],
        left_margin=left_m, right_margin=right_m, center_margin=center_m,
        side=side
    )
    
    left_x, left_z = left_surface
    right_x, right_z = right_surface
    
    if side == 'both':
        print(f"左斜面表面粒子数: {len(left_x)}, 右斜面表面粒子数: {len(right_x)}")
    elif side == 'left':
        print(f"左斜面表面粒子数: {len(left_x)}")
    else:
        print(f"右斜面表面粒子数: {len(right_x)}")
    
    # 安息角を計算
    print("安息角を計算中...")
    left_fit = calculate_repose_angle(left_x, left_z, side='left')
    right_fit = calculate_repose_angle(right_x, right_z, side='right')
    
    left_angle, left_slope, left_intercept, left_r2 = left_fit
    right_angle, right_slope, right_intercept, right_r2 = right_fit
    
    # 結果表示
    if side == 'both':
        # 平均安息角
        angles = [a for a in [left_angle, right_angle] if not np.isnan(a)]
        avg_angle = np.mean(angles) if angles else np.nan
        std_angle = np.std(angles) if len(angles) > 1 else 0.0
        
        print(f"左斜面安息角: {left_angle:.2f}° (R²={left_r2:.3f})")
        print(f"右斜面安息角: {right_angle:.2f}° (R²={right_r2:.3f})")
        print(f"平均安息角: {avg_angle:.2f}° ± {std_angle:.2f}°")
    elif side == 'left':
        avg_angle = left_angle if not np.isnan(left_angle) else np.nan
        std_angle = 0.0
        print(f"安息角 (左斜面): {left_angle:.2f}° (R²={left_r2:.3f})")
    else:  # side == 'right'
        avg_angle = right_angle if not np.isnan(right_angle) else np.nan
        std_angle = 0.0
        print(f"安息角 (右斜面): {right_angle:.2f}° (R²={right_r2:.3f})")
    
    # 出力ディレクトリを作成
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # プロット作成
    plot_file = output_path / f"{case_name}_repose_angle.png"
    plot_repose_angle(data, left_surface, right_surface, left_fit, right_fit, 
                     output_file=str(plot_file), side=side)
    
    # CSVに保存
    csv_file = output_path / f"{case_name}_results.csv"
    with open(csv_file, 'w') as f:
        f.write("parameter,value\n")
        f.write(f"case_name,{case_name}\n")
        f.write(f"num_particles,{data['num_particles']}\n")
        f.write(f"simulation_time,{data['time']}\n")
        f.write(f"measurement_side,{side}\n")
        if side in ['left', 'both']:
            f.write(f"left_angle_deg,{left_angle}\n")
            f.write(f"left_r_squared,{left_r2}\n")
        if side in ['right', 'both']:
            f.write(f"right_angle_deg,{right_angle}\n")
            f.write(f"right_r_squared,{right_r2}\n")
        f.write(f"average_angle_deg,{avg_angle}\n")
        f.write(f"std_angle_deg,{std_angle}\n")
        if x_range:
            f.write(f"x_range_min,{x_range[0]}\n")
            f.write(f"x_range_max,{x_range[1]}\n")
        if z_range:
            f.write(f"z_range_min,{z_range[0]}\n")
            f.write(f"z_range_max,{z_range[1]}\n")
    
    print(f"結果を保存しました: {csv_file}")
    
    return {
        'left_angle': left_angle,
        'right_angle': right_angle,
        'average_angle': avg_angle,
        'std_angle': std_angle,
        'left_r2': left_r2,
        'right_r2': right_r2,
        'side': side
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='DEMシミュレーションから安息角を測定')
    parser.add_argument('--data', '-d', default='data/graph11.d', 
                       help='入力データファイル (デフォルト: data/graph11.d)')
    parser.add_argument('--output', '-o', default='results',
                       help='出力ディレクトリ (デフォルト: results)')
    parser.add_argument('--name', '-n', default='repose_angle',
                       help='ケース名 (デフォルト: repose_angle)')
    
    parser.add_argument('--x-min', type=float, help='x座標の最小値')
    parser.add_argument('--x-max', type=float, help='x座標の最大値')
    parser.add_argument('--z-min', type=float, help='z座標の最小値')
    parser.add_argument('--z-max', type=float, help='z座標の最大値')
    parser.add_argument('--width', '-w', type=float, default=None,
                       help='容器の幅 (CSV入力の場合に推奨、デフォルト: 0.5)')
    parser.add_argument('--margin-center', type=float, default=0.1,
                       help='中央の解析除外マージン比率 (デフォルト: 0.1)')
    parser.add_argument('--margin-side', type=float, default=0.1,
                       help='左右両端の解析除外マージン比率 (デフォルト: 0.1)')
    parser.add_argument('--side', '-s', choices=['left', 'right', 'both'], default='both',
                       help='測定対象の斜面 (left=左壁引き抜き後, right=右壁引き抜き後, both=両斜面, デフォルト: both)')

    args = parser.parse_args()
    
    x_range = None
    if args.x_min is not None or args.x_max is not None:
        x_range = (args.x_min, args.x_max)
        
    z_range = None
    if args.z_min is not None or args.z_max is not None:
        z_range = (args.z_min, args.z_max)
    
    margins = (args.margin_side, args.margin_side, args.margin_center)
    
    try:
        results = analyze_repose_angle(args.data, args.output, args.name, x_range, z_range, 
                                       args.width, margins, args.side)
        print("\n=== 測定完了 ===")
        if args.side == 'both':
            print(f"平均安息角: {results['average_angle']:.2f}°")
        elif args.side == 'left':
            print(f"安息角 (左斜面): {results['left_angle']:.2f}°")
        else:
            print(f"安息角 (右斜面): {results['right_angle']:.2f}°")
    except Exception as e:
        print(f"エラーが発生しました: {e}", file=sys.stderr)
        sys.exit(1)
