#!/usr/bin/env python3
"""
粒子スナップショットプロッター

DEMシミュレーション結果（particles.csv）を読み込み、
指定した時間またはステップの粒子をプロットします。
粒子は半径で色分けされ、壁（コンテナ境界・斜面壁）も表示されます。

使用例:
    python scripts/plot_snapshot.py --file data/particles.csv --time 0.5 --output snapshot.png
    python scripts/plot_snapshot.py --step 100000 --width 0.5
    python scripts/plot_snapshot.py --time 0.01  # 画面表示
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
import matplotlib.cm as cm


def load_particles(filepath: str) -> pd.DataFrame:
    """
    particles.csvファイルを読み込む
    
    Parameters
    ----------
    filepath : str
        particles.csvファイルのパス
    
    Returns
    -------
    pd.DataFrame
        粒子データ
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"ファイルが見つかりません: {filepath}")
    
    df = pd.read_csv(filepath, skipinitialspace=True)
    # カラム名の空白を削除
    df.columns = df.columns.str.strip()
    return df


def load_walls(filepath: str) -> list:
    """
    walls.datファイルを読み込む
    
    Parameters
    ----------
    filepath : str
        walls.datファイルのパス
    
    Returns
    -------
    list
        壁データのリスト [(x1, z1, x2, z2), ...]
    """
    walls = []
    if not os.path.exists(filepath):
        return walls
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # 空行とコメント行をスキップ
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x1, z1, x2, z2 = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                    walls.append((x1, z1, x2, z2))
                except ValueError:
                    continue
    
    return walls


def get_snapshot_data(df: pd.DataFrame, time: float = None, step: int = None) -> pd.DataFrame:
    """
    指定した時間またはステップのスナップショットデータを取得
    
    Parameters
    ----------
    df : pd.DataFrame
        全粒子データ
    time : float, optional
        プロットする時間 [秒]
    step : int, optional
        プロットするステップ番号
    
    Returns
    -------
    pd.DataFrame
        スナップショットデータ
    """
    if step is not None:
        # ステップで抽出
        snapshot = df[df['step'] == step]
        if snapshot.empty:
            available_steps = df['step'].unique()
            # 最も近いステップを探す
            closest_step = min(available_steps, key=lambda x: abs(x - step))
            print(f"警告: ステップ {step} が見つかりません。最も近いステップ {closest_step} を使用します。")
            snapshot = df[df['step'] == closest_step]
    elif time is not None:
        # 時間で抽出（最も近い時間を選択）
        available_times = df['time'].unique()
        closest_time = min(available_times, key=lambda x: abs(x - time))
        if abs(closest_time - time) > 1e-10:
            print(f"情報: 時間 {time:.6f}s に最も近い時間 {closest_time:.6f}s を使用します。")
        snapshot = df[df['time'] == closest_time]
    else:
        # 最初のステップを使用
        first_step = df['step'].min()
        snapshot = df[df['step'] == first_step]
        print(f"情報: 時間/ステップが指定されていません。最初のステップ {first_step} を使用します。")
    
    return snapshot


def plot_snapshot(snapshot: pd.DataFrame, 
                  walls: list,
                  container_width: float,
                  container_height: float,
                  output_path: str = None,
                  title: str = None,
                  color_by: str = 'radius',
                  charge_mode: str = 'actual',
                  dpi: int = 150,
                  x_min: float = None,
                  x_max: float = None,
                  z_min: float = None,
                  z_max: float = None):
    """
    スナップショットをプロット
    
    Parameters
    ----------
    snapshot : pd.DataFrame
        スナップショットデータ
    walls : list
        斜面壁データ
    container_width : float
        コンテナ幅
    container_height : float
        コンテナ高さ (0なら上壁なし)
    output_path : str, optional
        出力ファイルパス (Noneなら画面表示)
    title : str, optional
        プロットタイトル
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 粒子データを取得
    x = snapshot['x'].values
    z = snapshot['z'].values
    radii = snapshot['radius'].values
    
    # 着色に使うスカラー値を決定
    if color_by == 'radius':
        color_values = radii
        cmap = cm.viridis
        colorbar_label = 'Radius [m]'
        vmin = float(np.nanmin(color_values))
        vmax = float(np.nanmax(color_values))
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            raise ValueError("エラー: radius の値が不正です (NaN/Inf)。")
        if vmin == vmax:
            eps = max(abs(vmin) * 1e-12, 1e-12)
            vmin, vmax = vmin - eps, vmax + eps
        norm = plt.Normalize(vmin, vmax)
    elif color_by == 'charge':
        if 'charge' not in snapshot.columns:
            raise ValueError("エラー: particles.csv に 'charge' 列がありません。--color-by radius を使うか、出力に charge を含めてください。")
        if charge_mode == 'actual':
            color_values = snapshot['charge'].values.astype(float)
        elif charge_mode == 'zero':
            # 発表用など「帯電なし」と誤解されないよう、可視化上は0Cに強制する
            color_values = np.zeros(len(snapshot), dtype=float)
        else:
            raise ValueError(f"エラー: charge_mode は 'actual' または 'zero' を指定してください (指定値: {charge_mode})")
        cmap = cm.coolwarm
        colorbar_label = 'Charge [C]'
        vmin = float(np.nanmin(color_values))
        vmax = float(np.nanmax(color_values))
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            raise ValueError("エラー: charge の値が不正です (NaN/Inf)。")
        v = max(abs(vmin), abs(vmax))
        if v == 0.0:
            # 全て0の場合でも中心色(白)で描けるように、適当な対称レンジを与える
            v = 1.0
            if charge_mode != 'zero':
                print("警告: charge が全て0のため、発散カラーマップの範囲を [-1, +1] に設定します。")
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-v, vmax=+v)
    else:
        raise ValueError(f"エラー: --color-by は 'radius' または 'charge' を指定してください (指定値: {color_by})")

    colors = cmap(norm(color_values))
    
    # 粒子を円としてプロット
    patches = []
    for xi, zi, ri in zip(x, z, radii):
        circle = Circle((xi, zi), ri)
        patches.append(circle)
    
    collection = PatchCollection(patches, facecolors=colors, edgecolors='black', linewidths=0.5, alpha=0.8)
    ax.add_collection(collection)
    
    # カラーバーを追加
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label=colorbar_label)
    
    # コンテナ壁を描画
    wall_color = 'gray'
    wall_lw = 2
    
    # 左壁 (x=0)
    ax.axvline(x=0, color=wall_color, linewidth=wall_lw, linestyle='-')
    
    # 右壁 (x=container_width)
    ax.axvline(x=container_width, color=wall_color, linewidth=wall_lw, linestyle='-')
    
    # 下壁 (z=0)
    ax.axhline(y=0, color=wall_color, linewidth=wall_lw, linestyle='-')
    
    # 上壁 (container_height > 0 の場合のみ)
    if container_height > 0:
        ax.axhline(y=container_height, color=wall_color, linewidth=wall_lw, linestyle='-')
    
    # 斜面壁を描画
    for wall in walls:
        x1, z1, x2, z2 = wall
        ax.plot([x1, x2], [z1, z2], color='brown', linewidth=wall_lw + 1, linestyle='-', label='斜面壁')
    
    # 軸の設定
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    
    # プロット範囲の設定
    if x_min is not None or x_max is not None or z_min is not None or z_max is not None:
        # ユーザー指定がある場合はそれを優先（指定されない側は自動）
        if x_min is None:
            x_min = float(np.nanmin(x - radii))
        if x_max is None:
            x_max = float(np.nanmax(x + radii))
        if z_min is None:
            z_min = float(np.nanmin(z - radii))
        if z_max is None:
            z_max = float(np.nanmax(z + radii))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(z_min, z_max)
    else:
        # 全体表示（従来の自動範囲）
        x_margin = max(radii.max() * 2, container_width * 0.05)
        z_margin = max(radii.max() * 2, 0.02)
        
        ax.set_xlim(-x_margin, container_width + x_margin)
        
        z_top = z.max() + radii.max() + z_margin
        if container_height > 0:
            z_top = max(z_top, container_height + z_margin)
        ax.set_ylim(-z_margin, z_top)
    
    # タイトル
    if title:
        ax.set_title(title)
    else:
        step_val = snapshot['step'].iloc[0]
        time_val = snapshot['time'].iloc[0]
        n_particles = len(snapshot)
        ax.set_title(f'Step {step_val}, Time = {time_val:.6f} s, N = {n_particles} particles')
    
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"スナップショットを保存しました: {output_path}")
    else:
        plt.show()
    
    plt.close(fig)


def list_available_steps(df: pd.DataFrame, max_display: int = 20):
    """
    利用可能なステップと時間を一覧表示
    
    Parameters
    ----------
    df : pd.DataFrame
        粒子データ
    max_display : int
        表示する最大数
    """
    steps = df.groupby('step')['time'].first().reset_index()
    steps = steps.sort_values('step')
    
    print(f"\n利用可能なステップ (全{len(steps)}個):")
    print("-" * 40)
    print(f"{'ステップ':>12} {'時間 [s]':>15}")
    print("-" * 40)
    
    if len(steps) <= max_display:
        for _, row in steps.iterrows():
            print(f"{int(row['step']):>12} {row['time']:>15.6e}")
    else:
        # 最初と最後の数個を表示
        half = max_display // 2
        for _, row in steps.head(half).iterrows():
            print(f"{int(row['step']):>12} {row['time']:>15.6e}")
        print(f"{'...':>12} {'...':>15}")
        for _, row in steps.tail(half).iterrows():
            print(f"{int(row['step']):>12} {row['time']:>15.6e}")
    
    print("-" * 40)


def main():
    parser = argparse.ArgumentParser(
        description='DEMシミュレーション結果の粒子スナップショットをプロット',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  python scripts/plot_snapshot.py --file data/particles.csv --time 0.5
  python scripts/plot_snapshot.py --step 100000 --width 0.5 --output snapshot.png
  python scripts/plot_snapshot.py --time 0.01 --color-by charge --output snapshot_charge.png
  python scripts/plot_snapshot.py --step 8273861 --color-by charge --dpi 400 --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 --output snapshot_charge_zoom.png
  python scripts/plot_snapshot.py --list  # 利用可能なステップを一覧表示
        """
    )
    
    parser.add_argument('--file', '-f', type=str, default='data/particles.csv',
                        help='particles.csvファイルのパス (デフォルト: data/particles.csv)')
    parser.add_argument('--time', '-t', type=float, default=None,
                        help='プロットする時間 [秒]')
    parser.add_argument('--step', '-s', type=int, default=None,
                        help='プロットするステップ番号')
    parser.add_argument('--walls', '-w', type=str, default='inputs/walls.dat',
                        help='walls.datファイルのパス (デフォルト: inputs/walls.dat)')
    parser.add_argument('--width', type=float, default=0.5,
                        help='コンテナ幅 [m] (デフォルト: 0.5)')
    parser.add_argument('--height', type=float, default=0.0,
                        help='コンテナ高さ [m] (デフォルト: 0.0, 0なら上壁なし)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='出力画像ファイルパス (指定なしなら画面表示)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='保存画像の解像度[dpi] (デフォルト: 150)')
    parser.add_argument('--color-by', type=str, default='radius', choices=['radius', 'charge'],
                        help='粒子の着色に使う量 (radius|charge, デフォルト: radius)')
    parser.add_argument('--charge-mode', type=str, default='actual', choices=['actual', 'zero'],
                        help="color-by=charge時の可視化モード (actual=実データ, zero=0C固定; デフォルト: actual)")
    parser.add_argument('--x-min', type=float, default=None, help='表示範囲のx最小値[m] (指定すると拡大表示)')
    parser.add_argument('--x-max', type=float, default=None, help='表示範囲のx最大値[m]')
    parser.add_argument('--z-min', type=float, default=None, help='表示範囲のz最小値[m]')
    parser.add_argument('--z-max', type=float, default=None, help='表示範囲のz最大値[m]')
    parser.add_argument('--list', '-l', action='store_true',
                        help='利用可能なステップを一覧表示')
    parser.add_argument('--title', type=str, default=None,
                        help='プロットのタイトル')
    
    args = parser.parse_args()
    
    # データ読み込み
    try:
        df = load_particles(args.file)
    except FileNotFoundError as e:
        print(f"エラー: {e}", file=sys.stderr)
        sys.exit(1)
    
    # ステップ一覧表示モード
    if args.list:
        list_available_steps(df)
        sys.exit(0)
    
    # 壁データの読み込み
    walls = load_walls(args.walls)
    if walls:
        print(f"斜面壁を読み込みました: {len(walls)} 本")
    
    # スナップショットデータの取得
    snapshot = get_snapshot_data(df, time=args.time, step=args.step)
    
    if snapshot.empty:
        print("エラー: スナップショットデータが空です。", file=sys.stderr)
        sys.exit(1)
    
    print(f"スナップショット: ステップ {snapshot['step'].iloc[0]}, "
          f"時間 {snapshot['time'].iloc[0]:.6e} s, "
          f"粒子数 {len(snapshot)}")
    
    # プロット
    plot_snapshot(
        snapshot=snapshot,
        walls=walls,
        container_width=args.width,
        container_height=args.height,
        output_path=args.output,
        title=args.title,
        color_by=args.color_by,
        charge_mode=args.charge_mode,
        dpi=args.dpi,
        x_min=args.x_min,
        x_max=args.x_max,
        z_min=args.z_min,
        z_max=args.z_max,
    )


if __name__ == '__main__':
    main()

