import os
import csv
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_velocity_data(path: str) -> Tuple[List[float], List[float]]:
    """
    CSVファイルから時間と数値解の速度データを読み込む
    Returns: (time, velocity_numeric)
    """
    t, vn = [], []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            # time,x_num,v_num,x_theory,v_theory
            if len(row) < 3:
                continue
            t.append(float(row[0]))
            vn.append(float(row[2]))  # v_num
    return t, vn


def plot_velocity(t: List[float], vn: List[float], output_path: str) -> None:
    """数値解の速度をプロットして保存"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(t, vn, label="数値解の速度 (v_num)", lw=1.5, color='#2E86AB')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    ax.set_xlabel("時間 [s]", fontsize=12)
    ax.set_ylabel("速度 [m/s]", fontsize=12)
    ax.set_title("減衰振動子の速度変化（数値解）", fontsize=14, fontweight='bold')
    ax.grid(True, ls=":", lw=0.5, alpha=0.7)
    ax.legend(fontsize=11)
    ax.set_xlim(0, 0.01)
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=200)
    print(f"グラフを保存しました: {output_path}")
    plt.close(fig)


def print_velocity_stats(t: List[float], vn: List[float]) -> None:
    """速度データの統計情報を表示"""
    print(f"\n=== 速度データの統計情報 ===")
    print(f"データ点数: {len(vn)}")
    print(f"時間範囲: {min(t):.6e} ~ {max(t):.6e} s")
    print(f"速度の最大値: {max(vn):.6e} m/s")
    print(f"速度の最小値: {min(vn):.6e} m/s")
    print(f"\n最初の10個の速度値:")
    for i in range(min(10, len(vn))):
        print(f"  t={t[i]:8.4f} s  →  v={vn[i]:12.6e} m/s")


def main() -> None:
    csv_path = os.path.join("PEM", "data", "damped_oscillator.csv")
    out_png = os.path.join("PEM", "data", "velocity_numeric.png")
    
    if not os.path.exists(csv_path):
        raise SystemExit(f"CSVが見つかりません: {csv_path}")
    
    # データ読み込み
    t, vn = read_velocity_data(csv_path)
    
    # 統計情報を表示
    print_velocity_stats(t, vn)
    
    # グラフをプロット
    plot_velocity(t, vn, out_png)


if __name__ == "__main__":
    main()

