import os
import csv
import math
from typing import List, Tuple, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path: str) -> Tuple[List[float], List[float], List[float], List[float], List[float]]:
    t, xn, vn, xt, vt = [], [], [], [], []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            # time,x_num,v_num,x_theory,v_theory
            if len(row) < 5:
                continue
            t.append(float(row[0]))
            xn.append(float(row[1]))
            vn.append(float(row[2]))
            xt.append(float(row[3]))
            vt.append(float(row[4]))
    return t, xn, vn, xt, vt


def calculate_restitution_from_velocities(xn: List[float], vn: List[float]) -> Tuple[float, float, float]:
    """
    ゼロクロス前後の速度から跳ね返り係数を計算
    Returns: (e, v_before, v_after)
    """
    crossings = find_zero_crossings(xn)
    if len(crossings) < 2:
        raise ValueError("ゼロクロスが2回見つかりませんでした")
    
    # 1回目のゼロクロス直前の速度（衝突前）
    idx_before = crossings[0]
    v_before = vn[idx_before]
    
    # 2回目のゼロクロス直後の速度（衝突後）
    idx_after = crossings[1] + 1 if (crossings[1] + 1) < len(vn) else crossings[1]
    v_after = vn[idx_after]
    
    # 跳ね返り係数 e = |v_after| / |v_before|
    if abs(v_before) > 0.0:
        e = abs(v_after) / abs(v_before)
    else:
        e = float('nan')
    
    return e, v_before, v_after


def find_zero_crossings(x: List[float]) -> List[int]:
    """x配列から符号変化（ゼロクロス）のインデックスを検出"""
    crossings = []
    for i in range(len(x) - 1):
        # x[i]とx[i+1]の符号が異なる場合
        if x[i] * x[i+1] < 0:
            crossings.append(i)
    return crossings


def extract_half_period_range(x: List[float]) -> Optional[Tuple[int, int]]:
    """
    1回目と2回目のゼロクロス間のインデックス範囲を返す
    Returns: (start_idx, end_idx) or None
    """
    crossings = find_zero_crossings(x)
    if len(crossings) < 2:
        return None
    return (crossings[0], crossings[1] + 1)  # +1 to include the end point


def plot_data(t: List[float], xn: List[float], vn: List[float], 
              xt: List[float], vt: List[float], output_path: str, title_suffix: str = "") -> None:
    """データをプロットして保存"""
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    ax[0].plot(t, xn, label="x_num", lw=1.2)
    ax[0].plot(t, xt, label="x_theory", lw=1.2, ls="--")
    ylabel = "Displacement x [m]"
    if title_suffix:
        ylabel += f" {title_suffix}"
    ax[0].set_ylabel(ylabel)
    ax[0].grid(True, ls=":", lw=0.5)
    ax[0].legend()

    ax[1].plot(t, vn, label="v_num", lw=1.0)
    ax[1].plot(t, vt, label="v_theory", lw=1.0, ls="--")
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Velocity v [m/s]")
    ax[1].grid(True, ls=":", lw=0.5)
    ax[1].legend()

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=200)
    print(f"保存しました: {output_path}")
    plt.close(fig)


def main() -> None:
    csv_path = os.path.join("PEM", "data", "damped_oscillator_rk4.csv")
    out_png_full = os.path.join("PEM", "data", "damped_oscillator_rk4.png")
    out_png_half = os.path.join("PEM", "data", "damped_oscillator_rk4_half_period.png")

    if not os.path.exists(csv_path):
        raise SystemExit(f"CSVが見つかりません: {csv_path}. 先に Fortran 実行を行ってください。")

    t, xn, vn, xt, vt = read_csv(csv_path)

    # 理論値 (damped_oscillator_rk4.f90で定義されている値)
    restitution_coeff_e = 0.25
    
    # 数値解から跳ね返り係数を算出
    try:
        e_numeric, v_before, v_after = calculate_restitution_from_velocities(xn, vn)
        print(f"跳ね返り係数 e = {e_numeric:.6f}  (v_before={v_before:.6e}, v_after={v_after:.6e})")
        
        # 理論値との相対誤差を計算
        relative_error = abs(e_numeric - restitution_coeff_e) / abs(restitution_coeff_e) * 100.0
        print(f"相対誤差 = {relative_error:.6f}%  (理論値 e_theory = {restitution_coeff_e})")
    except ValueError as e:
        print(f"跳ね返り係数の計算をスキップ: {e}")

    # 全体のプロット
    plot_data(t, xn, vn, xt, vt, out_png_full)

    # 半周期のプロット
    half_range = extract_half_period_range(xn)
    if half_range is not None:
        start_idx, end_idx = half_range
        t_half = t[start_idx:end_idx]
        xn_half = xn[start_idx:end_idx]
        vn_half = vn[start_idx:end_idx]
        xt_half = xt[start_idx:end_idx]
        vt_half = vt[start_idx:end_idx]
        
        plot_data(t_half, xn_half, vn_half, xt_half, vt_half, out_png_half, "(Half Period, RK4)")
        print(f"半周期の範囲: t = {t_half[0]:.6e} ~ {t_half[-1]:.6e} s")
    else:
        print("警告: ゼロクロスが2回見つかりませんでした。半周期プロットをスキップします。")


if __name__ == "__main__":
    main()




