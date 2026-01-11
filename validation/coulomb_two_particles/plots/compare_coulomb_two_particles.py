#!/usr/bin/env python3
"""
数値解と理論解を比較し、誤差統計とプロットを生成するスクリプト。
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# 日本語フォントの設定
# システムにインストールされている日本語フォントを検索
japanese_fonts = ['Noto Sans CJK JP', 'Noto Sans CJK', 'Takao Gothic', 
                  'IPAexGothic', 'IPAPGothic', 'VL PGothic', 'Yu Gothic']
font_found = False

try:
    # 利用可能なフォントを確認
    available_fonts = [f.name for f in fm.fontManager.ttflist]
    
    # 日本語フォントを優先順位で検索
    for font_name in japanese_fonts:
        if font_name in available_fonts:
            plt.rcParams['font.family'] = font_name
            font_found = True
            break
    
    # フォント名の部分一致で検索
    if not font_found:
        for jp_keyword in ['Noto', 'Takao', 'IPA', 'VL', 'Yu']:
            matching = [f for f in available_fonts if jp_keyword in f]
            if matching:
                plt.rcParams['font.family'] = matching[0]
                font_found = True
                break
    
    # フォントが見つからない場合の警告
    if not font_found:
        print("警告: 日本語フォントが見つかりませんでした。")
        print("日本語が正しく表示されない可能性があります。")
        print("利用可能なフォントの一部:", available_fonts[:10])
except Exception as e:
    print(f"警告: 日本語フォントの設定中にエラーが発生しました: {e}")


class CoulombTheory:
    """二粒子クーロン反発の解析解を扱うユーティリティ。"""

    def __init__(self, mass: float, charge: float, k_const: float, r0: float) -> None:
        if charge == 0.0:
            raise ValueError("電荷が0のケースは解析解比較の対象外です。")
        self.mass = mass
        self.charge = charge
        self.k_const = k_const
        self.r0 = r0
        self.tau = math.sqrt((mass * r0 ** 3) / (4.0 * k_const * charge ** 2))
        self.vel_coeff = math.sqrt(k_const * charge ** 2 / mass)

    def time_from_separation(self, r_val: float) -> float:
        if r_val <= self.r0:
            return 0.0
        ratio = r_val / self.r0
        alpha = math.sqrt(ratio)
        beta = math.sqrt(max(ratio - 1.0, 0.0))
        return self.tau * (alpha * beta + math.log(beta + alpha))

    def separation_from_time(self, time_val: float) -> float:
        if time_val <= 0.0:
            return self.r0

        lower = self.r0
        upper = self.r0 * (1.0 + (time_val / self.tau + 1.0) ** 2)
        while self.time_from_separation(upper) < time_val:
            upper *= 1.5

        for _ in range(80):
            mid = 0.5 * (lower + upper)
            if self.time_from_separation(mid) < time_val:
                lower = mid
            else:
                upper = mid
        return 0.5 * (lower + upper)

    def velocity(self, r_val: float) -> float:
        if r_val <= self.r0:
            return 0.0
        diff = max((1.0 / self.r0) - (1.0 / r_val), 0.0)
        return self.vel_coeff * math.sqrt(diff)

    def force_magnitude(self, r_val: float) -> float:
        """粒子間距離 r_val におけるクーロン力の大きさ |F| を返す。"""
        if r_val <= 0.0:
            return 0.0
        return self.k_const * self.charge ** 2 / (r_val ** 2)


def parse_metadata(header_line: str) -> Dict[str, float]:
    metadata: Dict[str, float] = {}
    if not header_line.startswith("#"):
        return metadata
    for token in header_line[1:].split(","):
        if "=" not in token:
            continue
        key, val = token.split("=", 1)
        metadata[key.strip()] = float(val)
    return metadata


def read_dataset(csv_path: Path) -> Tuple[Dict[str, float], np.ndarray]:
    metadata_line = ""
    rows: List[List[float]] = []
    with csv_path.open("r", encoding="utf-8") as fobj:
        for line in fobj:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                metadata_line = stripped
                continue
            if stripped.startswith("time"):
                continue
            rows.append([float(value) for value in stripped.split(",")])
    if not rows:
        raise RuntimeError(f"{csv_path} からデータを読み込めませんでした。")
    return parse_metadata(metadata_line), np.array(rows)


def build_theoretical_profiles(times: np.ndarray, theory: CoulombTheory) -> Tuple[np.ndarray, np.ndarray]:
    separation = np.zeros_like(times)
    velocity = np.zeros_like(times)
    force = np.zeros_like(times)
    for idx, t_val in enumerate(times):
        sep = theory.separation_from_time(float(t_val))
        separation[idx] = sep
        velocity[idx] = theory.velocity(sep)
        force[idx] = theory.force_magnitude(sep)
    return separation, velocity, force


def summarize_errors(numeric: np.ndarray, theory: np.ndarray, label: str) -> None:
    diff = numeric - theory
    max_abs = np.max(np.abs(diff))
    rms = math.sqrt(np.mean(diff ** 2))
    print(f"{label}: max |Δ| = {max_abs:.6e}, RMS = {rms:.6e}")


def plot_profiles(
    times: np.ndarray,
    sep_numeric: np.ndarray,
    sep_theory: np.ndarray,
    vel_numeric: np.ndarray,
    vel_theory: np.ndarray,
    force_numeric: np.ndarray,
    force_theory: np.ndarray,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7, 9))

    # 粒子間距離
    axes[0].plot(times, sep_numeric, label="数値解", color="#1f77b4")
    axes[0].plot(times, sep_theory, "--", label="理論解", color="#ff7f0e")
    axes[0].set_ylabel("粒子間距離 [m]")
    axes[0].legend()
    axes[0].grid(True, linestyle="--", alpha=0.4)

    # 粒子速度
    axes[1].plot(times, vel_numeric, label="数値解", color="#2ca02c")
    axes[1].plot(times, vel_theory, "--", label="理論解", color="#d62728")
    axes[1].set_ylabel("粒子速度 |v| [m/s]")
    axes[1].legend()
    axes[1].grid(True, linestyle="--", alpha=0.4)

    # クーロン力の大きさ
    axes[2].plot(times, force_numeric, label="数値解", color="#9467bd")
    axes[2].plot(times, force_theory, "--", label="理論解", color="#8c564b")
    axes[2].set_xlabel("時間 [s]")
    axes[2].set_ylabel("クーロン力 |F| [N]")
    axes[2].legend()
    axes[2].grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="二粒子クーロン検証: 数値解 vs 理論解")
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("data/coulomb_two_particles.csv"),
        help="シミュレーション出力CSV",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("plots/coulomb_two_particles_comparison.png"),
        help="比較プロットの保存パス",
    )
    args = parser.parse_args()

    metadata, dataset = read_dataset(args.data)
    required_keys = {"mass", "charge", "coulomb_constant", "separation0"}
    missing = required_keys - metadata.keys()
    if missing:
        raise RuntimeError(f"メタデータに {missing} が含まれていません。")

    theory = CoulombTheory(
        mass=metadata["mass"],
        charge=metadata["charge"],
        k_const=metadata["coulomb_constant"],
        r0=metadata["separation0"],
    )

    times = dataset[:, 0]
    x1_num = dataset[:, 1]
    x2_num = dataset[:, 2]
    separation_num = dataset[:, 5]
    vel_num = np.abs(dataset[:, 3])

    separation_theory, vel_theory, force_theory = build_theoretical_profiles(times, theory)
    x1_theory = -0.5 * separation_theory
    x2_theory = 0.5 * separation_theory

    # 数値シミュレーションからのクーロン力（大きさ）を再構成
    # F = k q^2 / r^2 を用いる。r が非常に小さい場合は NaN を入れておく。
    k_const = metadata["coulomb_constant"]
    charge = metadata["charge"]
    with np.errstate(divide="ignore", invalid="ignore"):
        force_num = k_const * charge ** 2 / (separation_num ** 2)
        force_num[~np.isfinite(force_num)] = np.nan

    summarize_errors(separation_num, separation_theory, "距離")
    summarize_errors(x1_num, x1_theory, "x1")
    summarize_errors(x2_num, x2_theory, "x2")
    summarize_errors(vel_num, vel_theory, "速度")
    summarize_errors(force_num, force_theory, "力")

    plot_profiles(
        times,
        separation_num,
        separation_theory,
        vel_num,
        vel_theory,
        force_num,
        force_theory,
        args.plot,
    )
    print(f"プロットを {args.plot} に保存しました。")


if __name__ == "__main__":
    main()

