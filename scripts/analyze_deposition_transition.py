#!/usr/bin/env python3
"""
堆積状態の遷移を「定量指標」で可視化するスクリプト

- particles.csv を読み込み、ステップをサンプリング
- 各ステップで安息角（src/analyze_repose_angle.py のロジックを流用）
- 自由表面高さ（x方向ビンごとの max(z+radius)）を時系列化
- 出力:
  - metrics.csv: step/time/angles/r2/n_particles
  - repose_angle_timeseries.png: 安息角の時系列
  - surface_height_heatmap.png: 自由表面高さの(時間×x)ヒートマップ
  - surface_height_profiles.png: 代表時刻の表面プロファイル重ね描き
  - surface_height.npz: 行列データ（後処理用）

例:
  python3 scripts/analyze_deposition_transition.py \
    --file results/.../particles.csv --width 1.0 --side right \
    --nsamples 240 --sampling log \
    --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 \
    --nbins 200 --outdir results/.../transition_metrics
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _select_steps(steps: np.ndarray, nsamples: int, sampling: str) -> np.ndarray:
    steps = np.asarray(steps, dtype=int)
    if steps.size == 0:
        return steps
    if nsamples >= steps.size:
        return steps

    if sampling == "linear":
        idx = np.linspace(0, steps.size - 1, nsamples)
    elif sampling == "log":
        t = np.logspace(-3, 0, nsamples)
        t = (t - t.min()) / (t.max() - t.min())
        idx = t * (steps.size - 1)
    else:
        raise ValueError(f"sampling must be linear|log, got: {sampling}")

    idx = np.unique(np.round(idx).astype(int))
    if idx.size < nsamples:
        fill = np.linspace(0, steps.size - 1, nsamples)
        idx = np.unique(np.concatenate([idx, np.round(fill).astype(int)]))
    idx = np.clip(idx, 0, steps.size - 1)
    idx = idx[:nsamples]
    return steps[idx]


def _surface_height_profile(
    x: np.ndarray,
    z: np.ndarray,
    r: np.ndarray,
    x_edges: np.ndarray,
) -> np.ndarray:
    """xビンごとの自由表面高さ=max(z+r)。粒子が無いbinはNaN。"""
    z_top = z + r
    height = np.full(len(x_edges) - 1, np.nan, dtype=float)
    for i in range(len(x_edges) - 1):
        mask = (x >= x_edges[i]) & (x < x_edges[i + 1])
        if np.any(mask):
            height[i] = float(np.nanmax(z_top[mask]))
    return height


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="堆積遷移の定量解析（安息角・自由表面高さ）")
    parser.add_argument("--file", "-f", required=True, help="particles.csv のパス")
    parser.add_argument("--width", type=float, required=True, help="コンテナ幅 [m]")
    parser.add_argument("--outdir", type=str, required=True, help="出力ディレクトリ")

    parser.add_argument("--step-min", type=int, default=None, help="使用するステップ範囲の最小値")
    parser.add_argument("--step-max", type=int, default=None, help="使用するステップ範囲の最大値")
    parser.add_argument("--nsamples", type=int, default=240, help="サンプル数（デフォルト: 240）")
    parser.add_argument("--sampling", choices=["linear", "log"], default="log", help="サンプリング方法")

    # 安息角計算の領域指定（指定がある場合はマージン自動除外を無効化）
    parser.add_argument("--x-min", type=float, default=None, help="解析範囲 x最小 [m]")
    parser.add_argument("--x-max", type=float, default=None, help="解析範囲 x最大 [m]")
    parser.add_argument("--z-min", type=float, default=None, help="解析範囲 z最小 [m]")
    parser.add_argument("--z-max", type=float, default=None, help="解析範囲 z最大 [m]")
    parser.add_argument("--margin-center", type=float, default=0.1, help="中央除外マージン比率")
    parser.add_argument("--margin-side", type=float, default=0.1, help="左右除外マージン比率")
    parser.add_argument("--side", choices=["left", "right", "both"], default="right", help="安息角の測定対象斜面")

    # 表面高さ
    parser.add_argument("--nbins", type=int, default=200, help="x方向ビン数（デフォルト: 200）")

    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # src/analyze_repose_angle.py を import
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "src"))
    import analyze_repose_angle as aor  # noqa: E402

    df = pd.read_csv(args.file, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    required = {"step", "time", "x", "z", "radius"}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"particles.csv に必要な列がありません: required={sorted(required)}")

    # NaN行を除去（破損したCSVファイル対策）
    n_before = len(df)
    df = df.dropna(subset=["step"])
    n_after = len(df)
    if n_before != n_after:
        print(f"警告: step列にNaNを含む{n_before - n_after}行を除去しました")

    steps = np.sort(df["step"].unique().astype(int))
    if args.step_min is not None:
        steps = steps[steps >= args.step_min]
    if args.step_max is not None:
        steps = steps[steps <= args.step_max]
    if steps.size == 0:
        raise ValueError("指定範囲にステップが存在しません。")

    selected_steps = _select_steps(steps, args.nsamples, args.sampling)
    grouped = df.groupby("step", sort=False)

    # 解析範囲
    x_range: Optional[Tuple[float, float]] = None
    if args.x_min is not None or args.x_max is not None:
        x_range = (args.x_min, args.x_max)
    z_range: Optional[Tuple[float, float]] = None
    if args.z_min is not None or args.z_max is not None:
        z_range = (args.z_min, args.z_max)

    left_m, right_m, center_m = (args.margin_side, args.margin_side, args.margin_center)
    if x_range is not None:
        left_m, right_m, center_m = (0.0, 0.0, 0.0)

    # 表面高さ用のxビンは解析範囲に合わせて固定
    x_min = x_range[0] if x_range and x_range[0] is not None else 0.0
    x_max = x_range[1] if x_range and x_range[1] is not None else args.width
    x_edges = np.linspace(x_min, x_max, args.nbins + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

    metrics = []
    surface_height = np.full((selected_steps.size, args.nbins), np.nan, dtype=float)

    print(f"時系列解析: {selected_steps.size} サンプル (steps {selected_steps[0]}..{selected_steps[-1]})")
    for i, step in enumerate(selected_steps):
        snap = grouped.get_group(step)
        time_val = float(snap["time"].iloc[0])
        x = snap["x"].to_numpy(dtype=float)
        z = snap["z"].to_numpy(dtype=float)
        r = snap["radius"].to_numpy(dtype=float)

        # 範囲フィルタ（analyze_repose_angle と合わせる）
        mask = np.ones_like(x, dtype=bool)
        if x_range:
            if x_range[0] is not None:
                mask &= x >= x_range[0]
            if x_range[1] is not None:
                mask &= x <= x_range[1]
        if z_range:
            if z_range[0] is not None:
                mask &= z >= z_range[0]
            if z_range[1] is not None:
                mask &= z <= z_range[1]

        fx, fz, fr = x[mask], z[mask], r[mask]

        left_surface, right_surface = aor.detect_surface_particles(
            fx,
            fz,
            fr,
            args.width,
            left_margin=left_m,
            right_margin=right_m,
            center_margin=center_m,
            side=args.side,
        )
        left_fit = aor.calculate_repose_angle(left_surface[0], left_surface[1], side="left")
        right_fit = aor.calculate_repose_angle(right_surface[0], right_surface[1], side="right")

        left_angle, _, _, left_r2 = left_fit
        right_angle, _, _, right_r2 = right_fit

        if args.side == "left":
            avg_angle = left_angle
        elif args.side == "right":
            avg_angle = right_angle
        else:
            angles = [a for a in [left_angle, right_angle] if not np.isnan(a)]
            avg_angle = float(np.mean(angles)) if angles else np.nan

        surface_height[i, :] = _surface_height_profile(x, z, r, x_edges)

        metrics.append(
            {
                "step": int(step),
                "time": time_val,
                "n_particles": int(len(snap)),
                "left_angle_deg": float(left_angle) if np.isfinite(left_angle) else np.nan,
                "right_angle_deg": float(right_angle) if np.isfinite(right_angle) else np.nan,
                "average_angle_deg": float(avg_angle) if np.isfinite(avg_angle) else np.nan,
                "left_r2": float(left_r2) if np.isfinite(left_r2) else np.nan,
                "right_r2": float(right_r2) if np.isfinite(right_r2) else np.nan,
            }
        )

        if i == 0 or i == selected_steps.size - 1 or i % max(1, selected_steps.size // 10) == 0:
            print(f"  [{i+1:4d}/{selected_steps.size}] step={int(step)} time={time_val:.6e} avg_angle={avg_angle:.3f}")

    metrics_df = pd.DataFrame(metrics).sort_values(["time", "step"])
    metrics_csv = outdir / "metrics.csv"
    metrics_df.to_csv(metrics_csv, index=False)
    print(f"metrics を保存しました: {metrics_csv}")

    # 角度時系列プロット
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(metrics_df["time"], metrics_df["average_angle_deg"], label="average_angle_deg", color="tab:blue")
    if args.side in ("both", "left"):
        ax.plot(metrics_df["time"], metrics_df["left_angle_deg"], label="left_angle_deg", alpha=0.5)
    if args.side in ("both", "right"):
        ax.plot(metrics_df["time"], metrics_df["right_angle_deg"], label="right_angle_deg", alpha=0.5)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("repose_angle [deg]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    angle_png = outdir / "repose_angle_timeseries.png"
    fig.savefig(angle_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 表面高さヒートマップ
    times = metrics_df["time"].to_numpy(dtype=float)
    # 表示用にtime順へ並べ替え（metrics_dfはtime/stepでsort済み）
    # ただし surface_height は selected_steps の順なので、同じ順序を作る
    order = np.argsort(times)
    times_sorted = times[order]
    surface_sorted = surface_height[order, :]

    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(
        surface_sorted,
        aspect="auto",
        origin="lower",
        extent=(x_centers[0], x_centers[-1], times_sorted[0], times_sorted[-1]),
        cmap="viridis",
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("time [s]")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("surface_height z_top [m]")
    fig.tight_layout()
    heat_png = outdir / "surface_height_heatmap.png"
    fig.savefig(heat_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 代表時刻の表面プロファイル（5本）
    k = min(5, times_sorted.size)
    pick = np.linspace(0, times_sorted.size - 1, k).round().astype(int)
    fig, ax = plt.subplots(figsize=(10, 4))
    for j in pick:
        ax.plot(x_centers, surface_sorted[j, :], label=f"t={times_sorted[j]:.3e}s")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z_top [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    prof_png = outdir / "surface_height_profiles.png"
    fig.savefig(prof_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # 行列データ保存（後処理用）
    npz_path = outdir / "surface_height.npz"
    np.savez_compressed(
        npz_path,
        steps=selected_steps.astype(int),
        times=times.astype(float),
        x_centers=x_centers.astype(float),
        surface_height=surface_height.astype(float),
    )
    print(f"表面高さ行列を保存しました: {npz_path}")
    print(f"プロットを保存しました: {angle_png}, {heat_png}, {prof_png}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())





















