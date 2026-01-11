#!/usr/bin/env python3
"""
堆積遷移の可視化（フレーム一括生成 + 任意で動画化）

particles.csv からステップをサンプリングし、plot_snapshot.py 相当の描画で
フレーム画像を出力します。ffmpeg が利用可能なら mp4 への結合も行えます。

例:
  python3 scripts/make_transition_movie.py \
    --file results/.../particles.csv --walls inputs/walls.dat --width 1.0 \
    --color-by charge --dpi 300 --nframes 240 --sampling log \
    --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 \
    --frames-dir results/.../frames_charge_zoom \
    --mp4 results/.../transition_charge_zoom.mp4 --fps 30
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd


def _select_steps(steps: np.ndarray, nframes: int, sampling: str) -> np.ndarray:
    steps = np.asarray(steps, dtype=int)
    if steps.size == 0:
        return steps
    if nframes >= steps.size:
        return steps

    if sampling == "linear":
        idx = np.linspace(0, steps.size - 1, nframes)
    elif sampling == "log":
        # 早い時刻を密、終盤を疎に（変化が大きい前半を見やすく）
        # logspace は端点を含むように 0..1 を指数分布させる
        t = np.logspace(-3, 0, nframes)  # 1e-3..1
        t = (t - t.min()) / (t.max() - t.min())  # 0..1
        idx = t * (steps.size - 1)
    else:
        raise ValueError(f"sampling must be linear|log, got: {sampling}")

    idx = np.unique(np.round(idx).astype(int))
    # uniqueで減った分は線形で埋める
    if idx.size < nframes:
        fill = np.linspace(0, steps.size - 1, nframes)
        idx = np.unique(np.concatenate([idx, np.round(fill).astype(int)]))
    idx = np.clip(idx, 0, steps.size - 1)
    idx = idx[:nframes]
    return steps[idx]


def _run_ffmpeg_mp4(frames_dir: Path, fps: int, mp4_path: Path) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("ffmpeg が見つかりません（PATHに無い）。フレーム画像は生成済みなので、別環境で結合してください。")

    pattern = str(frames_dir / "frame_%06d.png")
    cmd = [
        ffmpeg,
        "-y",
        "-framerate",
        str(fps),
        "-i",
        pattern,
        "-pix_fmt",
        "yuv420p",
        str(mp4_path),
    ]
    subprocess.run(cmd, check=True)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="堆積遷移のフレーム生成と動画化（任意）")
    parser.add_argument("--file", "-f", type=str, required=True, help="particles.csv のパス")
    parser.add_argument("--walls", "-w", type=str, default="inputs/walls.dat", help="walls.dat のパス")
    parser.add_argument("--width", type=float, required=True, help="コンテナ幅 [m]（plot_snapshotと同じ）")
    parser.add_argument("--height", type=float, default=0.0, help="コンテナ高さ [m]（0なら上壁なし）")

    parser.add_argument("--step-min", type=int, default=None, help="使用するステップ範囲の最小値")
    parser.add_argument("--step-max", type=int, default=None, help="使用するステップ範囲の最大値")
    parser.add_argument("--nframes", type=int, default=240, help="生成するフレーム数（デフォルト: 240）")
    parser.add_argument("--sampling", choices=["linear", "log"], default="log", help="ステップのサンプリング方法")

    parser.add_argument("--frames-dir", type=str, required=True, help="フレーム出力先ディレクトリ")
    parser.add_argument("--color-by", choices=["radius", "charge"], default="radius", help="粒子の着色に使う量")
    parser.add_argument("--charge-mode", choices=["actual", "zero"], default="actual",
                        help="color-by=charge時の可視化モード (actual=実データ, zero=0C固定)")
    parser.add_argument("--dpi", type=int, default=200, help="保存dpi（デフォルト: 200）")
    parser.add_argument("--x-min", type=float, default=None, help="表示範囲 x最小 [m]")
    parser.add_argument("--x-max", type=float, default=None, help="表示範囲 x最大 [m]")
    parser.add_argument("--z-min", type=float, default=None, help="表示範囲 z最小 [m]")
    parser.add_argument("--z-max", type=float, default=None, help="表示範囲 z最大 [m]")

    parser.add_argument("--mp4", type=str, default=None, help="mp4 出力パス（指定時はffmpegで結合）")
    parser.add_argument("--fps", type=int, default=30, help="動画fps（デフォルト: 30）")

    args = parser.parse_args(argv)

    particles_path = Path(args.file)
    frames_dir = Path(args.frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    # scripts/plot_snapshot.py を import するために scripts をsys.pathへ
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "scripts"))
    import plot_snapshot as ps  # noqa: E402

    df = ps.load_particles(str(particles_path))
    if "step" not in df.columns or "time" not in df.columns:
        raise ValueError("particles.csv に step/time 列がありません。")

    steps = np.sort(df["step"].unique().astype(int))
    if args.step_min is not None:
        steps = steps[steps >= args.step_min]
    if args.step_max is not None:
        steps = steps[steps <= args.step_max]
    if steps.size == 0:
        raise ValueError("指定範囲にステップが存在しません。--step-min/--step-max を見直してください。")

    selected_steps = _select_steps(steps, args.nframes, args.sampling)
    walls = ps.load_walls(args.walls)
    if walls:
        print(f"斜面壁を読み込みました: {len(walls)} 本")

    # step -> dataframe の参照を高速化
    grouped = df.groupby("step", sort=False)

    print(f"フレーム生成: {selected_steps.size} 枚 (steps {selected_steps[0]}..{selected_steps[-1]})")
    for i, step in enumerate(selected_steps, start=1):
        snapshot = grouped.get_group(step)
        out = frames_dir / f"frame_{i:06d}.png"
        ps.plot_snapshot(
            snapshot=snapshot,
            walls=walls,
            container_width=args.width,
            container_height=args.height,
            output_path=str(out),
            title=None,
            color_by=args.color_by,
            charge_mode=args.charge_mode,
            dpi=args.dpi,
            x_min=args.x_min,
            x_max=args.x_max,
            z_min=args.z_min,
            z_max=args.z_max,
        )
        if i == 1 or i == selected_steps.size or i % max(1, selected_steps.size // 10) == 0:
            t = float(snapshot["time"].iloc[0])
            print(f"  [{i:4d}/{selected_steps.size}] step={int(step)} time={t:.6e} -> {out.name}")

    if args.mp4:
        mp4_path = Path(args.mp4)
        mp4_path.parent.mkdir(parents=True, exist_ok=True)
        print("ffmpeg による mp4 結合を実行します...")
        _run_ffmpeg_mp4(frames_dir, args.fps, mp4_path)
        print(f"mp4 を保存しました: {mp4_path}")
    else:
        print("mp4 は未生成です（--mp4 を指定すると ffmpeg で結合します）。")
        print(f"例: ffmpeg -y -framerate {args.fps} -i {frames_dir}/frame_%06d.png -pix_fmt yuv420p out.mp4")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


