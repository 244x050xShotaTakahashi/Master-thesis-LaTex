#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rc掃引（AoR）結果を集計するスクリプト。

期待するディレクトリ構造（例）:
  results/aor_cutoff_sweep/rc_0.03/job_XXXX_task_Y/
    - repose_angle_results.csv
    - timing_report.csv         (dem_valid.f90 の profiler_write_csv)
    - timing.csv                (ジョブスクリプトの壁時計)

出力:
  - sweep_summary.csv（各ケースのAoR/時間内訳/ΔAoR）
  - 推奨rc（|ΔAoR|<=tol を満たす最小rc）を標準出力
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


RC_RE = re.compile(r"rc_(\d+(?:\.\d+)?)")


def read_kv_csv(path: Path) -> Dict[str, str]:
    """
    repose_angle_results.csv のような 'parameter,value' 形式を dict で返す。
    """
    out: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return out
        for row in reader:
            if len(row) < 2:
                continue
            k = row[0].strip()
            v = row[1].strip()
            if k:
                out[k] = v
    return out


def read_timing_report(path: Path) -> Dict[str, float]:
    """
    timing_report.csv を {name: total_seconds} で返す。
    """
    out: Dict[str, float] = {}
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("name") or "").strip()
            total_seconds = row.get("total_seconds")
            if not name or total_seconds is None:
                continue
            try:
                out[name] = float(total_seconds)
            except ValueError:
                continue
    return out


def read_wallclock_seconds(path: Path) -> Optional[float]:
    """
    timing.csv（case,rc_m,elapsed_seconds）から elapsed_seconds を読む。
    """
    try:
        with path.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                v = row.get("elapsed_seconds")
                if v is None:
                    continue
                return float(v)
    except Exception:
        return None
    return None


def parse_rc_from_path(case_dir: Path) -> Optional[float]:
    m = RC_RE.search(str(case_dir))
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


class SweepRow:
    def __init__(
        self,
        case_dir: Path,
        rc_m: Optional[float],
        aor_deg: Optional[float],
        sim_time_s: Optional[float],
        wallclock_s: Optional[float],
        coulomb_s: Optional[float],
        coulomb_pct: Optional[float],
        neighbor_s: Optional[float],
        integrate_s: Optional[float],
        output_s: Optional[float],
    ):
        self.case_dir = case_dir
        self.rc_m = rc_m
        self.aor_deg = aor_deg
        self.sim_time_s = sim_time_s
        self.wallclock_s = wallclock_s
        self.coulomb_s = coulomb_s
        self.coulomb_pct = coulomb_pct
        self.neighbor_s = neighbor_s
        self.integrate_s = integrate_s
        self.output_s = output_s


def safe_float(s: Optional[str]) -> Optional[float]:
    if s is None:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def collect_cases(root: Path) -> List[SweepRow]:
    rows: List[SweepRow] = []
    for aor_file in root.glob("**/repose_angle_results.csv"):
        case_dir = aor_file.parent
        aor_kv = read_kv_csv(aor_file)
        aor_deg = safe_float(aor_kv.get("average_angle_deg"))
        sim_time_s = safe_float(aor_kv.get("simulation_time"))
        rc_m = parse_rc_from_path(case_dir)

        timing_report = case_dir / "timing_report.csv"
        tmap = read_timing_report(timing_report) if timing_report.exists() else {}
        total_all = sum(tmap.values()) if tmap else None

        coulomb_s = tmap.get("coulomb_force")
        neighbor_s = tmap.get("neighbor_search_contact")
        integrate_s = tmap.get("integrate_leapfrog")
        output_s = tmap.get("output")
        coulomb_pct = (coulomb_s / total_all * 100.0) if (coulomb_s is not None and total_all and total_all > 0) else None

        wallclock_s = read_wallclock_seconds(case_dir / "timing.csv")

        rows.append(
            SweepRow(
                case_dir=case_dir,
                rc_m=rc_m,
                aor_deg=aor_deg,
                sim_time_s=sim_time_s,
                wallclock_s=wallclock_s,
                coulomb_s=coulomb_s,
                coulomb_pct=coulomb_pct,
                neighbor_s=neighbor_s,
                integrate_s=integrate_s,
                output_s=output_s,
            )
        )
    rows.sort(key=lambda r: (r.rc_m is None, r.rc_m if r.rc_m is not None else 1e9, str(r.case_dir)))
    return rows


def write_summary_csv(rows: List[SweepRow], out_csv: Path, ref_rc: float) -> None:
    # 参照AoR（ref_rc）を探す
    ref_aor = None
    for r in rows:
        if r.rc_m is not None and abs(r.rc_m - ref_rc) < 1e-12:
            ref_aor = r.aor_deg
            break

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "case_dir",
                "rc_m",
                "aor_deg",
                "delta_aor_deg_vs_ref",
                "sim_time_s",
                "wallclock_s",
                "coulomb_s",
                "coulomb_pct",
                "neighbor_s",
                "integrate_s",
                "output_s",
            ]
        )
        for r in rows:
            delta = None
            if ref_aor is not None and r.aor_deg is not None:
                delta = r.aor_deg - ref_aor
            writer.writerow(
                [
                    str(r.case_dir),
                    r.rc_m,
                    r.aor_deg,
                    delta,
                    r.sim_time_s,
                    r.wallclock_s,
                    r.coulomb_s,
                    r.coulomb_pct,
                    r.neighbor_s,
                    r.integrate_s,
                    r.output_s,
                ]
            )


def choose_rc(rows: List[SweepRow], ref_rc: float, tol_deg: float) -> Tuple[Optional[float], Optional[float]]:
    ref_aor = None
    for r in rows:
        if r.rc_m is not None and abs(r.rc_m - ref_rc) < 1e-12:
            ref_aor = r.aor_deg
            break
    if ref_aor is None:
        return None, None

    ok: List[float] = []
    for r in rows:
        if r.rc_m is None or r.aor_deg is None:
            continue
        if abs(r.aor_deg - ref_aor) <= tol_deg:
            ok.append(r.rc_m)
    if not ok:
        return None, ref_aor
    return min(ok), ref_aor


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=str, default="results/aor_cutoff_sweep", help="掃引結果のルートディレクトリ")
    p.add_argument("--out", type=str, default="results/aor_cutoff_sweep/sweep_summary.csv", help="集計CSV出力先")
    p.add_argument("--ref-rc", type=float, default=0.05, help="ΔAoRの参照にするrc")
    p.add_argument("--tol-deg", type=float, default=0.5, help="許容AoR差 [deg]")
    args = p.parse_args()

    root = Path(args.root)
    rows = collect_cases(root)
    if not rows:
        raise SystemExit(f"No repose_angle_results.csv found under: {root}")

    out_csv = Path(args.out)
    write_summary_csv(rows, out_csv, ref_rc=args.ref_rc)
    print(f"Saved: {out_csv}")

    rec_rc, ref_aor = choose_rc(rows, ref_rc=args.ref_rc, tol_deg=args.tol_deg)
    if rec_rc is None:
        print(f"[WARN] 推奨rcを決められませんでした（ref_rc={args.ref_rc}, tol={args.tol_deg}°）。")
        return

    print("----")
    print(f"ref_rc = {args.ref_rc} m, ref_aor = {ref_aor:.4f} deg")
    print(f"tolerance = ±{args.tol_deg} deg")
    print(f"recommended_rc = {rec_rc} m")


if __name__ == "__main__":
    main()















