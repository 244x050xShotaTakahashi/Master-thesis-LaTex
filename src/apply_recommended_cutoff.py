#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
掃引集計（sweep_summary.csv）から推奨rcを選び、指定inputの COULOMB_CUTOFF を更新する。

使い方（例）:
  python3 src/summarize_cutoff_sweep.py --root results/aor_cutoff_sweep --out results/aor_cutoff_sweep/sweep_summary.csv
  python3 src/apply_recommended_cutoff.py --summary results/aor_cutoff_sweep/sweep_summary.csv --input inputs/input_coulomb.dat --tol-deg 0.5 --ref-rc 0.05
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def read_summary(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_float(s: Optional[str]) -> Optional[float]:
    if s is None:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def choose_rc_from_summary(rows: List[Dict[str, str]], ref_rc: float, tol_deg: float) -> Tuple[Optional[float], Optional[float]]:
    # ref AoR
    ref_aor = None
    for r in rows:
        rc = parse_float(r.get("rc_m"))
        if rc is not None and abs(rc - ref_rc) < 1e-12:
            ref_aor = parse_float(r.get("aor_deg"))
            break
    if ref_aor is None:
        return None, None

    ok: List[float] = []
    for r in rows:
        rc = parse_float(r.get("rc_m"))
        aor = parse_float(r.get("aor_deg"))
        if rc is None or aor is None:
            continue
        if abs(aor - ref_aor) <= tol_deg:
            ok.append(rc)
    if not ok:
        return None, ref_aor
    return min(ok), ref_aor


def replace_coulomb_cutoff(text: str, rc: float) -> str:
    # 既存行があれば置換、なければ末尾に追記
    pat = re.compile(r"^COULOMB_CUTOFF\s+.*$", re.MULTILINE)
    new_line = f"COULOMB_CUTOFF              {rc:g}"
    if pat.search(text):
        return pat.sub(new_line, text, count=1)
    # 無ければ末尾
    if not text.endswith("\n"):
        text += "\n"
    return text + new_line + "\n"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="sweep_summary.csv のパス")
    p.add_argument("--input", required=True, help="更新する input .dat のパス")
    p.add_argument("--ref-rc", type=float, default=0.05)
    p.add_argument("--tol-deg", type=float, default=0.5)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    summary_path = Path(args.summary)
    input_path = Path(args.input)

    rows = read_summary(summary_path)
    rec_rc, ref_aor = choose_rc_from_summary(rows, ref_rc=args.ref_rc, tol_deg=args.tol_deg)
    if rec_rc is None:
        raise SystemExit(f"推奨rcを決められません（ref_rc={args.ref_rc}, tol={args.tol_deg}°）。")

    text = input_path.read_text(encoding="utf-8")
    new_text = replace_coulomb_cutoff(text, rec_rc)

    print(f"ref_rc={args.ref_rc} m, ref_aor={ref_aor:.4f} deg")
    print(f"tolerance=±{args.tol_deg} deg")
    print(f"recommended_rc={rec_rc:g} m")
    print(f"target_input={input_path}")
    if args.dry_run:
        print("[DRY RUN] inputは更新しません。")
        return

    input_path.write_text(new_text, encoding="utf-8")
    print("COULOMB_CUTOFF を更新しました。")


if __name__ == "__main__":
    main()















