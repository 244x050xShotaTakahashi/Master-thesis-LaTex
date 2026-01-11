#!/bin/bash
set -euo pipefail

# Phase 2 AoR の Race Condition 回避（EXPLICIT_POSITIONS_FILE 直参照）を確認する簡易チェッカー
#
# 使い方:
#   bash scripts/verify_phase2_aor_no_race.sh log/charge_dist_aor_stdout.<JOBID>_*.log
#   bash scripts/verify_phase2_aor_no_race.sh log/charge_dist_aor_stdout.6400586_*.log
#
# 期待:
#   - "粒子配置: 明示座標ファイルを使用(指定): <...>/phase1_filling/<dist>/filled_particles.dat" が出る
#   - "粒子配置: 充填状態ファイルを使用: inputs/filled_particles.dat" が出ない

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <stdout_log1> [stdout_log2 ...]" >&2
  exit 2
fi

bad=0

for f in "$@"; do
  if [ ! -f "$f" ]; then
    echo "[SKIP] not found: $f"
    continue
  fi

  echo "---- $f ----"

  if grep -q "粒子配置: 充填状態ファイルを使用: inputs/filled_particles.dat" "$f"; then
    echo "[NG] uses shared inputs/filled_particles.dat"
    bad=1
  else
    echo "[OK] does not use shared inputs/filled_particles.dat"
  fi

  if grep -q "粒子配置: 明示座標ファイルを使用(指定):" "$f"; then
    echo "[OK] uses EXPLICIT_POSITIONS_FILE (direct reference)"
    grep "粒子配置: 明示座標ファイルを使用(指定):" "$f" | head -n 1
  else
    echo "[NG] missing EXPLICIT_POSITIONS_FILE message"
    bad=1
  fi
done

if [ "$bad" -ne 0 ]; then
  echo ""
  echo "Result: NG"
  exit 1
fi

echo ""
echo "Result: OK"



