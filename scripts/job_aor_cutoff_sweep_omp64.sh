#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=1:t=64:c=64
#SBATCH -t 168:00:00
#SBATCH -o log/stdout.%A_%a.omp64.log
#SBATCH -e log/stderr.%A_%a.omp64.log
#SBATCH --array=1-4

# ==========================================
# AoR用: クーロンカットオフ rc 掃引（短縮ラン / AUTO_WALL_WITHDRAW）
#   - 64スレッド（OMP_NUM_THREADS=64）版
#
# 32スレッド版: job_aor_cutoff_sweep.sh
#
# 出力先を分けて比較しやすくする:
#   results/aor_cutoff_sweep_omp64/...
#
# 使い方:
#   sbatch job_aor_cutoff_sweep_omp64.sh
# ==========================================

set -euo pipefail

module load intel
ulimit -s unlimited

export OMP_NUM_THREADS=64
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

date

BASE_INPUT="inputs/input_coulomb_short.dat"
REPOSE_SIDE="right"

if [ ! -f "$BASE_INPUT" ]; then
  echo "Error: Base input file '$BASE_INPUT' not found."
  exit 1
fi

case ${SLURM_ARRAY_TASK_ID:-1} in
  1) RC="0.02"; CASE_NAME="rc_0.02" ;;
  2) RC="0.03"; CASE_NAME="rc_0.03" ;;
  3) RC="0.04"; CASE_NAME="rc_0.04" ;;
  4) RC="0.05"; CASE_NAME="rc_0.05" ;;
  *) echo "Error: Invalid array task ID: ${SLURM_ARRAY_TASK_ID:-unset}"; exit 1 ;;
esac

RUN_ID="job_${SLURM_ARRAY_JOB_ID:-local}_task_${SLURM_ARRAY_TASK_ID:-0}"
OUTPUT_DIR="results/aor_cutoff_sweep_omp64/${CASE_NAME}/${RUN_ID}"

TEMP_INPUT_DIR="inputs/aor_cutoff_sweep_omp64/${CASE_NAME}"
TEMP_INPUT="${TEMP_INPUT_DIR}/input_${RUN_ID}.dat"

mkdir -p "$OUTPUT_DIR" "$TEMP_INPUT_DIR" log

cp "$BASE_INPUT" "$TEMP_INPUT"

# COULOMB_CUTOFF を上書き
sed -i "s/^COULOMB_CUTOFF.*/COULOMB_CUTOFF              ${RC}/" "$TEMP_INPUT"

# 念のため、プロファイル出力のためサンプリングを無効化
if grep -q "^PROFILING_SAMPLE_INTERVAL" "$TEMP_INPUT"; then
  sed -i "s/^PROFILING_SAMPLE_INTERVAL.*/PROFILING_SAMPLE_INTERVAL   0/" "$TEMP_INPUT"
else
  echo "PROFILING_SAMPLE_INTERVAL   0" >> "$TEMP_INPUT"
fi

# CONTAINER_WIDTH を input から取得し、AoR解析の --width に渡す（不一致を防ぐ）
CONTAINER_WIDTH=$(awk '/^CONTAINER_WIDTH/{print $2; exit}' "$TEMP_INPUT")
if [ -z "${CONTAINER_WIDTH}" ]; then
  echo "Warning: failed to read CONTAINER_WIDTH from input; fallback to 2.0"
  CONTAINER_WIDTH="2.0"
fi

echo "=========================================="
echo "AoR cutoff sweep (AUTO_WALL_WITHDRAW) - OMP64"
echo "=========================================="
echo "Case Name       : $CASE_NAME"
echo "rc [m]          : $RC"
echo "Input File      : $TEMP_INPUT"
echo "Output Dir      : $OUTPUT_DIR"
echo "Repose Side     : $REPOSE_SIDE"
echo "Container Width : $CONTAINER_WIDTH"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "OMP_PLACES      : ${OMP_PLACES}"
echo "OMP_PROC_BIND   : ${OMP_PROC_BIND}"
echo "=========================================="

# ==========================================
# Phase 1: Compile
# ==========================================
echo ""
echo "[Phase 1] Compiling..."

BUILD_DIR="build/job_aor_cutoff_sweep_omp64/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module "$BUILD_DIR" -o "$EXE_NAME" "$SRC_FILE"

echo "Compilation successful."

# ==========================================
# Phase 2: Run
# ==========================================
echo ""
echo "[Phase 2] Running DEM simulation..."

START_TIME=$(date +%s)
srun "$EXE_NAME" "$TEMP_INPUT" "$OUTPUT_DIR"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "Simulation completed successfully."
echo "Elapsed time: ${ELAPSED} seconds"

echo "case,rc_m,elapsed_seconds,omp_threads" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${RC},${ELAPSED},${OMP_NUM_THREADS}" >> "${OUTPUT_DIR}/timing.csv"

rm -rf "$BUILD_DIR"
cp "$TEMP_INPUT" "${OUTPUT_DIR}/input_${CASE_NAME}.dat"
rm -f "$TEMP_INPUT"

# ==========================================
# Phase 3: AoR analysis
# ==========================================
echo ""
echo "[Phase 3] Measuring angle of repose..."

DATA_FILE="${OUTPUT_DIR}/particles.csv"
ANALYSIS_SCRIPT="src/analyze_repose_angle.py"

if [ -f "$DATA_FILE" ] && [ -f "$ANALYSIS_SCRIPT" ]; then
  python3 "$ANALYSIS_SCRIPT" \
    --data "$DATA_FILE" \
    --output "$OUTPUT_DIR" \
    --name "repose_angle" \
    --side "$REPOSE_SIDE" \
    --width "$CONTAINER_WIDTH"
  echo "Repose angle analysis completed."
else
  echo "Warning: particles.csv or analyze_repose_angle.py not found. Skipping."
fi

echo ""
echo "=========================================="
echo "Task completed: $CASE_NAME (OMP64)"
echo "rc [m]        : $RC"
echo "OMP Threads   : $OMP_NUM_THREADS"
echo "Elapsed time  : ${ELAPSED} seconds"
echo "Output Dir    : $OUTPUT_DIR"
echo "=========================================="
date


































