#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 168:00:00
#SBATCH -o log/charge_dist_aor_multirun_stdout.%A_%a.log
#SBATCH -e log/charge_dist_aor_multirun_stderr.%A_%a.log
#SBATCH --array=1-60

# ==========================================
# Phase 2 追加シード: 帯電分布×カットオフ距離の安息角計測
#
# - シード2-5 の安息角計測を実行（シード1は既存結果を使用）
# - 4シード × 3分布 × 5カットオフ = 60ケース
# - カットオフ: 0.02, 0.015, 0.01, 0.006, 0.003 m
#
# 使い方:
#   # Phase 1 完了後に投入
#   sbatch --dependency=afterok:<PHASE1_JOB_ID> job_charge_dist_aor_multirun.sh
# ==========================================

set -euo pipefail

module load intel
ulimit -s unlimited
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-32}

date

BASE_INPUT="inputs/charge_dist_study/input_aor_template.dat"
FILLING_BASE="results/charge_distribution_study/phase1_filling"
RESULTS_BASE="results/charge_distribution_study/phase2_aor"
REPOSE_SIDE="right"

if [ ! -f "$BASE_INPUT" ]; then
  echo "Error: Base input file '$BASE_INPUT' not found."
  exit 1
fi

# 乱数シードの配列（シード2-5、シード1は既存）
SEED_IDS=("seed2" "seed3" "seed4" "seed5")
NUM_SEEDS=${#SEED_IDS[@]}

# 分布の配列
DISTRIBUTIONS=("bimodal" "normal" "uniform")
NUM_DISTS=${#DISTRIBUTIONS[@]}

# カットオフ距離の配列
CUTOFFS=("0.02" "0.015" "0.01" "0.006" "0.003")
NUM_CUTOFFS=${#CUTOFFS[@]}

# Array ID から シード、分布、カットオフを決定
# Array ID 1-5:   seed2 × bimodal × (0.02, 0.015, 0.01, 0.006, 0.003)
# Array ID 6-10:  seed2 × normal × (0.02, 0.015, 0.01, 0.006, 0.003)
# Array ID 11-15: seed2 × uniform × (0.02, 0.015, 0.01, 0.006, 0.003)
# Array ID 16-20: seed3 × bimodal × ...
# ...

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

# 各シードには NUM_DISTS × NUM_CUTOFFS = 15 ケース
CASES_PER_SEED=$((NUM_DISTS * NUM_CUTOFFS))

SEED_IDX=$(( (TASK_ID - 1) / CASES_PER_SEED ))
REMAINING=$(( (TASK_ID - 1) % CASES_PER_SEED ))
DIST_IDX=$(( REMAINING / NUM_CUTOFFS ))
CUTOFF_IDX=$(( REMAINING % NUM_CUTOFFS ))

if [ $SEED_IDX -ge $NUM_SEEDS ] || [ $DIST_IDX -ge $NUM_DISTS ] || [ $CUTOFF_IDX -ge $NUM_CUTOFFS ]; then
  echo "Error: Invalid array task ID: $TASK_ID"
  echo "  SEED_IDX=$SEED_IDX, DIST_IDX=$DIST_IDX, CUTOFF_IDX=$CUTOFF_IDX"
  exit 1
fi

SEED_ID="${SEED_IDS[$SEED_IDX]}"
DIST_TYPE="${DISTRIBUTIONS[$DIST_IDX]}"
RC="${CUTOFFS[$CUTOFF_IDX]}"

FILLING_CASE="${DIST_TYPE}_${SEED_ID}"
CASE_NAME="${DIST_TYPE}_${SEED_ID}_rc_${RC}"

# Phase 1 で生成された充填ファイルを確認
FILLED_FILE="${FILLING_BASE}/${FILLING_CASE}/filled_particles.dat"
if [ ! -f "$FILLED_FILE" ]; then
  echo "Error: Filled particles file not found: $FILLED_FILE"
  echo "Please run Phase 1 (job_charge_dist_filling_multirun.sh) first."
  exit 1
fi

RUN_ID="job_${SLURM_ARRAY_JOB_ID:-local}_task_${SLURM_ARRAY_TASK_ID:-0}"
OUTPUT_DIR="${RESULTS_BASE}/${CASE_NAME}"

TEMP_INPUT_DIR="inputs/charge_dist_study/aor_multirun/${CASE_NAME}"
TEMP_INPUT="${TEMP_INPUT_DIR}/input_${RUN_ID}.dat"

mkdir -p "$OUTPUT_DIR" "$TEMP_INPUT_DIR" log

# 充填ファイルを inputs/ にコピー（シミュレーションが読み込めるように）
cp "$FILLED_FILE" "inputs/filled_particles.dat"

cp "$BASE_INPUT" "$TEMP_INPUT"

# COULOMB_CUTOFF を上書き
sed -i "s/^COULOMB_CUTOFF[[:space:]].*/COULOMB_CUTOFF              ${RC}/" "$TEMP_INPUT"

# CONTAINER_WIDTH を input から取得
CONTAINER_WIDTH=$(awk '/^CONTAINER_WIDTH/{print $2; exit}' "$TEMP_INPUT")
if [ -z "${CONTAINER_WIDTH}" ]; then
  echo "Warning: failed to read CONTAINER_WIDTH from input; fallback to 2.0"
  CONTAINER_WIDTH="2.0"
fi

echo "=========================================="
echo "Phase 2 Multirun: AoR Measurement"
echo "=========================================="
echo "Seed ID         : $SEED_ID"
echo "Distribution    : $DIST_TYPE"
echo "Cutoff [m]      : $RC"
echo "Case Name       : $CASE_NAME"
echo "Filled File     : $FILLED_FILE"
echo "Input File      : $TEMP_INPUT"
echo "Output Dir      : $OUTPUT_DIR"
echo "Repose Side     : $REPOSE_SIDE"
echo "Container Width : $CONTAINER_WIDTH"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "=========================================="

# ==========================================
# Compile
# ==========================================
echo ""
echo "[Phase 2.1] Compiling..."

BUILD_DIR="build/charge_dist_aor_multirun/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module "$BUILD_DIR" -o "$EXE_NAME" "$SRC_FILE"

echo "Compilation successful."

# ==========================================
# Run AoR Simulation
# ==========================================
echo ""
echo "[Phase 2.2] Running AoR simulation..."

START_TIME=$(date +%s)
srun "$EXE_NAME" "$TEMP_INPUT" "$OUTPUT_DIR"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "Simulation completed successfully."
echo "Elapsed time: ${ELAPSED} seconds"

# タイミング情報を保存
echo "case,distribution,seed_id,cutoff_m,elapsed_seconds" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${DIST_TYPE},${SEED_ID},${RC},${ELAPSED}" >> "${OUTPUT_DIR}/timing.csv"

# ビルドディレクトリをクリーンアップ
rm -rf "$BUILD_DIR"

# 一時入力ファイルを出力ディレクトリに保存
cp "$TEMP_INPUT" "${OUTPUT_DIR}/input_${CASE_NAME}.dat"
rm -f "$TEMP_INPUT"

# ==========================================
# AoR Analysis
# ==========================================
echo ""
echo "[Phase 2.3] Measuring angle of repose..."

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
echo "Phase 2 Multirun Task completed: $CASE_NAME"
echo "Seed ID       : $SEED_ID"
echo "Distribution  : $DIST_TYPE"
echo "Cutoff [m]    : $RC"
echo "Elapsed time  : ${ELAPSED} seconds"
echo "Output Dir    : $OUTPUT_DIR"
echo "=========================================="
date







