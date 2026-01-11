#!/bin/bash
#SBATCH -p gr20001b
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 168:00:00
#SBATCH -o log/charge_dist_aor_stdout.%A_%a.log
#SBATCH -e log/charge_dist_aor_stderr.%A_%a.log
#SBATCH --array=1-54

# ==========================================
# Phase 2: 帯電分布×カットオフ距離の安息角計測（複数シード対応）
#
# - Phase 1 で生成された充填ファイルを入力として使用
# - 3分布 × 6カットオフ × 3シード = 54ケースを並列実行
# - カットオフ: 0.05, 0.02, 0.015, 0.01, 0.006, 0.003 m
#
# Array ID割り当て (54ケース):
#   1-18:  bimodal × 6カットオフ × 3シード
#   19-36: normal × 6カットオフ × 3シード
#   37-54: uniform × 6カットオフ × 3シード
#
# 詳細:
#   1-6:   bimodal_seed1 × (0.05, 0.02, 0.015, 0.01, 0.006, 0.003)
#   7-12:  bimodal_seed2 × (0.05, 0.02, 0.015, 0.01, 0.006, 0.003)
#   13-18: bimodal_seed3 × (0.05, 0.02, 0.015, 0.01, 0.006, 0.003)
#   19-24: normal_seed1 × ...
#   ...
#
# 使い方:
#   # Phase 1 完了後に投入
#   sbatch --dependency=afterok:<PHASE1_JOB_ID> job_charge_dist_aor.sh
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

# カットオフ距離の配列
CUTOFFS=("0.05" "0.02" "0.015" "0.01" "0.006" "0.003")
NUM_CUTOFFS=${#CUTOFFS[@]}

# 帯電分布の配列
DISTRIBUTIONS=("bimodal" "normal" "uniform")
NUM_DISTS=${#DISTRIBUTIONS[@]}

# シード数
NUM_SEEDS=3

# Array ID から分布、シード、カットオフを決定
# 構造: [分布][シード][カットオフ]
# 1ケース = 1分布 × 1シード × 1カットオフ
# 1シードあたり = 6カットオフ
# 1分布あたり = 3シード × 6カットオフ = 18ケース
# 合計 = 3分布 × 18 = 54ケース

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
CASES_PER_SEED=$NUM_CUTOFFS
CASES_PER_DIST=$((NUM_SEEDS * CASES_PER_SEED))

DIST_IDX=$(( (TASK_ID - 1) / CASES_PER_DIST ))
REMAINING=$(( (TASK_ID - 1) % CASES_PER_DIST ))
SEED_IDX=$(( REMAINING / CASES_PER_SEED ))
CUTOFF_IDX=$(( REMAINING % CASES_PER_SEED ))

if [ $DIST_IDX -ge $NUM_DISTS ] || [ $SEED_IDX -ge $NUM_SEEDS ] || [ $CUTOFF_IDX -ge $NUM_CUTOFFS ]; then
  echo "Error: Invalid array task ID: $TASK_ID"
  echo "  DIST_IDX=$DIST_IDX, SEED_IDX=$SEED_IDX, CUTOFF_IDX=$CUTOFF_IDX"
  exit 1
fi

DIST_TYPE="${DISTRIBUTIONS[$DIST_IDX]}"
SEED_NUM=$((SEED_IDX + 1))
RC="${CUTOFFS[$CUTOFF_IDX]}"
CASE_NAME="${DIST_TYPE}_seed${SEED_NUM}_rc_${RC}"

# Phase 1 で生成された充填ファイルを確認（シード別）
FILLED_FILE="${FILLING_BASE}/${DIST_TYPE}_seed${SEED_NUM}/filled_particles.dat"
if [ ! -f "$FILLED_FILE" ]; then
  echo "Error: Filled particles file not found: $FILLED_FILE"
  echo "Please run Phase 1 (job_charge_dist_filling.sh) first."
  exit 1
fi

RUN_ID="job_${SLURM_ARRAY_JOB_ID:-local}_task_${SLURM_ARRAY_TASK_ID:-0}"
OUTPUT_DIR="${RESULTS_BASE}/${CASE_NAME}"

TEMP_INPUT_DIR="inputs/charge_dist_study/aor/${CASE_NAME}"
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
echo "Phase 2: AoR Measurement"
echo "=========================================="
echo "Distribution    : $DIST_TYPE"
echo "Seed Number     : $SEED_NUM"
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
# Phase 2.1: Compile
# ==========================================
echo ""
echo "[Phase 2.1] Compiling..."

BUILD_DIR="build/charge_dist_aor/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module "$BUILD_DIR" -o "$EXE_NAME" "$SRC_FILE"

echo "Compilation successful."

# ==========================================
# Phase 2.2: Run AoR Simulation
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
echo "case,distribution,seed_num,cutoff_m,elapsed_seconds" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${DIST_TYPE},${SEED_NUM},${RC},${ELAPSED}" >> "${OUTPUT_DIR}/timing.csv"

# ビルドディレクトリをクリーンアップ
rm -rf "$BUILD_DIR"

# 一時入力ファイルを出力ディレクトリに保存
cp "$TEMP_INPUT" "${OUTPUT_DIR}/input_${CASE_NAME}.dat"
rm -f "$TEMP_INPUT"

# ==========================================
# Phase 2.3: AoR Analysis
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
echo "Phase 2 Task completed: $CASE_NAME"
echo "Distribution  : $DIST_TYPE"
echo "Seed Number   : $SEED_NUM"
echo "Cutoff [m]    : $RC"
echo "Elapsed time  : ${ELAPSED} seconds"
echo "Output Dir    : $OUTPUT_DIR"
echo "=========================================="
date
