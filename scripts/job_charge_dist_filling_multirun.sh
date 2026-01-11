#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 168:00:00
#SBATCH -o log/charge_dist_filling_multirun_stdout.%A_%a.log
#SBATCH -e log/charge_dist_filling_multirun_stderr.%A_%a.log
#SBATCH --array=1-12

# ==========================================
# Phase 1 追加シード: 帯電分布ごとの充填状態ファイル生成
#
# - シード2-5 の充填を実行（シード1は既存結果を使用）
# - 4シード × 3分布 = 12ケース
# - カットオフ rc=0.05m（大きなカットオフ）で充填
#
# 使い方:
#   sbatch job_charge_dist_filling_multirun.sh
# ==========================================

set -euo pipefail

module load intel
ulimit -s unlimited
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-32}

date

BASE_INPUT="inputs/charge_dist_study/input_filling_template.dat"
RESULTS_BASE="results/charge_distribution_study/phase1_filling"

if [ ! -f "$BASE_INPUT" ]; then
  echo "Error: Base input file '$BASE_INPUT' not found."
  exit 1
fi

# 乱数シードの配列（シード2-5、シード1は既存）
SEEDS=("584288" "584289" "584290" "584291")
SEED_IDS=("seed2" "seed3" "seed4" "seed5")
NUM_SEEDS=${#SEEDS[@]}

# 分布の配列
DISTRIBUTIONS=("bimodal" "normal" "uniform")
NUM_DISTS=${#DISTRIBUTIONS[@]}

# Array ID からシードと分布を決定
# Array ID 1-3: seed2 × (bimodal, normal, uniform)
# Array ID 4-6: seed3 × (bimodal, normal, uniform)
# ...
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
SEED_IDX=$(( (TASK_ID - 1) / NUM_DISTS ))
DIST_IDX=$(( (TASK_ID - 1) % NUM_DISTS ))

if [ $SEED_IDX -ge $NUM_SEEDS ] || [ $DIST_IDX -ge $NUM_DISTS ]; then
  echo "Error: Invalid array task ID: $TASK_ID"
  exit 1
fi

SEED="${SEEDS[$SEED_IDX]}"
SEED_ID="${SEED_IDS[$SEED_IDX]}"
DIST_TYPE="${DISTRIBUTIONS[$DIST_IDX]}"
CASE_NAME="${DIST_TYPE}_${SEED_ID}"

RUN_ID="job_${SLURM_ARRAY_JOB_ID:-local}_task_${SLURM_ARRAY_TASK_ID:-0}"
OUTPUT_DIR="${RESULTS_BASE}/${CASE_NAME}"

TEMP_INPUT_DIR="inputs/charge_dist_study/multirun/${CASE_NAME}"
TEMP_INPUT="${TEMP_INPUT_DIR}/input_${RUN_ID}.dat"

mkdir -p "$OUTPUT_DIR" "$TEMP_INPUT_DIR" log

cp "$BASE_INPUT" "$TEMP_INPUT"

# 電荷分布タイプを上書き
sed -i "s/^CHARGE_DISTRIBUTION_TYPE.*/CHARGE_DISTRIBUTION_TYPE    ${DIST_TYPE}/" "$TEMP_INPUT"

# 乱数シードを上書き
sed -i "s/^RANDOM_SEED.*/RANDOM_SEED                 ${SEED}/" "$TEMP_INPUT"

# CONTAINER_WIDTH を input から取得
CONTAINER_WIDTH=$(awk '/^CONTAINER_WIDTH/{print $2; exit}' "$TEMP_INPUT")
if [ -z "${CONTAINER_WIDTH}" ]; then
  echo "Warning: failed to read CONTAINER_WIDTH from input; fallback to 2.0"
  CONTAINER_WIDTH="2.0"
fi

echo "=========================================="
echo "Phase 1 Multirun: Charge Distribution Filling"
echo "=========================================="
echo "Seed ID         : $SEED_ID"
echo "Random Seed     : $SEED"
echo "Distribution    : $DIST_TYPE"
echo "Case Name       : $CASE_NAME"
echo "Input File      : $TEMP_INPUT"
echo "Output Dir      : $OUTPUT_DIR"
echo "Container Width : $CONTAINER_WIDTH"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "=========================================="

# ==========================================
# Compile
# ==========================================
echo ""
echo "[Phase 1.1] Compiling..."

BUILD_DIR="build/charge_dist_filling_multirun/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module "$BUILD_DIR" -o "$EXE_NAME" "$SRC_FILE"

echo "Compilation successful."

# ==========================================
# Run Filling Simulation
# ==========================================
echo ""
echo "[Phase 1.2] Running filling simulation..."

START_TIME=$(date +%s)
srun "$EXE_NAME" "$TEMP_INPUT" "$OUTPUT_DIR"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "Simulation completed successfully."
echo "Elapsed time: ${ELAPSED} seconds"

# タイミング情報を保存
echo "case,distribution,seed_id,random_seed,elapsed_seconds" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${DIST_TYPE},${SEED_ID},${SEED},${ELAPSED}" >> "${OUTPUT_DIR}/timing.csv"

# ==========================================
# Post-processing
# ==========================================
echo ""
echo "[Phase 1.3] Post-processing..."

# filled_particles.dat を出力ディレクトリにコピー
if [ -f "inputs/filled_particles.dat" ]; then
  cp "inputs/filled_particles.dat" "${OUTPUT_DIR}/filled_particles.dat"
  echo "Copied filled_particles.dat to ${OUTPUT_DIR}/"
elif [ -f "${OUTPUT_DIR}/filled_particles.dat" ]; then
  echo "filled_particles.dat already in ${OUTPUT_DIR}/"
else
  echo "Warning: filled_particles.dat not found!"
fi

# ビルドディレクトリをクリーンアップ
rm -rf "$BUILD_DIR"

# 一時入力ファイルを出力ディレクトリに保存
cp "$TEMP_INPUT" "${OUTPUT_DIR}/input_${CASE_NAME}.dat"
rm -f "$TEMP_INPUT"

echo ""
echo "=========================================="
echo "Phase 1 Multirun Task completed: $CASE_NAME"
echo "Seed ID       : $SEED_ID"
echo "Random Seed   : $SEED"
echo "Distribution  : $DIST_TYPE"
echo "Elapsed time  : ${ELAPSED} seconds"
echo "Output Dir    : $OUTPUT_DIR"
echo "=========================================="
date







