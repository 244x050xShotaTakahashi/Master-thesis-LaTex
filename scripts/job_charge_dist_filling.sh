#!/bin/bash
#SBATCH -p gr10451a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 168:00:00
#SBATCH -o log/charge_dist_filling_stdout.%A_%a.log
#SBATCH -e log/charge_dist_filling_stderr.%A_%a.log
#SBATCH --array=1-3

# ==========================================
# Phase 1: 帯電分布ごとの充填状態ファイル生成
#
# - カットオフ rc=0.05m（大きなカットオフ）で充填
# - 3種類の帯電分布: bimodal, normal, uniform
# - AUTO_WALL_WITHDRAW=0 で壁引き抜きなし、静止で終了
# - filled_particles.dat を生成して Phase 2 で使用
#
# 使い方:
#   sbatch job_charge_dist_filling.sh
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

# 帯電分布の設定
case ${SLURM_ARRAY_TASK_ID:-1} in
  1) 
    DIST_TYPE="bimodal"
    CASE_NAME="bimodal"
    ;;
  2) 
    DIST_TYPE="normal"
    CASE_NAME="normal"
    ;;
  3) 
    DIST_TYPE="uniform"
    CASE_NAME="uniform"
    ;;
  *) 
    echo "Error: Invalid array task ID: ${SLURM_ARRAY_TASK_ID:-unset}"
    exit 1 
    ;;
esac

RUN_ID="job_${SLURM_ARRAY_JOB_ID:-local}_task_${SLURM_ARRAY_TASK_ID:-0}"
OUTPUT_DIR="${RESULTS_BASE}/${CASE_NAME}"

TEMP_INPUT_DIR="inputs/charge_dist_study/${CASE_NAME}"
TEMP_INPUT="${TEMP_INPUT_DIR}/input_${RUN_ID}.dat"

mkdir -p "$OUTPUT_DIR" "$TEMP_INPUT_DIR" log

cp "$BASE_INPUT" "$TEMP_INPUT"

# 電荷分布タイプを上書き
sed -i "s/^CHARGE_DISTRIBUTION_TYPE.*/CHARGE_DISTRIBUTION_TYPE    ${DIST_TYPE}/" "$TEMP_INPUT"

# CONTAINER_WIDTH を input から取得
CONTAINER_WIDTH=$(awk '/^CONTAINER_WIDTH/{print $2; exit}' "$TEMP_INPUT")
if [ -z "${CONTAINER_WIDTH}" ]; then
  echo "Warning: failed to read CONTAINER_WIDTH from input; fallback to 2.0"
  CONTAINER_WIDTH="2.0"
fi

echo "=========================================="
echo "Phase 1: Charge Distribution Filling"
echo "=========================================="
echo "Distribution    : $DIST_TYPE"
echo "Case Name       : $CASE_NAME"
echo "Input File      : $TEMP_INPUT"
echo "Output Dir      : $OUTPUT_DIR"
echo "Container Width : $CONTAINER_WIDTH"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "=========================================="

# ==========================================
# Phase 1: Compile
# ==========================================
echo ""
echo "[Phase 1.1] Compiling..."

BUILD_DIR="build/charge_dist_filling/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module "$BUILD_DIR" -o "$EXE_NAME" "$SRC_FILE"

echo "Compilation successful."

# ==========================================
# Phase 1: Run Filling Simulation
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
echo "case,distribution,elapsed_seconds" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${DIST_TYPE},${ELAPSED}" >> "${OUTPUT_DIR}/timing.csv"

# ==========================================
# Phase 1: Post-processing
# ==========================================
echo ""
echo "[Phase 1.3] Post-processing..."

# filled_particles.dat を出力ディレクトリにコピー（コードが生成する場所から）
# コードは use_explicit_positions=false の場合、inputs/filled_particles.dat に保存
# use_explicit_positions=true の場合、output_dir/filled_particles.dat に保存
# 今回は乱数生成なので inputs/filled_particles.dat から取得

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
echo "Phase 1 Task completed: $CASE_NAME"
echo "Distribution  : $DIST_TYPE"
echo "Elapsed time  : ${ELAPSED} seconds"
echo "Output Dir    : $OUTPUT_DIR"
echo "=========================================="
date













