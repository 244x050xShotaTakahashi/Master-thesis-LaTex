#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 48:00:00
#SBATCH -o log/stdout.%A_%a.log
#SBATCH -e log/stderr.%A_%a.log
#SBATCH --array=1-9

# ==========================================
# クーロン力カットオフ半径スウィープ（アレイジョブ版）
# 
# 比較ケース:
#   1. クーロン力OFF（ベースライン）
#   --- 固定電荷（DEFAULT_CHARGE使用）---
#   2. カットオフ 0.001 m
#   3. カットオフ 0.01 m
#   4. カットオフ 0.1 m
#   5. カットオフ 1.0 m
#   --- 二峰性電荷分布（charges.dat使用）---
#   6. カットオフ 0.001 m + bimodal
#   7. カットオフ 0.01 m + bimodal
#   8. カットオフ 0.1 m + bimodal
#   9. カットオフ 1.0 m + bimodal
#
# 使用方法:
#   sbatch job_coulomb_comparison.sh
# ==========================================

# Intelコンパイラのロード
module load intel

# スタックサイズ制限解除
ulimit -s unlimited

# OpenMPスレッド数設定
export OMP_NUM_THREADS=32

date

# ==========================================
# 設定エリア
# ==========================================

# ベース入力ファイル（テンプレート）
BASE_INPUT="inputs/input_coulomb.dat"

# 二峰性電荷分布ファイル
BIMODAL_CHARGE_FILE="inputs/dist_radii_study/charges.dat"

# 安息角測定の対象斜面
REPOSE_SIDE="right"

# 容器幅
CONTAINER_WIDTH=3.0

# ==========================================
# アレイジョブごとの設定
# ==========================================

case $SLURM_ARRAY_TASK_ID in
    1)
        CASE_NAME="no_coulomb"
        ENABLE_COULOMB=0
        CUTOFF_RADIUS=0.1
        USE_BIMODAL=0
        ;;
    2)
        CASE_NAME="cutoff_0.001_fixed"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.001
        USE_BIMODAL=0
        ;;
    3)
        CASE_NAME="cutoff_0.01_fixed"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.01
        USE_BIMODAL=0
        ;;
    4)
        CASE_NAME="cutoff_0.1_fixed"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.1
        USE_BIMODAL=0
        ;;
    5)
        CASE_NAME="cutoff_1.0_fixed"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=1.0
        USE_BIMODAL=0
        ;;
    6)
        CASE_NAME="cutoff_0.001_bimodal"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.001
        USE_BIMODAL=1
        ;;
    7)
        CASE_NAME="cutoff_0.01_bimodal"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.01
        USE_BIMODAL=1
        ;;
    8)
        CASE_NAME="cutoff_0.1_bimodal"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=0.1
        USE_BIMODAL=1
        ;;
    9)
        CASE_NAME="cutoff_1.0_bimodal"
        ENABLE_COULOMB=1
        CUTOFF_RADIUS=1.0
        USE_BIMODAL=1
        ;;
    *)
        echo "Error: Invalid array task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

# ==========================================
# 入力ファイルの生成
# ==========================================

if [ ! -f "$BASE_INPUT" ]; then
    echo "Error: Base input file '$BASE_INPUT' not found."
    exit 1
fi

if [ "$USE_BIMODAL" -eq 1 ] && [ ! -f "$BIMODAL_CHARGE_FILE" ]; then
    echo "Error: Bimodal charge file '$BIMODAL_CHARGE_FILE' not found."
    exit 1
fi

# 出力ディレクトリの決定
RUN_ID="job_${SLURM_ARRAY_JOB_ID}_task_${SLURM_ARRAY_TASK_ID}"
OUTPUT_DIR="results/coulomb_comparison/${CASE_NAME}/${RUN_ID}"

# 一時入力ファイルの生成（inputs/ 以下に配置してプログラムが正しく読み込めるようにする）
TEMP_INPUT_DIR="inputs/coulomb_comparison/${CASE_NAME}"
TEMP_INPUT="${TEMP_INPUT_DIR}/input_${RUN_ID}.dat"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_INPUT_DIR"

# ベース入力ファイルをコピーして設定を上書き
cp "$BASE_INPUT" "$TEMP_INPUT"

# sed で設定を上書き
sed -i "s/^ENABLE_COULOMB_FORCE.*/ENABLE_COULOMB_FORCE        $ENABLE_COULOMB/" "$TEMP_INPUT"
sed -i "s/^COULOMB_CUTOFF.*/COULOMB_CUTOFF              $CUTOFF_RADIUS/" "$TEMP_INPUT"
sed -i "s/^COULOMB_USE_CELL.*/COULOMB_USE_CELL            1/" "$TEMP_INPUT"

# 二峰性電荷分布を使用する場合
if [ "$USE_BIMODAL" -eq 1 ]; then
    # CHARGE_DISTRIBUTION_TYPE を file に設定
    if grep -q "^CHARGE_DISTRIBUTION_TYPE" "$TEMP_INPUT"; then
        sed -i "s/^CHARGE_DISTRIBUTION_TYPE.*/CHARGE_DISTRIBUTION_TYPE    file/" "$TEMP_INPUT"
    else
        echo "CHARGE_DISTRIBUTION_TYPE    file" >> "$TEMP_INPUT"
    fi
    
    # CHARGE_LIST_FILE を設定（パスをクォートで囲む - Fortranの/終端問題を回避）
    if grep -q "^CHARGE_LIST_FILE" "$TEMP_INPUT"; then
        sed -i "s|^CHARGE_LIST_FILE.*|CHARGE_LIST_FILE            \"$BIMODAL_CHARGE_FILE\"|" "$TEMP_INPUT"
    else
        echo "CHARGE_LIST_FILE            \"$BIMODAL_CHARGE_FILE\"" >> "$TEMP_INPUT"
    fi
fi

echo "=========================================="
echo "クーロン力カットオフ半径スウィープ（アレイジョブ）"
echo "=========================================="
echo "Job Array ID    : $SLURM_ARRAY_JOB_ID"
echo "Task ID         : $SLURM_ARRAY_TASK_ID"
echo "Case Name       : $CASE_NAME"
echo "Coulomb Force   : $ENABLE_COULOMB"
echo "Cutoff Radius   : $CUTOFF_RADIUS m"
echo "Bimodal Charge  : $USE_BIMODAL"
echo "Input File      : $TEMP_INPUT"
echo "Output Dir      : $OUTPUT_DIR"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "=========================================="

# ==========================================
# Phase 1: コンパイル
# ==========================================
echo ""
echo "[Phase 1] Compiling..."

BUILD_DIR="build/job_coulomb_comparison/${CASE_NAME}/${RUN_ID}"
mkdir -p "$BUILD_DIR"

SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

ifort -qopenmp -O3 -xHost -module $BUILD_DIR -o $EXE_NAME $SRC_FILE

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi
echo "Compilation successful."

# ==========================================
# Phase 2: DEMシミュレーション実行
# ==========================================
echo ""
echo "[Phase 2] Running DEM simulation..."
echo "  Case: $CASE_NAME"

START_TIME=$(date +%s)

srun "$EXE_NAME" "$TEMP_INPUT" "$OUTPUT_DIR"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

SIM_EXIT_CODE=$?
if [ $SIM_EXIT_CODE -ne 0 ]; then
    echo "Simulation failed with exit code: $SIM_EXIT_CODE"
    rm -rf "$BUILD_DIR"
    exit 1
fi

echo "Simulation completed successfully."
echo "Elapsed time: ${ELAPSED} seconds"

# 実行時間をファイルに記録
echo "case,elapsed_seconds" > "${OUTPUT_DIR}/timing.csv"
echo "${CASE_NAME},${ELAPSED}" >> "${OUTPUT_DIR}/timing.csv"

# ビルドディレクトリの削除
rm -rf "$BUILD_DIR"

# 一時入力ファイルを出力ディレクトリにコピー（記録用）
cp "$TEMP_INPUT" "${OUTPUT_DIR}/input_${CASE_NAME}.dat"

# 一時入力ファイルの削除
rm -f "$TEMP_INPUT"

# ==========================================
# Phase 3: 安息角測定
# ==========================================
echo ""
echo "[Phase 3] Measuring angle of repose..."

DATA_FILE="${OUTPUT_DIR}/particles.csv"

if [ ! -f "$DATA_FILE" ]; then
    echo "Warning: Data file '$DATA_FILE' not found. Skipping repose angle analysis."
else
    ANALYSIS_SCRIPT="src/analyze_repose_angle.py"
    
    if [ -f "$ANALYSIS_SCRIPT" ]; then
        python3 $ANALYSIS_SCRIPT \
            --data "$DATA_FILE" \
            --output "$OUTPUT_DIR" \
            --name "repose_angle" \
            --side "$REPOSE_SIDE" \
            --width "$CONTAINER_WIDTH"
        
        if [ $? -eq 0 ]; then
            echo "Repose angle analysis completed."
        else
            echo "Warning: Repose angle analysis failed."
        fi
    else
        echo "Warning: Analysis script '$ANALYSIS_SCRIPT' not found."
    fi
fi

# ==========================================
# 完了
# ==========================================
echo ""
echo "=========================================="
echo "Task $SLURM_ARRAY_TASK_ID ($CASE_NAME) completed."
echo "Cutoff Radius   : $CUTOFF_RADIUS m"
echo "Bimodal Charge  : $USE_BIMODAL"
echo "Elapsed time    : ${ELAPSED} seconds"
echo "=========================================="
date


