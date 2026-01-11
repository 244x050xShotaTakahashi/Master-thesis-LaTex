#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 24:00:00
#SBATCH -o log/stdout.%j.log
#SBATCH -e log/stderr.%j.log

# ==========================================
# クーロン力（静電気力）シミュレーション ジョブスクリプト
# 
# 処理フロー:
#   1. 粒子を容器内に充填・堆積（静止まで）
#   2. 指定ステップで片側壁を引き抜き
#   3. 粒子が崩落して再び静止するまで計算
#   4. 安息角を自動測定
#
# 使用方法:
#   sbatch job_coulomb.sh
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

# 入力ファイル
INPUT_FILE="inputs/input_coulomb.dat"

# 安息角測定の対象斜面 (left=左壁引き抜き後, right=右壁引き抜き後)
REPOSE_SIDE="right"

# 容器幅 (安息角測定時に使用、入力ファイルのCONTAINER_WIDTHと一致させる)
CONTAINER_WIDTH=3.0

# ==========================================
# 入力ファイルのチェック
# ==========================================

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# 出力ディレクトリの決定
BASENAME=$(basename "$INPUT_FILE" .dat)
RUN_ID="job_${SLURM_JOB_ID}"
OUTPUT_DIR_BASE="results/coulomb_study/${BASENAME}"
OUTPUT_DIR="${OUTPUT_DIR_BASE}/${RUN_ID}"

echo "=========================================="
echo "クーロン力シミュレーション"
echo "=========================================="
echo "Job ID          : $SLURM_JOB_ID"
echo "Input File      : $INPUT_FILE"
echo "Output Dir      : $OUTPUT_DIR"
echo "Repose Side     : $REPOSE_SIDE"
echo "Container Width : $CONTAINER_WIDTH m"
echo "OMP Threads     : $OMP_NUM_THREADS"
echo "=========================================="

# 出力ディレクトリの作成
mkdir -p "$OUTPUT_DIR"

# ==========================================
# Phase 1: コンパイル
# ==========================================
echo ""
echo "[Phase 1] Compiling..."

BUILD_DIR_ROOT="build/job_coulomb_study"
BUILD_DIR="${BUILD_DIR_ROOT}/${BASENAME}/${RUN_ID}"
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
echo "[Phase 2] Running DEM simulation with Coulomb force..."
echo "  - 充填・堆積 → 壁引き抜き → 崩落・再静止"

srun "$EXE_NAME" "$INPUT_FILE" "$OUTPUT_DIR"

SIM_EXIT_CODE=$?
if [ $SIM_EXIT_CODE -ne 0 ]; then
    echo "Simulation failed with exit code: $SIM_EXIT_CODE"
    rm -rf "$BUILD_DIR"
    exit 1
fi
echo "Simulation completed successfully."

# ビルドディレクトリの削除
rm -rf "$BUILD_DIR"

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
        echo "Running: python $ANALYSIS_SCRIPT --data $DATA_FILE --output $OUTPUT_DIR --name repose_angle --side $REPOSE_SIDE --width $CONTAINER_WIDTH"
        
        python3 $ANALYSIS_SCRIPT \
            --data "$DATA_FILE" \
            --output "$OUTPUT_DIR" \
            --name "repose_angle" \
            --side "$REPOSE_SIDE" \
            --width "$CONTAINER_WIDTH"
        
        if [ $? -eq 0 ]; then
            echo "Repose angle analysis completed."
            echo ""
            echo "Results saved to:"
            echo "  - ${OUTPUT_DIR}/repose_angle_results.csv"
            echo "  - ${OUTPUT_DIR}/repose_angle_repose_angle.png"
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
echo "All tasks completed."
echo "=========================================="
date


