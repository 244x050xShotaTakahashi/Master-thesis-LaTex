#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=1:t=32:c=32
#SBATCH -t 96:00:00
#SBATCH -o log/stdout.%A_%a.log
#SBATCH -e log/stderr.%A_%a.log

# set -x

# Intelコンパイラのロード
module load intel

# スタックサイズ制限解除
ulimit -s unlimited

# OpenMPスレッド数設定
export OMP_NUM_THREADS=32

date

# ==========================================
# 設定エリア: 実行したい入力ファイルを列挙
# ==========================================
INPUT_FILES="inputs/input_AoR.dat"

INPUT_FILE=${INPUT_FILES}

# 入力ファイルのチェック
if [ -z "$INPUT_FILE" ]; then
    echo "Error: No input file defined"
    exit 1
fi
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# 出力ディレクトリの決定
BASENAME=$(basename "$INPUT_FILE" .dat)
OUTPUT_DIR="results/${BASENAME}"

echo "=========================================="
echo "Input File   : $INPUT_FILE"
echo "Output Dir   : $OUTPUT_DIR"
echo "=========================================="

echo "Compiling..."

# --- 変更点: コンパイル用の一時ディレクトリを作成 ---
# ジョブごとに固有のディレクトリを作成し、そこでコンパイルすることで競合を防ぐ
BUILD_DIR="build/job_AoR"
mkdir -p $BUILD_DIR

# ソースコードのパス（絶対パスまたは相対パスを適切に解決）
SRC_FILE="src/dem_valid.f90"
EXE_NAME="${BUILD_DIR}/dem_valid"

# 一時ディレクトリ内でコンパイルを実行
# -module オプションで .mod ファイルの出力先を指定することも可能ですが、
# ここではシンプルに一時ディレクトリを作業場所として扱います。
# ただし、ifortはデフォルトでカレントディレクトリに .mod を吐くため、
# -module オプションで出力先を指定するのが確実です。

ifort -qopenmp -O3 -xHost -module $BUILD_DIR -o $EXE_NAME $SRC_FILE

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "Running..."
# 実行
srun ./$EXE_NAME "$INPUT_FILE" "$OUTPUT_DIR"

# 終了後に一時ディレクトリごと削除（実行ファイルと .mod ファイルを含む）
rm -rf $BUILD_DIR

date