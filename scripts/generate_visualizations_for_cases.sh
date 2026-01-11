#!/bin/bash
# charge_distribution_study/phase2_aor の各ケースに対して可視化を生成するバッチスクリプト
#
# 使用方法:
#   bash scripts/generate_visualizations_for_cases.sh

set -e  # エラー時に停止

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="$REPO_ROOT/results/charge_distribution_study/phase2_aor"

# コンテナ幅（すべてのケースで2.0m）
CONTAINER_WIDTH=2.0

# 解析範囲（参照ケースと同じ設定）
X_MIN=0.20
X_MAX=0.50
Z_MIN=0.00
Z_MAX=0.06

# walls.datファイルのパス
WALLS_FILE="$REPO_ROOT/inputs/walls.dat"

echo "=========================================="
echo "可視化生成バッチスクリプト"
echo "=========================================="
echo "結果ディレクトリ: $RESULTS_DIR"
echo "コンテナ幅: $CONTAINER_WIDTH m"
echo "解析範囲: x=[$X_MIN, $X_MAX], z=[$Z_MIN, $Z_MAX]"
echo ""

# 各ケースディレクトリを処理
for case_dir in "$RESULTS_DIR"/*/; do
    if [ ! -d "$case_dir" ]; then
        continue
    fi
    
    case_name=$(basename "$case_dir")
    echo "----------------------------------------"
    echo "処理中: $case_name"
    echo "----------------------------------------"
    
    particles_file="$case_dir/particles.csv"
    
    # particles.csvの存在確認
    if [ ! -f "$particles_file" ]; then
        echo "警告: particles.csvが見つかりません。スキップします: $particles_file"
        continue
    fi
    
    cd "$case_dir"
    
    # 1. 分布データファイルの生成（charges_step*.dat, radii_step*.dat）
    echo "[1/4] 分布データファイルを生成中..."
    python3 "$SCRIPT_DIR/extract_distribution_from_csv.py" \
        --file "particles.csv" \
        --output-dir "." || {
        echo "エラー: 分布データファイルの生成に失敗しました"
        continue
    }
    
    # 最終ステップを取得（extract_distribution_from_csv.pyの出力から）
    final_step=$(python3 -c "
import pandas as pd
df = pd.read_csv('particles.csv', skipinitialspace=True)
df.columns = df.columns.str.strip()
print(int(df['step'].max()))
")
    
    charges_file="charges_step${final_step}.dat"
    radii_file="radii_step${final_step}.dat"
    
    # 2. 分布プロットの生成
    echo "[2/4] 分布プロットを生成中..."
    if [ -f "$charges_file" ] && [ -f "$radii_file" ]; then
        python3 "$REPO_ROOT/inputs/dist_radii_study/plot_distribution.py" \
            --radius-file "$radii_file" \
            --charge-file "$charges_file" \
            --output "distribution_step${final_step}.png" \
            --no-show \
            --dpi 150 || {
            echo "警告: 分布プロットの生成に失敗しました"
        }
    else
        echo "警告: 分布データファイルが見つかりません。分布プロットをスキップします"
    fi
    
    # 3. スナップショットプロットの生成
    echo "[3/4] スナップショットプロットを生成中..."
    
    # 通常のスナップショット（半径で着色）
    python3 "$SCRIPT_DIR/plot_snapshot.py" \
        --file "particles.csv" \
        --step "$final_step" \
        --width "$CONTAINER_WIDTH" \
        --walls "$WALLS_FILE" \
        --output "snapshot.png" \
        --dpi 150 \
        --color-by radius || {
        echo "警告: スナップショット（半径）の生成に失敗しました"
    }
    
    # 電荷で着色したスナップショット
    if [ -f "$charges_file" ]; then
        python3 "$SCRIPT_DIR/plot_snapshot.py" \
            --file "particles.csv" \
            --step "$final_step" \
            --width "$CONTAINER_WIDTH" \
            --walls "$WALLS_FILE" \
            --output "snapshot_charge.png" \
            --dpi 150 \
            --color-by charge || {
            echo "警告: スナップショット（電荷）の生成に失敗しました"
        }
        
        # 拡大表示（参照ケースと同じ範囲）
        python3 "$SCRIPT_DIR/plot_snapshot.py" \
            --file "particles.csv" \
            --step "$final_step" \
            --width "$CONTAINER_WIDTH" \
            --walls "$WALLS_FILE" \
            --output "snapshot_charge_zoom.png" \
            --dpi 400 \
            --color-by charge \
            --x-min "$X_MIN" \
            --x-max "$X_MAX" \
            --z-min "$Z_MIN" \
            --z-max "$Z_MAX" || {
            echo "警告: スナップショット（電荷、拡大）の生成に失敗しました"
        }
    fi
    
    # 4. トランジションメトリクスの生成
    echo "[4/4] トランジションメトリクスを生成中..."
    transition_dir="transition_metrics"
    mkdir -p "$transition_dir"
    
    python3 "$SCRIPT_DIR/analyze_deposition_transition.py" \
        --file "particles.csv" \
        --width "$CONTAINER_WIDTH" \
        --outdir "$transition_dir" \
        --nsamples 240 \
        --sampling log \
        --x-min "$X_MIN" \
        --x-max "$X_MAX" \
        --z-min "$Z_MIN" \
        --z-max "$Z_MAX" \
        --nbins 200 \
        --side right || {
        echo "警告: トランジションメトリクスの生成に失敗しました"
    }
    
    echo "完了: $case_name"
    echo ""
done

echo "=========================================="
echo "すべてのケースの処理が完了しました"
echo "=========================================="

