# 可視化生成スクリプト

`charge_distribution_study/phase2_aor` の各ケースに対して、参照ケース（`aor_cutoff_sweep/rc_0.02/job_6292283_task_1/`）と同様の可視化を生成するためのスクリプトです。

## 生成される可視化

各ケースディレクトリに以下のファイルが生成されます：

1. **分布プロット**
   - `distribution_step*.png`: 半径・電荷の分布ヒストグラムとKDE

2. **スナップショットプロット**
   - `snapshot.png`: 最終ステップの粒子配置（半径で着色）
   - `snapshot_charge.png`: 最終ステップの粒子配置（電荷で着色）
   - `snapshot_charge_zoom.png`: 電荷で着色した拡大表示（x=[0.20, 0.50], z=[0.00, 0.06]）

3. **トランジションメトリクス** (`transition_metrics/` ディレクトリ)
   - `metrics.csv`: 各ステップの安息角などのメトリクス
   - `repose_angle_timeseries.png`: 安息角の時系列プロット
   - `surface_height_heatmap.png`: 自由表面高さのヒートマップ
   - `surface_height_profiles.png`: 代表時刻の表面プロファイル
   - `surface_height.npz`: 表面高さデータ（後処理用）

## 使用方法

### すべてのケースに対して可視化を生成

```bash
cd /LARGE0/gr20001/b37581/DEM-valid
bash scripts/generate_visualizations_for_cases.sh
```

### 個別のケースに対して可視化を生成

```bash
cd /LARGE0/gr20001/b37581/DEM-valid
CASE_DIR="results/charge_distribution_study/phase2_aor/uniform_rc_0.02"

# 1. 分布データファイルの生成
python3 scripts/extract_distribution_from_csv.py \
    --file "$CASE_DIR/particles.csv" \
    --output-dir "$CASE_DIR"

# 2. 分布プロットの生成
FINAL_STEP=$(python3 -c "import pandas as pd; df = pd.read_csv('$CASE_DIR/particles.csv', skipinitialspace=True); df.columns = df.columns.str.strip(); print(int(df['step'].max()))")
python3 inputs/dist_radii_study/plot_distribution.py \
    --radius-file "$CASE_DIR/radii_step${FINAL_STEP}.dat" \
    --charge-file "$CASE_DIR/charges_step${FINAL_STEP}.dat" \
    --output "$CASE_DIR/distribution_step${FINAL_STEP}.png" \
    --no-show

# 3. スナップショットプロットの生成
python3 scripts/plot_snapshot.py \
    --file "$CASE_DIR/particles.csv" \
    --step "$FINAL_STEP" \
    --width 2.0 \
    --walls inputs/walls.dat \
    --output "$CASE_DIR/snapshot.png" \
    --color-by radius

python3 scripts/plot_snapshot.py \
    --file "$CASE_DIR/particles.csv" \
    --step "$FINAL_STEP" \
    --width 2.0 \
    --walls inputs/walls.dat \
    --output "$CASE_DIR/snapshot_charge.png" \
    --color-by charge

python3 scripts/plot_snapshot.py \
    --file "$CASE_DIR/particles.csv" \
    --step "$FINAL_STEP" \
    --width 2.0 \
    --walls inputs/walls.dat \
    --output "$CASE_DIR/snapshot_charge_zoom.png" \
    --dpi 400 \
    --color-by charge \
    --x-min 0.20 --x-max 0.50 \
    --z-min 0.00 --z-max 0.06

# 4. トランジションメトリクスの生成
python3 scripts/analyze_deposition_transition.py \
    --file "$CASE_DIR/particles.csv" \
    --width 2.0 \
    --outdir "$CASE_DIR/transition_metrics" \
    --nsamples 240 \
    --sampling log \
    --x-min 0.20 --x-max 0.50 \
    --z-min 0.00 --z-max 0.06 \
    --nbins 200 \
    --side right
```

## スクリプトの説明

### `extract_distribution_from_csv.py`

`particles.csv`から最終ステップの電荷・半径データを抽出して、`charges_step*.dat`と`radii_step*.dat`ファイルを生成します。

**使用方法:**
```bash
python3 scripts/extract_distribution_from_csv.py \
    --file particles.csv \
    --output-dir .
```

### `generate_visualizations_for_cases.sh`

すべてのケースに対して上記の可視化を一括生成するバッチスクリプトです。

**設定:**
- コンテナ幅: 2.0 m
- 解析範囲: x=[0.20, 0.50], z=[0.00, 0.06]
- サンプリング: 240サンプル（対数スケール）

## 注意事項

- `particles.csv`に`charge`列が含まれている必要があります
- 各ケースディレクトリに`particles.csv`が存在する必要があります
- `inputs/walls.dat`が存在する必要があります



