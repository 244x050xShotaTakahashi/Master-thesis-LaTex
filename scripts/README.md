# DEM-valid スクリプト

DEMシミュレーション結果の可視化・解析用Pythonスクリプト集です。

## 必要なライブラリ

```bash
pip install numpy pandas matplotlib
```

---

## plot_snapshot.py

DEMシミュレーション結果（`particles.csv`）を読み込み、任意の時間またはステップの粒子スナップショットをプロットします。

### 機能

- 粒子を円として描画（半径で色分け）
- コンテナ壁（左・右・下・上）を描画
- 斜面壁（`walls.dat`）を描画
- 時間またはステップで任意の時点を選択
- 画像ファイル出力または画面表示

### 基本的な使い方

```bash
# ヘルプを表示
python scripts/plot_snapshot.py --help

# 利用可能なステップの一覧を表示
python scripts/plot_snapshot.py --list

# 時間を指定してプロット（画面表示）
python scripts/plot_snapshot.py --time 0.125

# ステップを指定してプロット（画面表示）
python scripts/plot_snapshot.py --step 250000

# 画像ファイルとして保存
python scripts/plot_snapshot.py --time 0.1 --output snapshot.png
```

### オプション一覧

| オプション | 短縮形 | デフォルト | 説明 |
|-----------|--------|-----------|------|
| `--file` | `-f` | `data/particles.csv` | 入力ファイルパス |
| `--time` | `-t` | - | プロットする時間 [秒] |
| `--step` | `-s` | - | プロットするステップ番号 |
| `--walls` | `-w` | `inputs/walls.dat` | 斜面壁定義ファイルパス |
| `--width` | - | `0.5` | コンテナ幅 [m] |
| `--height` | - | `0.0` | コンテナ高さ [m] (0=上壁なし) |
| `--output` | `-o` | - | 出力画像ファイルパス (指定なし=画面表示) |
| `--color-by` | - | `radius` | 着色に使う量 (`radius`/`charge`) |
| `--dpi` | - | `150` | 保存画像の解像度[dpi] |
| `--x-min` | - | - | 表示範囲のx最小値[m]（指定すると拡大表示） |
| `--x-max` | - | - | 表示範囲のx最大値[m] |
| `--z-min` | - | - | 表示範囲のz最小値[m] |
| `--z-max` | - | - | 表示範囲のz最大値[m] |
| `--list` | `-l` | - | 利用可能なステップを一覧表示 |
| `--title` | - | - | プロットのカスタムタイトル |

### 使用例

```bash
# 基本的な使用（時間指定、画面表示）
python scripts/plot_snapshot.py --time 0.05

# 別のディレクトリの結果を読み込み
python scripts/plot_snapshot.py --file results/run1/particles.csv --step 100000

# コンテナ幅と高さを指定
python scripts/plot_snapshot.py --time 0.1 --width 0.8 --height 0.3

# 斜面壁ファイルを指定
python scripts/plot_snapshot.py --time 0.1 --walls inputs/custom_walls.dat

# PNG画像として保存
python scripts/plot_snapshot.py --step 500000 --output final_state.png

# 帯電量で着色し、拡大＆高解像度で保存
python scripts/plot_snapshot.py --step 8273861 --color-by charge --dpi 400 \
  --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 \
  --output snapshot_charge_zoom.png

# カスタムタイトル付きで出力
python scripts/plot_snapshot.py --time 0.2 --title "Simulation at t=0.2s" --output result.png
```

### 入力ファイル形式

#### particles.csv

DEMシミュレーション (`dem_valid.f90`) が出力するCSVファイル：

```csv
step,time,id,x,z,radius,vx,vz,omega,angle,mass,charge
1,5.0000000E-07,1,1.0010000E-02,1.0010000E-02,1.0000000E-02,...
```

| カラム | 説明 |
|--------|------|
| `step` | ステップ番号 |
| `time` | 時刻 [s] |
| `id` | 粒子ID |
| `x`, `z` | 粒子中心の座標 [m] |
| `radius` | 粒子半径 [m] |
| `vx`, `vz` | 速度 [m/s] |
| `omega` | 角速度 [rad/s] |
| `angle` | 回転角度 [rad] |
| `mass` | 質量 [kg] |
| `charge` | 電荷 [C] |

#### walls.dat（オプション）

斜面壁定義ファイル：

```
# コメント行（#または!で始まる）
x1 z1 x2 z2    # 壁の始点(x1,z1)から終点(x2,z2)
0.05 0.00 0.05 0.50    # 例: x=0.05の位置に垂直壁
```

### 出力例

プロットには以下が含まれます：

- **粒子**: 円として描画、半径に応じて色分け（viridisカラーマップ）
- **カラーバー**: 半径のスケール
- **コンテナ壁**: 灰色の線
- **斜面壁**: 茶色の線
- **タイトル**: ステップ番号、時刻、粒子数
- **軸ラベル**: x [m], z [m]

---

## 関連ファイル

- `src/dem_valid.f90` - DEMシミュレーション本体
- `src/animate_pem.py` - アニメーション作成スクリプト
- `src/analyze_repose_angle.py` - 安息角解析スクリプト

### 読み込みデータを指定して GIF/MP4 を作る: animate_pem.py

`particles.csv`（CSV）または `graph11.d`（旧形式）からアニメーションを生成します。

```bash
# CSV + 帯電(charge)で着色してGIF
python3 src/animate_pem.py \
  --file results/.../particles.csv \
  --walls inputs/walls.dat \
  --output pem_charge.gif \
  --format auto \
  --color-by charge --charge-mode actual \
  --frame-step 10 --max-frames 300 --fps 30

# 旧形式(graph11.d)も auto 判定で動きます（既存の位置引数呼び出しも互換で残しています）
python3 src/animate_pem.py results/.../graph11.d out.gif
```

---

## 堆積遷移を見る（動画 + 定量）

### 1) フレーム生成＋動画化: make_transition_movie.py

`particles.csv` からステップをサンプリングしてフレーム画像を生成し、任意で `ffmpeg` により mp4 に結合します。

```bash
python3 scripts/make_transition_movie.py \
  --file results/.../particles.csv --walls inputs/walls.dat --width 1.0 \
  --color-by charge --dpi 300 --nframes 240 --sampling log \
  --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 \
  --frames-dir results/.../frames_charge_zoom \
  --mp4 results/.../transition_charge_zoom.mp4 --fps 30
```

### 2) 指標の時系列: analyze_deposition_transition.py

同じようにステップをサンプリングし、安息角と自由表面（x方向ビンごとの上端高さ）を時系列化します。

```bash
python3 scripts/analyze_deposition_transition.py \
  --file results/.../particles.csv --width 1.0 --side right \
  --nsamples 240 --sampling log \
  --x-min 0.20 --x-max 0.50 --z-min 0.00 --z-max 0.06 \
  --nbins 200 --outdir results/.../transition_metrics
```

出力（例）:
- `metrics.csv`（step/time/安息角/R^2など）
- `repose_angle_timeseries.png`
- `surface_height_heatmap.png`（時間×xの自由表面高さ）
- `surface_height_profiles.png`
- `surface_height.npz`（行列データ）















