# 2次元斜面検証用DEMシミュレータ

## 概要

このプログラムは、個別要素法（DEM: Discrete Element Method）を用いて、斜面上での球形粒子の運動（滑り・転がり）をシミュレートし、理論解と比較検証するものです。

## 理論的背景

### 斜面上の球の運動

斜面角度 θ、摩擦係数 μ の斜面上を運動する半径 r、質量 m、慣性モーメント I = (2/5)mr² の球について：

#### 臨界摩擦係数

球が滑るか転がるかを分ける臨界摩擦係数：

```
μ_c = (2/7) tan(θ)
```

#### 滑り条件 (μ < μ_c)

粒子が滑りながら斜面を下る場合：

- 加速度: `a = g(sin(θ) - μ cos(θ))`
- 速度: `v = at`
- 角速度: `ω = 0` （滑っているため回転しない）

#### 転がり条件 (μ ≥ μ_c)

粒子が滑らずに転がる場合：

- 加速度: `a = (5/7) g sin(θ) = g sin(θ) / (1 + I/(mr²))`
- 速度: `v = at`
- 角速度: `ω = v/r`

### 接触力モデル

#### 法線方向力

```
F_n = k δ + c v_n
```

- k: ばね定数（入力値を使用）
- δ: オーバーラップ量
- c: 粘性係数（反発係数 e から計算）
- v_n: 法線方向相対速度

粘性係数の計算：
```
c = -2 ln(e) √(mk) / √(ln²(e) + π²)
```

#### 接線方向力（クーロン摩擦）

```
F_t = min(|F_t_elastic|, μ F_n) × sign(F_t_elastic)
```

- F_t_elastic: 弾性力（履歴を考慮）
- μ: 摩擦係数

滑り判定：
- |F_t_elastic| ≤ μ F_n → 転がり（弾性力が支配的）
- |F_t_elastic| > μ F_n → 滑り（摩擦力が上限）

## プログラム構成

### ファイル構成

```
TOY-PEM/
├── particle_2d_slope.f90          # メインプログラム
├── input/
│   ├── slope_input_sliding.dat    # 入力ファイル（滑り条件）
│   └── slope_input_rolling.dat    # 入力ファイル（転がり条件）
├── data/
│   └── slope_trace.csv            # 出力データ（実行後に生成）
├── plots/
│   ├── plot_slope_validation.py   # プロットスクリプト
│   ├── slope_velocity.png         # 速度比較プロット（実行後に生成）
│   ├── slope_trajectory.png       # 軌跡プロット（実行後に生成）
│   └── slope_angular_velocity.png # 角速度プロット（実行後に生成）
└── README_slope.md                # このファイル
```

### プログラムの主要部分

- **モジュール**:
  - `slope_constants_mod`: 物理定数
  - `slope_parameters_mod`: シミュレーションパラメータ
  - `particle_data_mod`: 粒子の状態変数

- **主要サブルーチン**:
  - `read_input_file`: 入力ファイル読み込み
  - `init_sub`: パラメータ初期化
  - `wcont_sub`: 斜面接触判定
  - `actf_sub`: 接触力計算
  - `nposit_leapfrog_sub`: 蛙飛び法による時間積分
  - `calculate_theoretical_velocity`: 理論解計算
  - `gfout_sub`: データ出力

## 使用方法

### 1. コンパイル

```bash
cd TOY-PEM
gfortran -o particle_2d_slope particle_2d_slope.f90
```

または最適化オプション付き：
```bash
gfortran -O3 -o particle_2d_slope particle_2d_slope.f90
```

### 2. 実行

#### 滑り条件での実行

```bash
./particle_2d_slope input/slope_input_sliding.dat
```

または引数なし（デフォルトで slope_input.dat を探す）：
```bash
./particle_2d_slope slope_input_sliding.dat
```

#### 転がり条件での実行

```bash
./particle_2d_slope input/slope_input_rolling.dat
```

### 3. 結果の可視化

```bash
python3 plots/plot_slope_validation.py
```

出力される図：
- `plots/slope_velocity.png`: 速度の数値解と理論解の比較
- `plots/slope_trajectory.png`: 粒子の軌跡
- `plots/slope_angular_velocity.png`: 角速度の時系列

## 入力パラメータ

### slope_input.dat の設定項目

| パラメータ | 説明 | 単位 | 推奨値 |
|-----------|------|------|--------|
| TIME_STEP | 時間刻み | s | 1.0d-4 |
| MAX_STEPS | 最大ステップ数 | - | 50000 |
| OUTPUT_INTERVAL | 出力間隔 | ステップ | 100 |
| SLOPE_ANGLE | 斜面角度 | rad | 0.785398 (45°) |
| FRICTION_COEFF | 摩擦係数 | - | 0.2 (滑り) / 0.5 (転がり) |
| STIFFNESS_K | ばね定数 | N/m | 1.0d3 |
| RESTITUTION_COEFF | 反発係数 | - | 0.8 |
| PARTICLE_RADIUS | 粒子半径 | m | 1.0d-2 |
| PARTICLE_DENSITY | 粒子密度 | kg/m³ | 2500.0 |
| INITIAL_X | 初期x座標 | m | 1.0d-2 |
| INITIAL_Z | 初期z座標 | m | 1.0d-2 |

## 検証結果の解釈

### 精度評価

相対誤差が1%以内であれば、シミュレーションは理論解を良好に再現していると判断できます。

#### 滑り条件 (μ = 0.2, θ = 45°)

- 理論加速度: a = g(sin(45°) - 0.2×cos(45°)) ≈ 5.55 m/s²
- 期待される挙動: 粒子が滑りながら加速、角速度はほぼゼロ

#### 転がり条件 (μ = 0.5, θ = 45°)

- 理論加速度: a = (5/7)g sin(45°) ≈ 4.95 m/s²
- 期待される挙動: 粒子が転がりながら加速、ω = v/r の関係が成立

### 誤差の原因

以下の要因により、理論解との誤差が生じる可能性があります：

1. **時間刻みの大きさ**: 時間刻みが大きすぎると数値誤差が増加
2. **ばね定数の設定**: ばね定数が小さいと接触面での貫入が大きくなる
3. **初期条件**: 粒子が完全に斜面に接触するまでに時間がかかる
4. **接線剛性の設定**: 接線剛性が不適切だと滑り・転がりの判定が正確でない

## 数値積分法

本プログラムでは**蛙飛び法（Leapfrog法）**を使用しています。

### 蛙飛び法の特徴

- 時間2次精度
- シンプレクティック積分法（エネルギー保存性が良い）
- 位置と速度が半ステップずれた時刻で評価される

### 更新手順

1. **初期化**: v(0) → v(Δt/2)
2. **位置更新**: x(t+Δt) = x(t) + v(t+Δt/2) Δt
3. **力計算**: F(t+Δt) を計算
4. **速度更新**: v(t+3Δt/2) = v(t+Δt/2) + a(t+Δt) Δt

## 参考文献

1. Cundall, P. A., & Strack, O. D. (1979). A discrete numerical model for granular assemblies. *Géotechnique*, 29(1), 47-65.

2. 個別要素法の紹介 - Qiita
   - https://qiita.com/fujitagodai4/items/c936fd82247e46731289
   - https://qiita.com/fujitagodai4/items/172ea4a5a056fc90057e

3. Goldsmith, W. (2001). *Impact: the theory and physical behaviour of colliding solids*. Courier Corporation.

## トラブルシューティング

### コンパイルエラー

- gfortran がインストールされているか確認
- Fortran 90/95 以降に対応したコンパイラを使用

### 実行時エラー

- `data` ディレクトリが存在するか確認（自動作成されない場合は手動作成）
- 入力ファイルのパスが正しいか確認

### 結果が理論解と大きく異なる

1. 時間刻みを小さくする (TIME_STEP を 1.0d-5 程度に)
2. ばね定数を大きくする (STIFFNESS_K を 1.0d4 以上に)
3. 初期位置を斜面に近づける

## ライセンス

このプログラムは教育・研究目的で自由に使用できます。

## 作成者

個別要素法検証プログラム開発チーム

## 更新履歴

- 2025-10-21: 初版作成



