# DEM Validation Tests Organization

このディレクトリには、個別要素法(DEM)の計算妥当性を検証するためのFortranトイモデルと、その結果データおよび可視化用Pythonスクリプトが整理されています。

## ディレクトリ構造

```
validation/
├── docs/                    # 全体ドキュメント
│   ├── README.md            # TOY-PEMメインREADME
│   └── IMPLEMENTATION_STATUS.md  # 実装状況
├── damped_oscillator/       # 減衰振動の検証
│   ├── src/                 # Fortranソースコード
│   ├── bin/                 # コンパイル済みバイナリ
│   ├── data/                # 出力データとグラフ
│   └── plots/               # Pythonプロットスクリプト
├── wall1d/                  # 壁-粒子接触力の検証
│   ├── src/
│   ├── bin/
│   ├── data/
│   └── plots/
├── particle1d/              # 1次元粒子衝突の検証
│   ├── src/
│   ├── bin/
│   ├── data/
│   └── plots/
└── slope2d/                 # 2次元斜面の検証
    ├── src/                 # TOY-PEMとPEM-SLOPEのソースコード
    ├── bin/                 # 複数のバイナリ
    ├── data/                # 両プロジェクトの出力データ
    ├── input/               # 入力ファイル
    ├── plots/               # プロットスクリプトとグラフ
    └── docs/                # slope2d固有のドキュメント
```

## 各検証テストの説明

### 1. damped_oscillator（減衰振動）
- **目的**: 蛙飛び法による数値積分の妥当性検証
- **検証内容**: マスバネダンパ系の理論解との比較
- **主要ファイル**:
  - `damped_oscillator.f90`: 基本実装
  - `damped_oscillator_rk4.f90`: Runge-Kutta法実装

### 2. wall1d（壁接触）
- **目的**: 粒子-壁間の接触力モデルの検証
- **検証内容**: オーバーラップの理論解との比較
- **主要ファイル**:
  - `wall1d_contact.f90`: 1次元壁接触シミュレーション

### 3. particle1d（粒子衝突）
- **目的**: 粒子間衝突の接触力モデルの検証
- **検証内容**: 1次元での粒子衝突シミュレーション
- **主要ファイル**:
  - `particle_1dcollision.f90`: 粒子衝突シミュレーション
  - `animate_particle_collision.py`: アニメーション生成

### 4. slope2d（斜面シミュレーション）
- **目的**: 2次元での斜面上の粒子運動の検証
- **検証内容**: 摩擦・転がり抵抗を含む複雑な力学系の検証
- **主要ファイル**:
  - `particle_2d_slope.f90`: TOY-PEM版
  - `pem2d_slope_kv.f90`: PEM-SLOPE版
- **注意**: 実装状況ドキュメント参照（数値安定性の課題あり）

## 使用方法

各カテゴリのディレクトリに移動して:

```bash
# コンパイル (例: damped_oscillator)
cd validation/damped_oscillator/src
ifort -O2 -o ../bin/damped_oscillator damped_oscillator.f90

# 実行
cd ../bin
./damped_oscillator

# プロット
cd ../plots
python3 plot_damped_oscillator.py
```

## 元のフォルダ

- **TOY-PEM**: `/home/b/b37581/TOY-PEM/`
- **PEM-SLOPE**: `/home/b/b37581/PEM-SLOPE/`

元のフォルダは保持されています（コピーされたため削除されていません）。

## 整理日

2025-11-10
