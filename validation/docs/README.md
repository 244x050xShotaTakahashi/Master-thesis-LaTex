# TOY-PEM: 簡易粒子要素法検証プログラム

PEMシミュレータの妥当性検証のための簡易1次元検証プログラム群です。

## プログラム一覧

### 1. damped_oscillator.f90
- **目的**: 減衰振動の運動方程式を蛙飛び法で解き、理論解との一致を確認
- **特徴**: 
  - マスバネダンパ系の数値積分
  - 反発係数eから粘性係数cを計算
  - 理論解（減衰振動の解析解）との比較
- **出力**: `PEM/data/damped_oscillator.csv`

### 2. wall1d_contact.f90
- **目的**: 粒子-壁間の接触力を個別要素法（DEM）でモデル化し、理論解との一致を確認
- **特徴**:
  - 粒子が右壁に衝突する1次元シミュレーション
  - **与えられたm, k, eから直接パラメータを計算**（Hertz理論は使用せず、問題の切り分けを容易にする）
  - 接触力 = 弾性力（k×δ）+ 粘性力（c×δ̇）
  - 粘性係数c: `c = -2*ln(e)*sqrt(m*k/(ln(e)^2 + π^2))`
  - オーバーラップの理論解（減衰振動）との比較
- **パラメータ**:
  - 質量 m = 1.0 kg
  - 剛性 k = 1.0 N/m
  - 反発係数 e = 0.25
  - 粘性係数 c = 0.807 N·s/m（計算値）
  - 粒子半径 radius = 1.0 m
  - 壁位置 wall_x = 10.0 m
- **出力**: `PEM/data/overlap_validation.csv`
- **可視化**: `PEM/tools/plot_overlap_validation.py`

## コンパイルと実行

```bash
cd TOY-PEM

# damped_oscillatorのコンパイルと実行
ifort -o damped_oscillator damped_oscillator.f90 -O2
./damped_oscillator

# wall1d_contactのコンパイルと実行
ifort -o wall1d_contact wall1d_contact.f90 -O2
./wall1d_contact

# 可視化（wall1d_contactの結果）
cd ../PEM
python3 tools/plot_overlap_validation.py
```

## 検証結果の評価基準

- **RMSE < 1e-6 m**: 数値解が理論解と良好に一致
- **MaxAbsErr < 1e-6 m**: 最大誤差が許容範囲内
- グラフで数値解と理論解が重なっていることを確認

## 実装の特徴

1. **簡潔性**: 各プログラムは300-500行以内に収まり、理解しやすい
2. **検証可能性**: 理論解との直接比較により、実装の正しさを確認可能
3. **PEMとの整合性**: `pem_simulator.f90`の接触力計算と同じアルゴリズムを使用
4. **蛙飛び法**: PEMシミュレータと同じ数値積分法を使用し、一貫性を確保

## 参考

- オーバーラップの理論解: 減衰振動の解析解（underdamped, critical, overdamped）
- 粘性係数の計算: `c = -2*ln(e)*sqrt(m*k/(ln(e)^2 + π^2))`
- 接触力モデル: 線形バネダンパモデル `F = k*δ + c*δ̇`
  - Hertz接触理論は使用せず、与えられたm, k, eから直接計算
  - これにより問題の切り分けが容易になり、数値積分法の検証に集中できる

