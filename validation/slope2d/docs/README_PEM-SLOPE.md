# PEM-SLOPE 斜面・接線KV比較メモ

## 比較のねらい
- 接線方向のKVモデル（剛性 `KT` と減衰 `c_t`）が生む過渡振動が、等価1自由度系（有効質量 m_eq=2/7 m）の減衰振動理論と一致するかを検証。
- 一般化速度 u' = v_t + r*omega と接線力 Ft の理論/数値を同一CSVに出力します。

## 実行設定の推奨
- クーロン限界によるクリップの影響を避けるため、比較時は `FRICTION_COEFF` を十分大きく設定（例: 10）。
- 例の入力 (input/slope_input.dat) 設定例:
  - `FRICTION_COEFF 10`
  - `KT 5.0e4` （任意）
  - `E_TANGENT 0.8` （任意。小さいほど減衰が大きい）
  - `TIME_STEP 1.0e-5`（既定）

## 出力CSV
- ファイル: PEM-SLOPE/data/pem_slope_trace.csv
- 追加列:
  - `u_generalized`: 数値の一般化速度 u'
  - `u_generalized_theory`: 理論の一般化変位 u
  - `vg_theory`: 理論の一般化速度 u'
  - `Ft_theory`: 理論の接線力 Ft

## 解析の観点
- アンダーダンピング条件（E_TANGENT < 1 かつ十分な KT）では、数値の Ft 振動周波数・減衰が理論の omega_n=sqrt(KT/m_eq), zeta=c_t/(2*sqrt(KT*m_eq)) と一致するか確認。
- 定常で Ft が (2/7) m g sin(theta) に近づくことを確認（剛体転がりの理論と整合）。
# 二次元PEM（KV＋定数k）斜面検証プログラム
# ビルド例:
#  gfortran -O3 -Wall -o pem_slope_kv src/pem2d_slope_kv.f90
