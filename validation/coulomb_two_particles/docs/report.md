# クーロン力妥当性検証レポート

## 目的

`src/dem_valid.f90` に実装されている `coulomb_force_sub` の力学挙動を、最も単純な二粒子系で数値解と解析解を突き合わせることで検証する。

## 手順概要

1. `validation/coulomb_two_particles/src/coulomb_two_particles.f90`  
   - モジュール本体をそのまま再利用し、粒子数を 2 個に限定した軽量ドライバを実装。  
   - 時間積分は本体と同じ蛙飛び法（速度半ステップシフト）で実装し、出力時に 0.5Δt 補正後の速度を CSV (`data/coulomb_two_particles.csv`) に書き出す。

2. 解析解の整理 (`docs/derivation.md`)  
   - 相対距離 `r(t)` が `r'' = 2kq²/(m r²)` に従うこと、エネルギー積分から導かれる
     ```
     t(r) = τ [ √(r/r0) √(r/r0 - 1) + ln( √(r/r0 - 1) + √(r/r0) ) ]
     ```
     （`τ = √(m r0³ / (4kq²))`）を利用することをまとめた。

3. 比較スクリプト (`plots/compare_coulomb_two_particles.py`)  
   - CSV からタイムシリーズを読み込み、上記解析式を数値反転（単調増加性を利用した二分探索）して `r(t)` を取得。  
   - 解析・数値それぞれの距離 / 粒子位置 / 速度を比較し、最大誤差と RMS、並びに可視化プロットを生成。

## 実行方法

```bash
cd /home/b/b37581/DEM-valid/validation/coulomb_two_particles
make          # ドライバをビルド
./bin/coulomb_two_particles
python3 plots/compare_coulomb_two_particles.py
```

出力:
- `data/coulomb_two_particles.csv` … 数値積分データ
- `plots/coulomb_two_particles_comparison.png` … 距離・速度プロット
- 標準出力 … 誤差統計

## 結果サマリ

当初の `coulomb_force_sub` 実装には、力ベクトルの計算において符号の誤り（同符号電荷に対して引力が働く設定となっていた）が存在した。これを修正（反発力となるよう符号を反転）した結果、数値解は理論解と極めて良好に一致した。

修正後のスクリプト出力:

```
距離: max |Δ| = 1.57e-06 m, RMS = 9.04e-07 m
x1  : max |Δ| = 7.84e-07 m, RMS = 4.52e-07 m
x2  : max |Δ| = 7.84e-07 m, RMS = 4.52e-07 m
速度: max |Δ| = 8.65e-06 m/s, RMS = 3.97e-06 m/s
```

これは数値積分の離散化誤差の範囲内と考えられ、クーロン力計算ルーチン `coulomb_force_sub` が物理的に正しく機能していることが確認された。本修正は `src/dem_valid.f90` にも適用済みである。
