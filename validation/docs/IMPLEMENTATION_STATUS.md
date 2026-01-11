# 2次元斜面検証用DEMシミュレータ - 実装状況

## 実装完了項目

### ✅ ファイル構成
- `particle_2d_slope.f90`: メインプログラム (約460行)
- `input/slope_input_sliding.dat`: 滑り条件用入力ファイル
- `input/slope_input_rolling.dat`: 転がり条件用入力ファイル
- `input/slope_input.dat`: デフォルト入力ファイル
- `plots/plot_slope_validation.py`: プロット用Pythonスクリプト
- `README_slope.md`: 詳細な使用説明書

### ✅ プログラム構造
- モジュール構成:
  - `slope_constants_mod`: 物理定数
  - `slope_parameters_mod`: シミュレーションパラメータ
  - `particle_data_mod`: 粒子の状態変数

- サブルーチン:
  - `read_input_file`: 入力ファイル読み込み
  - `init_sub`: パラメータ初期化
  - `wcont_sub`: 斜面接触判定
  - `actf_sub`: 接触力計算
  - `nposit_leapfrog_sub`: 蛙飛び法による時間積分
  - `calculate_theoretical_velocity`: 理論解計算
  - `gfout_sub`: データ出力

### ✅ 理論的実装
- 斜面接触判定（法線ベクトルによる距離計算）
- 法線方向接触力（ばね-ダンパーモデル）
- 接線方向接触力（クーロン摩擦）
- 蛙飛び法による時間積分
- 理論解の計算（滑り・転がり条件の判定）

## ⚠️ 既知の問題

### 数値的不安定性
**症状**: シミュレーション実行後、速度が指数関数的に発散し、NaNまたはInfinityになる

**発生タイミング**: 粒子が斜面に接触した直後（t≈0.02s）

**推定原因**:
1. 接触力計算のロジックに誤りがある
2. 力の向きの符号が不適切
3. 接線方向力の履歴管理が不適切
4. 粘性係数または剛性の設定が不適切

**試行した対策** (いずれも効果なし):
- 時間刻みの縮小 (1e-4 → 1e-5)
- ばね定数の調整 (1e3 → 1e2)
- 初期位置の調整（斜面から離す）
- 接触力計算の単純化（履歴なしモデル）
- 法線方向相対速度の符号修正

## 🔧 必要な修正

### 高優先度

1. **接触力計算の検証**
   - `actf_sub`の法線・接線方向力の計算ロジックを再確認
   - 力の向き（符号）の整合性チェック
   - `pem_simulator.f90`の`actf_sub`との詳細な比較

2. **デバッグ出力の追加**
   ```fortran
   ! actf_sub内に追加
   if (it_step < 100) then
       write(*,*) 'DEBUG: overlap=', overlap, ' F_n=', force_normal, ' F_t=', force_tangent
       write(*,*) '       v_n=', v_normal, ' v_t=', v_tangent
   end if
   ```

3. **参照実装との比較**
   - `particle_1dcollision.f90`は安定動作している
   - 同じ接触力計算ロジックを2D版に正確に移植

### 中優先度

4. **接線方向力モデルの改善**
   - 現在: 粘性のみ（履歴なし）
   - 改善案: 弾性+粘性（ただし履歴管理を慎重に）
   - `pem_simulator.f90`のロジックを参考

5. **パラメータの最適化**
   - 安定条件: `dt < 0.1 * sqrt(m/k)`
   - 現在のパラメータでの安定条件を計算
   - レイリー時間刻みの確認

### 低優先度

6. **数値積分法の検証**
   - 蛙飛び法の実装が正しいか確認
   - オイラー法での動作確認（単純化のため）

7. **可視化の改善**
   - リアルタイムプロット機能
   - アニメーション出力

## 📋 今後の実装手順

1. **段階的デバッグ**
   ```bash
   # Step 1: 重力のみ（接触なし）で落下を確認
   # Step 2: 法線方向力のみ（接線力なし）で反発を確認
   # Step 3: 接線方向力を追加
   ```

2. **簡略版の作成**
   - まず1次元（垂直落下）で動作確認
   - 次に2次元（斜面）に拡張

3. **理論値との比較**
   - 各ステップで理論値と比較
   - 誤差が許容範囲か確認

## 参考

### 動作する類似プログラム
- `particle_1dcollision.f90`: 1次元粒子衝突（安定動作）
- `PEM/src/pem_simulator.f90`: 2次元多粒子DEM（安定動作）

### 参考文献
- https://qiita.com/fujitagodai4/items/c936fd82247e46731289
- https://qiita.com/fujitagodai4/items/172ea4a5a056fc90057e

## コンパイルと実行

```bash
# コンパイル
cd TOY-PEM
ifort -O2 -o particle_2d_slope particle_2d_slope.f90

# 実行（現在は数値的に不安定）
./particle_2d_slope input/slope_input_sliding.dat

# データ確認
head -20 data/slope_trace.csv
```

## 作成日

2025-10-21

## 更新履歴

- 2025-10-21: 初版作成、数値的不安定性の問題を記録



