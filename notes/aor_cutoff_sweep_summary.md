### 目的（このチャットの要点）
- **安息角（AoR）計測の計算時間を短縮**するため、クーロン力のカットオフ半径 \(r_c\) を **短縮ランで収束判定して最小化**する。
- 許容基準: **AoR差 ±0.5°以内**（同一シード・同一初期条件で比較）。

### 重要な前提（実装/解析の仕様）
- `src/analyze_repose_angle.py` は `particles.csv` の **最終stepのみ**を読み、表面粒子→回帰でAoRを計算する。
  - したがって、**終了時点の最終状態が `particles.csv` に必ず出力されていること**が重要。
- DEMコード `src/dem_valid.f90` には **AUTO_WALL_WITHDRAW（充填静止→引き抜き→再静止で終了）**が既にあるため、短縮ランが組みやすい。

### 実行手順（rc掃引 → 集計 → 推奨rc反映）

#### 1) rc掃引（短縮ラン / AUTO_WALL_WITHDRAW）

```bash
cd /LARGE0/gr20001/b37581/DEM-valid
sbatch job_aor_cutoff_sweep.sh
```

- 掃引点: **0.02 / 0.03 / 0.04 / 0.05 m**（`job_aor_cutoff_sweep.sh` の配列ジョブ）
- 出力先: `results/aor_cutoff_sweep/rc_0.0X/job_${ARRAY_JOB_ID}_task_${TASK_ID}/`

#### 2) 結果集計（AoR差±0.5°で最小rcを推奨）

```bash
cd /LARGE0/gr20001/b37581/DEM-valid
python3 src/summarize_cutoff_sweep.py \
  --root results/aor_cutoff_sweep \
  --out results/aor_cutoff_sweep/sweep_summary.csv \
  --ref-rc 0.05 \
  --tol-deg 0.5
```

- 参照: `ref_rc=0.05` の AoR（ΔAoR をこれ基準で計算）
- 出力: `results/aor_cutoff_sweep/sweep_summary.csv`

#### 3) 推奨rcを input に反映（既定は dry-run）

```bash
cd /LARGE0/gr20001/b37581/DEM-valid
bash postprocess_cutoff_sweep.sh results/aor_cutoff_sweep inputs/input_coulomb.dat 0.05 0.5
```

実際に更新する場合（`APPLY=1`）:
```bash
cd /LARGE0/gr20001/b37581/DEM-valid
APPLY=1 bash postprocess_cutoff_sweep.sh results/aor_cutoff_sweep inputs/input_coulomb.dat 0.05 0.5
```

### 変更点と理由（何を変えたか / なぜ必要か）

#### 1) 短縮掃引用 input を追加
- 追加: `inputs/input_coulomb_short.dat`
- 内容:
  - `AUTO_WALL_WITHDRAW 1`（短縮ラン）
  - `CHARGE_DISTRIBUTION_TYPE bimodal`（コード内で電荷分布を生成し、ファイルを読まない）
  - `COULOMB_CUTOFF` はジョブ側で上書きして掃引
  - `PROFILING_SAMPLE_INTERVAL 0`（後述の timing_report.csv を確実に出す）

**理由**: 同一条件の短縮ランを簡単に再現でき、rcだけを変えてAoR収束判定できるようにするため。

#### 2) 終了時に `particles.csv` の最終状態を必ず出力
- 変更: `src/dem_valid.f90`
- 変更内容:
  - シミュレーションループを抜けた直後に **`gfout_sub(it_step, current_time, rmax_particle_radius)` を追加**

**理由**: AoR解析が `particles.csv` の **最終step** だけを読むため、`OUTPUT_INTERVAL` を粗くしても最終状態が欠落しないようにする（I/O削減が可能）。

#### 3) `timing_report.csv` が「サンプリング無効化のタイミング」で出ない問題を修正
- 変更: `src/dem_valid.f90`（`profiler_write_csv`）
- 変更内容:
  - 従来: `profiling_enabled=false` だと **即return**
  - 修正: `profile_entry_count>0` なら書き出し可能に変更（終了時点でサンプリングOFFでも書ける）

**理由**: `PROFILING_SAMPLE_INTERVAL` を使うと終了時に profiling が OFF のままになり、`timing_report.csv` が出ないケースがあったため。

#### 4) rc掃引用の配列ジョブを追加
- 追加: `job_aor_cutoff_sweep.sh`
- 変更点:
  - `inputs/input_coulomb_short.dat` をコピーし、`COULOMB_CUTOFF` を 0.02/0.03/0.04/0.05 に上書き
  - inputから `CONTAINER_WIDTH` を抽出し、AoR解析（`--width`）に渡す

**理由**: rc掃引を「同一条件・自動終了・自動解析」で回し、比較がブレないようにする。

#### 5) 集計・推奨rc選定・input反映スクリプトを追加
- 追加:
  - `src/summarize_cutoff_sweep.py`（AoR/時間内訳/ΔAoR を CSV に集約）
  - `src/apply_recommended_cutoff.py`（集計から推奨rcを選び、inputの `COULOMB_CUTOFF` を更新）
  - `postprocess_cutoff_sweep.sh`（上記をまとめて実行。既定は dry-run）

**理由**: rc決定を手作業にせず、**AoR差±0.5°で最小rcを機械的に選定**できるようにするため。

### 出力ファイル（チェックポイント）
- 各ケース（例: `results/aor_cutoff_sweep/rc_0.03/...`）
  - `particles.csv`（AoR解析の入力）
  - `repose_angle_results.csv`（AoR結果）
  - `timing_report.csv`（`coulomb_force` 等の内訳時間）
  - `timing.csv`（ジョブの壁時計）
- 集計:
  - `results/aor_cutoff_sweep/sweep_summary.csv`


