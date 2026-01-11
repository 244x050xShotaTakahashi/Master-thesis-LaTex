! メインプログラムおよびモジュール: 2次元粒子要素法 (球体モデル)

! モジュール: シミュレーション定数 (配列サイズ、数学定数)
module simulation_constants_mod
    implicit none
    integer, parameter :: ni_max = 50000  ! ni: 最大粒子数
    integer, parameter :: nj_max = 40    ! nj: 粒子ごとの最大接触点数 (粒子間10 + 壁(軸/斜面)30)
    integer, parameter :: nc_max = 500000 ! nc: グリッド内の最大セル数 
    real(8), parameter :: PI_VAL = 3.141592653589793d0 ! pi: 円周率
    real(8), parameter :: GRAVITY_ACCEL = 9.80665d0    ! g: 重力加速度
end module simulation_constants_mod

! モジュール: シミュレーション制御パラメータと材料物性値
module simulation_parameters_mod
    use simulation_constants_mod, only: ni_max
    implicit none
    
    ! シミュレーション制御パラメータ
    real(8) :: time_step                      ! dt: 時間刻み
    real(8) :: friction_coeff_particle        ! fri: 粒子間摩擦係数
    real(8) :: friction_coeff_wall            ! frw: 壁-粒子間摩擦係数
    real(8) :: rolling_friction_coeff_particle ! 粒子間転がり摩擦係数
    real(8) :: rolling_friction_coeff_wall     ! 壁-粒子間転がり摩擦係数
    real(8) :: young_modulus_particle         ! e: 粒子のヤング率
    real(8) :: young_modulus_wall             ! ew: 壁のヤング率
    real(8) :: poisson_ratio_particle         ! po: 粒子のポアソン比
    real(8) :: poisson_ratio_wall             ! pow: 壁のポアソン比
    real(8) :: shear_to_normal_stiffness_ratio ! so: せん断弾性係数と法線方向弾性係数の比
    real(8) :: particle_density               ! de: 粒子の密度
    real(8) :: reference_overlap            ! 参照食い込み量 δ_ref（平均半径の5%）
    logical :: stop_when_static             ! 静止検出時に計算を停止するか（1=停止, 0=続行）
    real(8) :: kinetic_energy_threshold     ! 運動エネルギー閾値 [J]（静止判定用）
    integer :: min_steps_before_static_check ! 静止判定を有効にするまでに必要な最小ステップ数
    
    ! 粒子生成パラメータ
    real(8) :: particle_radius_large          ! r1: 大きな粒子の半径
    real(8) :: particle_radius_small          ! r2: 小さな粒子の半径
    real(8) :: container_width                ! w: 容器の幅
    real(8) :: container_height               ! h: 容器の高さ (0.0=制限なし)
    integer :: particle_gen_layers            ! ipz: 初期粒子生成層数
    integer :: random_seed                    ! 乱数シード
    
    ! セル法アルゴリズム制御パラメータ
    logical :: disable_cell_algorithm         ! セル法アルゴリズムを無効化するフラグ
    real(8) :: cell_size_override            ! セルサイズの手動設定値 (0.0=自動計算)
    
    ! 出力制御パラメータ
    integer :: output_interval                ! 出力間隔 [ステップ]
    integer :: max_calculation_steps          ! 最大計算ステップ数
    
    ! 明示座標入力の制御
    logical :: use_explicit_positions         ! 明示座標ファイルの有無で切替
    character(len=256) :: positions_file      ! 明示座標ファイルパス
    character(len=256) :: explicit_positions_file = ''  ! 入力ファイルで指定する明示座標ファイル（指定時は最優先）
    
    ! クーロン力関連パラメータ
    real(8) :: coulomb_constant               ! k: クーロン定数 [N⋅m²/C²]
    real(8) :: default_charge                 ! デフォルトの粒子電荷 [C]
    logical :: enable_coulomb_force           ! クーロン力の有効化フラグ
    real(8) :: coulomb_cutoff = 0.1d0         ! クーロン力カットオフ半径 [m] (0=カットオフなし)
    real(8) :: coulomb_cutoff_radius_mult = 0.0d0  ! カットオフを平均半径の倍数で指定 (0=絶対値使用)
    real(8) :: coulomb_softening = 0.0d0      ! ソフトニングパラメータ δ [m] (発散回避用)
    logical :: coulomb_shift_force = .true.   ! シフトドフォースを使用するか（カットオフで力を連続に）
    logical :: coulomb_use_cell = .true.      ! セル法を使用するか（false=旧来の全対全計算）
    character(len=256) :: output_dir          ! 出力ディレクトリパス
    logical :: enable_profiling = .true.      ! プロファイル機能の有効化
    integer :: profiling_sample_interval = 0  ! サンプリング計測間隔 (0=サンプリングなし、N=Nステップごとに有効化)

    ! OpenMP 実行の再現性重視モード（1=OpenMP並列を避けて決定的に実行）
    logical :: deterministic_omp = .false.
    
    ! 半径分布設定
    character(len=32) :: radius_distribution_type = 'fixed'  ! 'fixed', 'file', 'uniform', 'normal'
    character(len=256) :: radius_list_file = ''              ! 半径リストファイルパス
    real(8) :: radius_uniform_min = 0.005d0                  ! 一様分布の最小値 [m]
    real(8) :: radius_uniform_max = 0.010d0                  ! 一様分布の最大値 [m]
    real(8) :: radius_normal_mean = 0.0075d0                 ! 正規分布の平均値 [m]
    real(8) :: radius_normal_std = 0.001d0                   ! 正規分布の標準偏差 [m]
    
    ! 電荷分布設定
    character(len=32) :: charge_distribution_type = 'fixed'  ! 'fixed', 'file', 'uniform', 'normal'
    character(len=256) :: charge_list_file = ''              ! 電荷リストファイルパス
    real(8) :: charge_uniform_min = 0.0d0                    ! 一様分布の最小値 [C]
    real(8) :: charge_uniform_max = 1.0d-9                   ! 一様分布の最大値 [C]
    real(8) :: charge_normal_mean = 1.0d-9                   ! 正規分布の平均値 [C]
    real(8) :: charge_normal_std = 0.2d-9                    ! 正規分布の標準偏差 [C]
    
    ! 指数分布パラメータ
    real(8) :: radius_exponential_mean = 0.0075d0            ! 指数分布の平均値 [m]
    real(8) :: charge_exponential_mean = 1.0d-9              ! 指数分布の平均値 [C]
    
    ! 二峰性分布パラメータ（半径）
    real(8) :: radius_bimodal_mean1 = 0.005d0                ! 第1ピークの平均値 [m]
    real(8) :: radius_bimodal_std1 = 0.001d0                 ! 第1ピークの標準偏差 [m]
    real(8) :: radius_bimodal_mean2 = 0.010d0                ! 第2ピークの平均値 [m]
    real(8) :: radius_bimodal_std2 = 0.001d0                 ! 第2ピークの標準偏差 [m]
    real(8) :: radius_bimodal_ratio = 0.5d0                  ! 第1ピークの割合 [-]
    
    ! 二峰性分布パラメータ（電荷）
    real(8) :: charge_bimodal_mean1 = 1.0d-9                 ! 第1ピークの平均値 [C]
    real(8) :: charge_bimodal_std1 = 0.2d-9                  ! 第1ピークの標準偏差 [C]
    real(8) :: charge_bimodal_mean2 = -1.0d-9                ! 第2ピークの平均値 [C]
    real(8) :: charge_bimodal_std2 = 0.2d-9                  ! 第2ピークの標準偏差 [C]
    real(8) :: charge_bimodal_ratio = 0.5d0                  ! 第1ピークの割合 [-]
    
    ! 壁引き抜きパラメータ
    logical :: enable_wall_withdraw = .false.                ! 壁引き抜きの有効化
    integer :: wall_withdraw_step = 0                        ! 引き抜き開始ステップ
    integer :: withdraw_wall_id = 0                          ! 引き抜く壁 (1=左,2=下,3=右,4=上, 0=無効)
    integer :: withdraw_sloped_wall_id = 0                   ! 引き抜く斜面壁ID (walls.datの壁番号, 0=無効)
    logical :: left_wall_active = .true.                     ! 左壁アクティブ状態
    logical :: bottom_wall_active = .true.                   ! 下壁アクティブ状態
    logical :: right_wall_active = .true.                    ! 右壁アクティブ状態
    logical :: top_wall_active = .true.                      ! 上壁アクティブ状態
    logical :: wall_withdraw_done = .false.                  ! 壁引き抜き実行済みフラグ
    logical :: filled_particles_saved = .false.              ! 充填状態保存済みフラグ
    
    ! 自動壁引き抜き制御パラメータ
    logical :: auto_wall_withdraw = .false.                  ! 自動壁引き抜きモードの有効化
    integer :: simulation_phase = 1                          ! シミュレーションフェーズ (1=充填, 2=崩落, 3=再堆積)
    logical :: energy_risen_after_withdraw = .false.         ! 壁引き抜き後に運動エネルギーが上昇したか
    real(8) :: peak_kinetic_energy = 0.0d0                   ! 崩落中の最大運動エネルギー（デバッグ用）
    
    save
end module simulation_parameters_mod

! モジュール: 粒子固有データ (物理特性、運動学、力)
module particle_data_mod
    use simulation_constants_mod, only: ni_max, nj_max
    implicit none

    ! 物理特性
    real(8), dimension(ni_max) :: radius         ! rr(ni): 粒子半径
    real(8), dimension(ni_max) :: mass           ! wei(ni): 粒子質量
    real(8), dimension(ni_max) :: moment_inertia ! pmi(ni): 粒子の慣性モーメント
    real(8), dimension(ni_max) :: charge         ! q(ni): 粒子の電荷 [C]

    ! 位置と向き
    real(8), dimension(ni_max) :: x_coord        ! x0(ni): 粒子中心のx座標
    real(8), dimension(ni_max) :: z_coord        ! z0(ni): 粒子中心のz座標 (原文ではy、コードではz)
    real(8), dimension(ni_max) :: rotation_angle ! qq(ni): 粒子の回転変位 (角度)

    ! 速度 (並進および回転)
    real(8), dimension(ni_max) :: x_vel          ! u0(ni): 粒子のx方向速度
    real(8), dimension(ni_max) :: z_vel          ! v0(ni): 粒子のz方向速度
    real(8), dimension(ni_max) :: rotation_vel   ! f0(ni): 粒子の回転速度

    ! 合力とモーメント
    real(8), dimension(ni_max) :: x_force_sum    ! xf(ni): 粒子に働くx方向の合力
    real(8), dimension(ni_max) :: z_force_sum    ! zf(ni): 粒子に働くz方向の合力
    real(8), dimension(ni_max) :: moment_sum     ! mf(ni): 粒子に働くモーメント

    ! 接触力の成分と接触相手のインデックス
    real(8), dimension(ni_max, nj_max) :: normal_force_contact  ! en(ni,nj): 法線方向接触力
    real(8), dimension(ni_max, nj_max) :: shear_force_contact   ! es(ni,nj): せん断方向接触力
    integer, dimension(ni_max, nj_max) :: contact_partner_idx ! je(ni,nj): 接触点番号配列 (接触している粒子/壁のインデックスを格納)
    real(8), dimension(ni_max, nj_max) :: previous_overlap      ! 前ステップのオーバーラップ（接触開始検出用）

    ! 現時間ステップにおける増分変位 (common/dpm/ より)
    real(8), dimension(ni_max) :: x_disp_incr    ! u(ni): x方向変位増分
    real(8), dimension(ni_max) :: z_disp_incr    ! v(ni): z方向変位増分
    real(8), dimension(ni_max) :: rot_disp_incr  ! f(ni): 回転変位増分
    
    save
end module particle_data_mod

! モジュール: セル格子システムデータ
module cell_system_mod
    use simulation_constants_mod, only: ni_max, nc_max
    implicit none

    integer :: num_particles          ! n: 粒子数
    integer :: cells_x_dir            ! idx: x方向のセル数
    integer :: cells_z_dir            ! idz: z方向のセル数 (使用状況から推測)

    real(8) :: cell_size              ! c: セルの幅/サイズ

    ! 各セルに属する粒子を連結リストで保持する（セル先頭 index → next → ...）
    integer, dimension(nc_max) :: cell_head        ! そのセルで最初に登録された粒子インデックス (空=0)
    integer, dimension(ni_max) :: particle_cell_next ! 次の粒子インデックス (0 ならリスト終端)

    ! 後方互換用に「最後に登録された粒子」を保持（デバッグ用途）
    integer, dimension(nc_max) :: cell_particle_map

    integer, dimension(ni_max) :: particle_cell_idx ! 粒子iが格納されているセル番号

    save
end module cell_system_mod

! モジュール: 斜面壁データ
module wall_data_mod
    implicit none
    integer, parameter :: nw_max = 128
    integer :: num_walls = 0
    real(8), dimension(nw_max) :: wall_x_start
    real(8), dimension(nw_max) :: wall_z_start
    real(8), dimension(nw_max) :: wall_x_end
    real(8), dimension(nw_max) :: wall_z_end
    real(8), dimension(nw_max) :: wall_length
    real(8), dimension(nw_max) :: wall_tangent_x
    real(8), dimension(nw_max) :: wall_tangent_z
    real(8), dimension(nw_max) :: wall_normal_x
    real(8), dimension(nw_max) :: wall_normal_z
    logical, dimension(nw_max) :: wall_active = .true.     ! 壁のアクティブ状態
    character(len=256) :: walls_file = 'inputs/walls.dat'
    logical :: walls_file_exists = .false.
    save
end module wall_data_mod

! モジュール: プロファイリングユーティリティ
module profiling_mod
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter :: max_profile_entries = 512
    integer, parameter :: profiler_tick_kind = int64

    type profile_entry_t
        character(len=64) :: name = ''
        integer(int64) :: total_ticks = 0_int64
        integer(int64) :: max_ticks = 0_int64
        integer :: call_count = 0
    end type profile_entry_t

    type(profile_entry_t), dimension(max_profile_entries) :: profile_entries
    integer :: profile_entry_count = 0
    integer :: profile_clock_rate = 1
    logical :: profiling_enabled = .false.

contains

    subroutine profiler_init(enable_flag)
        logical, intent(in), optional :: enable_flag
        integer :: i

        profiling_enabled = .true.
        if (present(enable_flag)) profiling_enabled = enable_flag

        call system_clock(count_rate=profile_clock_rate)
        if (profile_clock_rate <= 0) profile_clock_rate = 1

        profile_entry_count = 0
        do i = 1, max_profile_entries
            profile_entries(i)%name = ''
            profile_entries(i)%total_ticks = 0_int64
            profile_entries(i)%max_ticks = 0_int64
            profile_entries(i)%call_count = 0
        end do
    end subroutine profiler_init

    subroutine profiler_time_now(ticks_out)
        integer(int64), intent(out) :: ticks_out
        call system_clock(count=ticks_out)
    end subroutine profiler_time_now

    subroutine profiler_set_enabled(flag)
        logical, intent(in) :: flag
        integer :: i

        profiling_enabled = flag
        if (.not. profiling_enabled) then
            profile_entry_count = 0
            do i = 1, max_profile_entries
                profile_entries(i)%name = ''
                profile_entries(i)%total_ticks = 0_int64
                profile_entries(i)%max_ticks = 0_int64
                profile_entries(i)%call_count = 0
            end do
        end if
    end subroutine profiler_set_enabled

    subroutine profiler_touch(name)
        character(len=*), intent(in) :: name
        integer :: idx
        character(len=64) :: key
        integer :: i

        if (.not. profiling_enabled) return

        key = adjustl(name)
        key = trim(key)

        !$omp critical(prof_registry)
        idx = 0
        do i = 1, profile_entry_count
            if (trim(profile_entries(i)%name) == key) then
                idx = i
                exit
            end if
        end do

        if (idx == 0) then
            if (profile_entry_count < max_profile_entries) then
                profile_entry_count = profile_entry_count + 1
                profile_entries(profile_entry_count)%name = key
            end if
        end if
        !$omp end critical(prof_registry)
    end subroutine profiler_touch

    subroutine profiler_start(name, token)
        character(len=*), intent(in) :: name
        integer(int64), intent(out) :: token

        if (.not. profiling_enabled) then
            token = 0_int64
            return
        end if
        if (.false.) call profiler_touch(name)
        call profiler_time_now(token)
    end subroutine profiler_start

    subroutine profiler_stop(name, token)
        character(len=*), intent(in) :: name
        integer(int64), intent(in) :: token
        integer(int64) :: now_ticks, elapsed

        if (.not. profiling_enabled) return

        call profiler_time_now(now_ticks)
        elapsed = now_ticks - token
        if (elapsed < 0_int64) return
        call profiler_add(name, elapsed)
    end subroutine profiler_stop

    subroutine profiler_add(name, elapsed_ticks)
        character(len=*), intent(in) :: name
        integer(int64), intent(in) :: elapsed_ticks
        integer :: idx
        character(len=64) :: key
        integer :: i
        logical :: skip_update
        
        if (.not. profiling_enabled) return
        if (elapsed_ticks < 0_int64) return
        
        key = adjustl(name)
        key = trim(key)
        skip_update = .false.
        
        !$omp critical(prof_registry)
        idx = 0
        do i = 1, profile_entry_count
            if (trim(profile_entries(i)%name) == key) then
                idx = i
                exit
            end if
        end do
        
        if (idx == 0) then
            if (profile_entry_count < max_profile_entries) then
                profile_entry_count = profile_entry_count + 1
                idx = profile_entry_count
                profile_entries(idx)%name = key
            else
                skip_update = .true.
            end if
        end if
        
        if (.not. skip_update) then
            profile_entries(idx)%total_ticks = profile_entries(idx)%total_ticks + elapsed_ticks
            profile_entries(idx)%call_count = profile_entries(idx)%call_count + 1
            if (elapsed_ticks > profile_entries(idx)%max_ticks) then
                profile_entries(idx)%max_ticks = elapsed_ticks
            end if
        end if
        !$omp end critical(prof_registry)
        
        if (skip_update) return
    end subroutine profiler_add

    real(8) function profiler_ticks_to_seconds(ticks) result(seconds)
        integer(int64), intent(in) :: ticks
        seconds = real(ticks, 8) / real(profile_clock_rate, 8)
    end function profiler_ticks_to_seconds

    real(8) function profiler_total_seconds(name) result(total_sec)
        character(len=*), intent(in) :: name
        integer :: i
        character(len=64) :: key

        total_sec = 0.0d0
        if (.not. profiling_enabled) return

        key = adjustl(name)
        key = trim(key)

        do i = 1, profile_entry_count
            if (trim(profile_entries(i)%name) == key) then
                total_sec = profiler_ticks_to_seconds(profile_entries(i)%total_ticks)
                return
            end if
        end do
    end function profiler_total_seconds

    subroutine profiler_write_csv(file_path)
        character(len=*), intent(in) :: file_path
        integer :: unit_num, ios, i
        real(8) :: total_all, avg_time, percent_share
        real(8) :: entry_total, entry_max

        ! profiling_enabled は「計測のオン/オフ（サンプリング）」にも使われるため、
        ! 終了時点で false でも、既に蓄積した結果は書き出せるようにする。
        if (profile_entry_count <= 0) return

        total_all = 0.0d0
        do i = 1, profile_entry_count
            total_all = total_all + profiler_ticks_to_seconds(profile_entries(i)%total_ticks)
        end do
        if (total_all <= 0.0d0) total_all = 1.0d0

        open(newunit=unit_num, file=trim(file_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'プロファイルCSVを開けません: ', trim(file_path)
            return
        end if

        write(unit_num,'(A)') 'name,call_count,total_seconds,average_seconds,max_seconds,percent_total'
        do i = 1, profile_entry_count
            entry_total = profiler_ticks_to_seconds(profile_entries(i)%total_ticks)
            entry_max = profiler_ticks_to_seconds(profile_entries(i)%max_ticks)
            if (profile_entries(i)%call_count > 0) then
                avg_time = entry_total / real(profile_entries(i)%call_count, 8)
            else
                avg_time = 0.0d0
            end if
            percent_share = (entry_total / total_all) * 100.0d0
            write(unit_num,'(A,",",I0,",",ES14.6,",",ES14.6,",",ES14.6,",",ES12.4)') &
                trim(profile_entries(i)%name), profile_entries(i)%call_count, entry_total, avg_time, entry_max, percent_share
        end do

        close(unit_num)
    end subroutine profiler_write_csv

    logical function profiler_is_enabled() result(flag)
        flag = profiling_enabled
    end function profiler_is_enabled

end module profiling_mod

! メインプログラム
program two_dimensional_dem
    use profiling_mod
    use omp_lib
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    use wall_data_mod
    implicit none

    integer :: it_step, static_judge_flag          ! static_judge_flag: 静止判定フラグ
    integer :: i ! ループカウンタ用にiを宣言
    real(8) :: current_time                        ! 現在時刻
    real(8) :: rmax_particle_radius                ! fpositから返される最大粒子半径
    logical :: leapfrog_initialized                ! 蛙飛び法の初期化フラグ
    
    ! 計算時間計測用変数
    integer :: start_time, end_time, clock_rate
    real(8) :: elapsed_time
    
    ! 蛙飛び法用速度記録変数（位置更新時の速度v(t+Δt/2)）
    real(8) :: z_vel_at_position_update
    integer :: interval_clock_last, interval_clock_curr
    real(8) :: neighbor_time_last, neighbor_time_now
    real(8) :: interval_elapsed, interval_neighbor_elapsed, interval_avg_step
    integer(profiler_tick_kind) :: neighbor_block_token
    integer(profiler_tick_kind) :: coulomb_token, integrate_token, output_token

    ! 計算時間計測開始
    call profiler_init(.true.)
    call system_clock(start_time, clock_rate)
    if (clock_rate <= 0) clock_rate = 1
    interval_clock_last = start_time
    neighbor_time_last = 0.0d0
    
    ! inputファイルからパラメータを読み込み
    call read_input_file
    call profiler_set_enabled(enable_profiling)

    ! OpenMP の再現性重視モード:
    ! - 並列実行による加算順序揺らぎ・スケジューリング揺らぎを避けるため、スレッド数を 1 に固定する。
    ! - 併せて OpenMP の動的スレッド調整を無効化する。
    if (deterministic_omp) then
        call omp_set_dynamic(.false.)
        call omp_set_num_threads(1)
        write(*,*) '[OMP] deterministic_omp=1: OpenMPスレッド数を1に固定'
    end if
    
    ! 初期位置と初期条件の設定
    call fposit_sub(rmax_particle_radius)
    call inmat_sub
    call init_sub
    
    ! クーロン力ソフトニングのデフォルト値設定（未指定の場合は最大半径の5%）
    if (coulomb_softening <= 0.0d0 .and. enable_coulomb_force) then
        coulomb_softening = rmax_particle_radius * 0.05d0
        write(*,*) 'クーロン力ソフトニングを自動設定: ', coulomb_softening, ' [m]'
    end if
    
    ! クーロン力カットオフの倍数指定からの計算（平均半径を使用）
    if (coulomb_cutoff_radius_mult > 0.0d0 .and. enable_coulomb_force) then
        block
            real(8) :: mean_radius
            if (num_particles > 0) then
                mean_radius = sum(radius(1:num_particles)) / real(num_particles, 8)
            else
                mean_radius = (particle_radius_large + particle_radius_small) / 2.0d0
            end if
            coulomb_cutoff = mean_radius * coulomb_cutoff_radius_mult
            write(*,*) '平均半径: ', mean_radius, ' [m]'
            write(*,*) 'カットオフを平均半径の', coulomb_cutoff_radius_mult, '倍に設定: ', coulomb_cutoff, ' [m]'
        end block
    end if
    
    ! クーロン力パラメータの表示
    if (enable_coulomb_force) then
        write(*,*) '=== クーロン力設定 ==='
        write(*,*) 'カットオフ半径: ', coulomb_cutoff, ' [m]'
        write(*,*) 'ソフトニング δ: ', coulomb_softening, ' [m]'
        write(*,*) 'シフトドフォース: ', coulomb_shift_force
        write(*,*) 'セル法使用: ', coulomb_use_cell
        write(*,*) '======================'
        
        ! 初期状態で新旧実装の検証を実行（デバッグ用）
        ! call verify_coulomb_implementations()
    end if

    current_time = 0.0d0
    leapfrog_initialized = .false.

    ! 各ステップの繰り返し計算
    do it_step = 1, max_calculation_steps
        current_time = current_time + time_step
        
        ! サンプリング計測: N ステップごとにプロファイリングを有効化
        if (profiling_sample_interval > 0) then
            if (mod(it_step, profiling_sample_interval) == 1) then
                call profiler_set_enabled(.true.)
            else
                call profiler_set_enabled(.false.)
            end if
        end if
        
        ! 壁引き抜き処理（手動モードのみ）
        if (enable_wall_withdraw .and. .not. wall_withdraw_done .and. .not. auto_wall_withdraw) then
            if (it_step >= wall_withdraw_step) then
                ! コンテナ壁の引き抜き（withdraw_wall_id > 0 の場合）
                if (withdraw_wall_id > 0) then
                    select case (withdraw_wall_id)
                        case (1)
                            left_wall_active = .false.
                            write(*,*) '左壁を引き抜きました。ステップ: ', it_step
                        case (2)
                            bottom_wall_active = .false.
                            write(*,*) '下壁を引き抜きました。ステップ: ', it_step
                        case (3)
                            right_wall_active = .false.
                            write(*,*) '右壁を引き抜きました。ステップ: ', it_step
                        case (4)
                            top_wall_active = .false.
                            write(*,*) '上壁を引き抜きました。ステップ: ', it_step
                    end select
                end if
                ! 斜面壁の引き抜き（withdraw_sloped_wall_id > 0 の場合）
                if (withdraw_sloped_wall_id > 0 .and. withdraw_sloped_wall_id <= num_walls) then
                    wall_active(withdraw_sloped_wall_id) = .false.
                    write(*,*) '斜面壁', withdraw_sloped_wall_id, 'を引き抜きました。ステップ: ', it_step
                end if
                wall_withdraw_done = .true.
            end if
        end if
        
        if (.not. leapfrog_initialized) then
            ! 初回のみ: v(0) → v(Δt/2) への変換
            call profiler_start('neighbor_search_contact', neighbor_block_token)
            call ncel_sub
            
            ! 全粒子の合力をクリア
            !$omp parallel do private(i)
            do i = 1, num_particles
                x_force_sum(i) = 0.0d0
                z_force_sum(i) = 0.0d0
                moment_sum(i) = 0.0d0
            end do
            !$omp end parallel do
            
            if (deterministic_omp) then
                do i = 1, num_particles
                    call wcont_sub(i)
                    call pcont_sub(i, rmax_particle_radius)
                end do
            else
                !$omp parallel do schedule(dynamic, 16) private(i)
                do i = 1, num_particles
                    ! 粒子と壁との接触力計算
                    call wcont_sub(i)
                    ! 粒子間の接触力計算
                    call pcont_sub(i, rmax_particle_radius)
                end do
                !$omp end parallel do
            end if

            call profiler_stop('neighbor_search_contact', neighbor_block_token)
            
            ! クーロン力の計算（セル法 or 全対全を選択）
            call profiler_start('coulomb_force', coulomb_token)
            if (coulomb_use_cell) then
                call coulomb_force_cell_cutoff_sub()
            else
                call coulomb_force_sub_full_pairs()
            end if
            call profiler_stop('coulomb_force', coulomb_token)
            
            call profiler_start('integrate_leapfrog', integrate_token)
            call nposit_leapfrog_sub(static_judge_flag, 0)
            call profiler_stop('integrate_leapfrog', integrate_token)
            leapfrog_initialized = .true.
        else
            ! 通常ループ: 蛙飛び法のメインステップ
            ! フェーズ1: 位置更新
            call profiler_start('integrate_leapfrog', integrate_token)
            call nposit_leapfrog_sub(static_judge_flag, 1)
            call profiler_stop('integrate_leapfrog', integrate_token)
            
            ! 新しい位置で力を計算
            call profiler_start('neighbor_search_contact', neighbor_block_token)
            call ncel_sub
            
            ! 全粒子の合力をクリア
            !$omp parallel do private(i)
            do i = 1, num_particles
                x_force_sum(i) = 0.0d0
                z_force_sum(i) = 0.0d0
                moment_sum(i) = 0.0d0
            end do
            !$omp end parallel do
            
            if (deterministic_omp) then
                do i = 1, num_particles
                    call wcont_sub(i)
                    call pcont_sub(i, rmax_particle_radius)
                end do
            else
                !$omp parallel do schedule(dynamic, 16) private(i)
                do i = 1, num_particles
                    ! 粒子と壁との接触力計算
                    call wcont_sub(i)
                    ! 粒子間の接触力計算
                    call pcont_sub(i, rmax_particle_radius)
                end do
                !$omp end parallel do
            end if

            call profiler_stop('neighbor_search_contact', neighbor_block_token)
            
            ! クーロン力の計算（セル法 or 全対全を選択）
            call profiler_start('coulomb_force', coulomb_token)
            if (coulomb_use_cell) then
                call coulomb_force_cell_cutoff_sub()
            else
                call coulomb_force_sub_full_pairs()
            end if
            call profiler_stop('coulomb_force', coulomb_token)
            
            ! フェーズ2: 速度更新（新しい位置での力を使用）
            call profiler_start('integrate_leapfrog', integrate_token)
            call nposit_leapfrog_sub(static_judge_flag, 2)
            call profiler_stop('integrate_leapfrog', integrate_token)
        end if

        ! ========== フェーズ制御による静止判定 ==========
        if (auto_wall_withdraw .and. enable_wall_withdraw) then
            ! 自動壁引き抜きモード
            select case (simulation_phase)
            case (1)  ! 充填段階
                if (static_judge_flag == 1 .and. it_step >= min_steps_before_static_check) then
                    ! 充填完了 → 壁引き抜き
                    if (.not. filled_particles_saved) then
                        call save_filled_particles_sub
                        filled_particles_saved = .true.
                    end if
                    
                    ! 壁引き抜き実行
                    if (withdraw_wall_id > 0) then
                        select case (withdraw_wall_id)
                            case (1); left_wall_active = .false.
                            case (2); bottom_wall_active = .false.
                            case (3); right_wall_active = .false.
                            case (4); top_wall_active = .false.
                        end select
                        write(*,*) '自動壁引き抜き: 壁ID=', withdraw_wall_id, ', ステップ=', it_step
                    end if
                    if (withdraw_sloped_wall_id > 0 .and. withdraw_sloped_wall_id <= num_walls) then
                        wall_active(withdraw_sloped_wall_id) = .false.
                        write(*,*) '自動壁引き抜き: 斜面壁ID=', withdraw_sloped_wall_id, ', ステップ=', it_step
                    end if
                    wall_withdraw_done = .true.
                    simulation_phase = 2  ! 崩落段階へ遷移
                    write(*,*) 'フェーズ遷移: 充填 → 崩落'
                end if
                
            case (2)  ! 崩落段階（運動エネルギー上昇を待つ）
                if (static_judge_flag == 0) then  ! 運動エネルギー > 閾値
                    energy_risen_after_withdraw = .true.
                    simulation_phase = 3  ! 再堆積段階へ遷移
                    write(*,*) 'フェーズ遷移: 崩落 → 再堆積（崩落開始を確認）'
                end if
                
            case (3)  ! 再堆積段階（静止検出で終了）
                if (static_judge_flag == 1) then
                    write(*,*) '安息角計測完了。静止状態に到達しました。時刻: ', current_time
                    write(*,*) '総ステップ数: ', it_step
                    goto 200
                end if
            end select
        else
            ! 従来モード（手動壁引き抜き）
            if (static_judge_flag == 1) then
                if (it_step >= min_steps_before_static_check) then
                    if (enable_wall_withdraw .and. .not. wall_withdraw_done .and. .not. filled_particles_saved) then
                        if (it_step < wall_withdraw_step) then
                            call save_filled_particles_sub
                            filled_particles_saved = .true.
                        end if
                    end if
                    
                    if (stop_when_static) then
                        write(*,*) '静止状態に到達しました。時刻: ', current_time
                        goto 200
                    end if
                end if
            end if
        end if

        ! 計算状況の出力
        if (mod(it_step, 10000) == 0) then
            write(*, '(A,F10.6,A,F12.6,A,F12.6)') 'Time= ', current_time, &
                                                 ' Z0(N)= ', z_coord(num_particles), &
                                                 ' V0(N)= ', z_vel(num_particles)
        end if

        ! グラフィック用データの出力
        if (it_step == 1 .or. mod(it_step, output_interval) == 0) then
            call profiler_start('output', output_token)
            call gfout_sub(it_step, current_time, rmax_particle_radius)
            call profiler_stop('output', output_token)
        end if

        if (enable_profiling .and. mod(it_step, 50000) == 0) then
            call system_clock(interval_clock_curr)
            interval_elapsed = real(interval_clock_curr - interval_clock_last, 8) / real(clock_rate, 8)
            if (interval_elapsed < 0.0d0) interval_elapsed = 0.0d0
            interval_avg_step = interval_elapsed / 50000.0d0
            interval_clock_last = interval_clock_curr
            neighbor_time_now = profiler_total_seconds('neighbor_search_contact')
            interval_neighbor_elapsed = neighbor_time_now - neighbor_time_last
            neighbor_time_last = neighbor_time_now
            write(*,'(A,I0,A,F12.4,A,F12.6,A,F12.6)') &
                '[Profile] 50000 step@', it_step, ' : Δt= ', interval_elapsed, ' s, Δt/step= ', interval_avg_step, &
                ' s, neighbor= ', interval_neighbor_elapsed
        end if
    end do

200 continue ! シミュレーションループ脱出用のラベル

    ! AoR解析（particles.csv の最終stepを使用）向けに、終了時の最終状態を必ず出力する
    call gfout_sub(it_step, current_time, rmax_particle_radius)

    ! バックアップデータの出力
    call bfout_sub
    call profiler_write_csv(trim(output_dir)//'/timing_report.csv')

    close(10) ! particles.csv
    close(11) ! contacts.csv
    close(13) ! backl.d

    ! 計算時間計測終了
    call system_clock(end_time)
    elapsed_time = real(end_time - start_time) / real(clock_rate)
    
    write(*,*) '================================='
    write(*,*) 'シミュレーション実行結果'
    write(*,*) '================================='
    write(*,*) '粒子数: ', num_particles
    write(*,*) '計算ステップ数: ', it_step
    write(*,*) '実行時間: ', elapsed_time, ' 秒'
    write(*,*) '1ステップあたりの平均時間: ', elapsed_time / real(it_step), ' 秒'
    write(*,*) 'コンテナ幅: ', container_width
    if (container_height > 0.0d0) then
        write(*,*) 'コンテナ高さ: ', container_height, ' (上壁あり)'
    else
        write(*,*) 'コンテナ高さ: 制限なし (上壁なし)'
    end if
    
    if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
        write(*,*) 'セル法アルゴリズム: 無効化'
        write(*,*) 'セルサイズ: ', cell_size
    else
        write(*,*) 'セル法アルゴリズム: 有効'
        write(*,*) 'セルサイズ: ', cell_size
        write(*,*) 'セル数 (X方向): ', cells_x_dir
        write(*,*) 'セル数 (Z方向): ', cells_z_dir
        write(*,*) '総セル数: ', cells_x_dir * cells_z_dir
    end if

    write(*,*) '================================='
    
    stop
contains

    !===============================================================
    ! 粒子の直径リスト (m) をファイルから読み込み，半径配列 r_part(:) を設定する
    !   - filename : 入力ファイル名
    !   - r_part   : 半径 [m] を格納する配列（呼び出し側で allocate 済み）
    !   - Nmax     : r_part の最大要素数 (= size(r_part))
    !   - Nread    : 実際に読み込んだ粒子数（出力）
    !
    ! 入力ファイル仕様：
    !   - テキストファイル
    !   - 1 列目に直径 [m]（実数）
    !   - ヘッダ行や # で始まるコメント行、空行はスキップ
    !===============================================================
    subroutine load_particle_radii_from_file(filename, r_part, Nmax, Nread)
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        character(len=*), intent(in)  :: filename
        real(8),          intent(out) :: r_part(:)
        integer,          intent(in)  :: Nmax
        integer,          intent(out) :: Nread

        integer :: unit
        integer :: ios_line, ios_val
        character(len=256) :: line
        real(8) :: d_tmp   ! 読み込んだ直径 [m]
        integer(profiler_tick_kind) :: prof_token

        call profiler_start('load_particle_radii_from_file', prof_token)
        Nread = 0

        ! ファイルを開く
        open(newunit=unit, file=filename, status='old', action='read', &
            iostat=ios_line)
        if (ios_line /= 0) then
            write(*,*) '*** ERROR[load_particle_radii_from_file]: ', &
                        'cannot open file: ', trim(filename)
            call profiler_stop('load_particle_radii_from_file', prof_token)
            stop
        end if

        do
            ! 1 行ずつ文字列として読み込む
            read(unit, '(A)', iostat=ios_line) line
            if (ios_line /= 0) exit   ! EOF or read error

            ! 空行はスキップ
            if (len_trim(line) == 0) cycle

            ! コメント行（先頭が #）はスキップ
            if (line(1:1) == '#') cycle

            ! 文字列から実数（直径）を読み取る
            read(line, *, iostat=ios_val) d_tmp
            if (ios_val /= 0) then
                ! 数値に変換できない行（ヘッダなど）はスキップ
                cycle
            end if

            ! 粒子数カウント
            Nread = Nread + 1
            if (Nread > Nmax) then
                write(*,*) '*** ERROR[load_particle_radii_from_file]: ', &
                            'too many particles in file: ', trim(filename)
                write(*,*) '    Nmax = ', Nmax
                call profiler_stop('load_particle_radii_from_file', prof_token)
                stop
            end if

            ! 直径 -> 半径に変換して格納
            r_part(Nread) = 0.5d0 * d_tmp
        end do

        close(unit)

        if (Nread == 0) then
            write(*,*) '*** ERROR[load_particle_radii_from_file]: ', &
                        'no valid data found in file: ', trim(filename)
            call profiler_stop('load_particle_radii_from_file', prof_token)
            stop
        end if

        call profiler_stop('load_particle_radii_from_file', prof_token)
    end subroutine load_particle_radii_from_file

    !===============================================================
    ! 粒子の電荷リスト [C] をファイルから読み込み，配列 q_part(:) を設定する
    !   - filename : 入力ファイル名
    !   - q_part   : 電荷 [C] を格納する配列
    !   - Nmax     : q_part の最大要素数 (= size(q_part))
    !   - Nread    : 実際に読み込んだ粒子数（出力）
    !
    ! 入力ファイル仕様：
    !   - テキストファイル
    !   - 1 列目に電荷 [C]（実数）
    !   - ヘッダ行や # で始まるコメント行、空行はスキップ
    !===============================================================
    subroutine load_particle_charges_from_file(filename, q_part, Nmax, Nread)
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        character(len=*), intent(in)  :: filename
        real(8),          intent(out) :: q_part(:)
        integer,          intent(in)  :: Nmax
        integer,          intent(out) :: Nread

        integer :: unit
        integer :: ios_line, ios_val
        character(len=256) :: line
        real(8) :: q_tmp   ! 読み込んだ電荷 [C]
        integer(profiler_tick_kind) :: prof_token

        call profiler_start('load_particle_charges_from_file', prof_token)
        Nread = 0

        ! ファイルを開く
        open(newunit=unit, file=filename, status='old', action='read', &
            iostat=ios_line)
        if (ios_line /= 0) then
            write(*,*) '*** ERROR[load_particle_charges_from_file]: ', &
                        'cannot open file: ', trim(filename)
            call profiler_stop('load_particle_charges_from_file', prof_token)
            stop
        end if

        do
            ! 1 行ずつ文字列として読み込む
            read(unit, '(A)', iostat=ios_line) line
            if (ios_line /= 0) exit   ! EOF or read error

            ! 空行はスキップ
            if (len_trim(line) == 0) cycle

            ! コメント行（先頭が #）はスキップ
            if (line(1:1) == '#') cycle

            ! 文字列から実数（電荷）を読み取る
            read(line, *, iostat=ios_val) q_tmp
            if (ios_val /= 0) then
                ! 数値に変換できない行（ヘッダなど）はスキップ
                cycle
            end if

            ! 粒子数カウント
            Nread = Nread + 1
            if (Nread > Nmax) then
                write(*,*) '*** ERROR[load_particle_charges_from_file]: ', &
                            'too many particles in file: ', trim(filename)
                write(*,*) '    Nmax = ', Nmax
                call profiler_stop('load_particle_charges_from_file', prof_token)
                stop
            end if

            ! 電荷を格納
            q_part(Nread) = q_tmp
        end do

        close(unit)

        if (Nread == 0) then
            write(*,*) '*** ERROR[load_particle_charges_from_file]: ', &
                        'no valid data found in file: ', trim(filename)
            call profiler_stop('load_particle_charges_from_file', prof_token)
            stop
        end if

        call profiler_stop('load_particle_charges_from_file', prof_token)
    end subroutine load_particle_charges_from_file

    !> inputファイルからパラメータを読み込むサブルーチン
    subroutine read_input_file
        use wall_data_mod
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        character(len=256) :: line, keyword
        character(len=256) :: input_filename
        integer :: ios, unit_num
        real(8) :: value
        integer(profiler_tick_kind) :: prof_token
        
        call profiler_start('read_input_file', prof_token)
        ! inputファイル名の決定（コマンドライン引数または固定名）
        if (command_argument_count() > 0) then
            call get_command_argument(1, input_filename)
            ! 相対パスの場合、inputsフォルダを追加
            if (input_filename(1:1) /= '/' .and. input_filename(1:6) /= 'inputs') then
                input_filename = 'inputs/' // trim(input_filename)
            end if

            ! 第2引数がある場合は出力ディレクトリとして取得
            if (command_argument_count() >= 2) then
                call get_command_argument(2, output_dir)
            else
                output_dir = "data"
            end if
        else
            input_filename = "inputs/input_valid.dat"
            output_dir = "data"
        end if
        
        write(*,*) 'inputファイルを読み込み中: ', trim(input_filename)
        write(*,*) '出力ディレクトリ: ', trim(output_dir)
        
        unit_num = 20
        open(unit=unit_num, file=input_filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: inputファイルを開けません: ', trim(input_filename)
            call profiler_stop('read_input_file', prof_token)
            stop
        end if
        
        ! デフォルト値の設定
        time_step = 5.0d-7
        friction_coeff_particle = 0.25d0
        friction_coeff_wall = 0.17d0
        rolling_friction_coeff_particle = 0.0d0
        rolling_friction_coeff_wall = 0.0d0
        young_modulus_particle = 4.9d9
        young_modulus_wall = 3.9d9
        poisson_ratio_particle = 0.23d0
        poisson_ratio_wall = 0.25d0
        particle_density = 2.48d3
        particle_radius_large = 1.0d-2
        particle_radius_small = 5.0d-3
        container_width = 5.0d-1
        container_height = 0.0d0  ! 0.0=制限なし（上壁なし）
        particle_gen_layers = 30
        random_seed = 584287
        disable_cell_algorithm = .false.
        stop_when_static = .true.
        kinetic_energy_threshold = 1.0d-6
        min_steps_before_static_check = 100000
        cell_size_override = 0.0d0
        output_interval = 50000
        max_calculation_steps = 2000000
        enable_coulomb_force = .false.
        coulomb_constant = 8.99d9  ! クーロン定数 k = 1/(4πε₀) [N⋅m²/C²]
        default_charge = 0.0d0      ! デフォルト電荷 [C]
        coulomb_cutoff = 0.1d0      ! カットオフ半径 [m]（0でカットオフなし）
        coulomb_cutoff_radius_mult = 0.0d0  ! カットオフを半径の倍数で指定（0=絶対値使用）
        coulomb_softening = 0.0d0   ! ソフトニング δ [m]（後で平均半径の5%をデフォルトに設定）
        coulomb_shift_force = .true. ! シフトドフォース（カットオフで力連続）
        coulomb_use_cell = .true.   ! セル法使用（false=旧来の全対全）
        enable_profiling = .true.
        deterministic_omp = .false.
        explicit_positions_file = ''
        
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! コメント行と空行をスキップ
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            
            ! キーワードを読み取り
            read(line, *, iostat=ios) keyword
            if (ios /= 0) cycle
            
            select case (trim(keyword))
                ! 数値パラメータ
                case ('TIME_STEP')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) time_step = value
                case ('FRICTION_COEFF_PARTICLE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) friction_coeff_particle = value
                case ('FRICTION_COEFF_WALL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) friction_coeff_wall = value
                case ('ROLLING_FRICTION_COEFF_PARTICLE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) rolling_friction_coeff_particle = value
                case ('ROLLING_FRICTION_COEFF_WALL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) rolling_friction_coeff_wall = value
                case ('YOUNG_MODULUS_PARTICLE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) young_modulus_particle = value
                case ('YOUNG_MODULUS_WALL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) young_modulus_wall = value
                case ('POISSON_RATIO_PARTICLE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) poisson_ratio_particle = value
                case ('POISSON_RATIO_WALL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) poisson_ratio_wall = value
                case ('PARTICLE_DENSITY')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) particle_density = value
                case ('PARTICLE_RADIUS_LARGE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) particle_radius_large = value
                case ('PARTICLE_RADIUS_SMALL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) particle_radius_small = value
                case ('CONTAINER_WIDTH')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) container_width = value
                case ('CONTAINER_HEIGHT')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) container_height = value
                case ('PARTICLE_GEN_LAYERS')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) particle_gen_layers = int(value)
                case ('RANDOM_SEED')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) random_seed = int(value)
                case ('DISABLE_CELL_ALGORITHM')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) disable_cell_algorithm = (int(value) == 1)
                case ('STOP_WHEN_STATIC')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) stop_when_static = (int(value) == 1)
                case ('KINETIC_ENERGY_THRESHOLD')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) kinetic_energy_threshold = value
                case ('CELL_SIZE_OVERRIDE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) cell_size_override = value
                case ('MAX_CALCULATION_STEPS')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) max_calculation_steps = int(value)
                case ('MIN_STEPS_BEFORE_STATIC_CHECK')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) min_steps_before_static_check = int(value)
                case ('ENABLE_COULOMB_FORCE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) enable_coulomb_force = (int(value) == 1)
                case ('COULOMB_CONSTANT')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_constant = value
                case ('DEFAULT_CHARGE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) default_charge = value
                case ('COULOMB_CUTOFF')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_cutoff = value
                case ('COULOMB_CUTOFF_RADIUS_MULT')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_cutoff_radius_mult = value
                case ('COULOMB_SOFTENING')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_softening = value
                case ('COULOMB_SHIFT_FORCE')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_shift_force = (int(value) == 1)
                case ('COULOMB_USE_CELL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) coulomb_use_cell = (int(value) == 1)
                case ('OUTPUT_INTERVAL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) output_interval = int(value)
                case ('ENABLE_PROFILING')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) enable_profiling = (int(value) == 1)
                case ('PROFILING_SAMPLE_INTERVAL')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) profiling_sample_interval = int(value)
                case ('DETERMINISTIC_OMP')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) deterministic_omp = (int(value) == 1)

                ! 明示座標ファイル（直接参照用）
                case ('EXPLICIT_POSITIONS_FILE')
                    read(line, *, iostat=ios) keyword, explicit_positions_file
                    if (ios /= 0) then
                        write(*,*) '警告: EXPLICIT_POSITIONS_FILE の読み込みに失敗: ', trim(line)
                        explicit_positions_file = ''
                    end if
                
                ! 半径分布パラメータ
                case ('RADIUS_DISTRIBUTION_TYPE')
                    read(line, *, iostat=ios) keyword, radius_distribution_type
                case ('RADIUS_LIST_FILE')
                    read(line, *, iostat=ios) keyword, radius_list_file
                case ('RADIUS_UNIFORM_MIN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_uniform_min = value
                case ('RADIUS_UNIFORM_MAX')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_uniform_max = value
                case ('RADIUS_NORMAL_MEAN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_normal_mean = value
                case ('RADIUS_NORMAL_STD')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_normal_std = value
                case ('RADIUS_EXPONENTIAL_MEAN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_exponential_mean = value
                case ('RADIUS_BIMODAL_MEAN1')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_bimodal_mean1 = value
                case ('RADIUS_BIMODAL_STD1')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_bimodal_std1 = value
                case ('RADIUS_BIMODAL_MEAN2')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_bimodal_mean2 = value
                case ('RADIUS_BIMODAL_STD2')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_bimodal_std2 = value
                case ('RADIUS_BIMODAL_RATIO')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) radius_bimodal_ratio = value
                
                ! 電荷分布パラメータ
                case ('CHARGE_DISTRIBUTION_TYPE')
                    read(line, *, iostat=ios) keyword, charge_distribution_type
                case ('CHARGE_LIST_FILE')
                    read(line, *, iostat=ios) keyword, charge_list_file
                case ('CHARGE_UNIFORM_MIN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_uniform_min = value
                case ('CHARGE_UNIFORM_MAX')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_uniform_max = value
                case ('CHARGE_NORMAL_MEAN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_normal_mean = value
                case ('CHARGE_NORMAL_STD')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_normal_std = value
                case ('CHARGE_EXPONENTIAL_MEAN')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_exponential_mean = value
                case ('CHARGE_BIMODAL_MEAN1')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_bimodal_mean1 = value
                case ('CHARGE_BIMODAL_STD1')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_bimodal_std1 = value
                case ('CHARGE_BIMODAL_MEAN2')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_bimodal_mean2 = value
                case ('CHARGE_BIMODAL_STD2')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_bimodal_std2 = value
                case ('CHARGE_BIMODAL_RATIO')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) charge_bimodal_ratio = value
                
                ! 壁引き抜きパラメータ
                case ('ENABLE_WALL_WITHDRAW')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) enable_wall_withdraw = (int(value) == 1)
                case ('WALL_WITHDRAW_STEP')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) wall_withdraw_step = int(value)
                case ('WITHDRAW_WALL_ID')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) withdraw_wall_id = int(value)
                case ('WITHDRAW_SLOPED_WALL_ID')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) withdraw_sloped_wall_id = int(value)
                case ('AUTO_WALL_WITHDRAW')
                    read(line, *, iostat=ios) keyword, value
                    if (ios == 0) auto_wall_withdraw = (int(value) == 1)
                
                case default
                    write(*,*) '警告: 不明なキーワード: ', trim(keyword)
            end select
        end do
        
        close(unit_num)
        
        write(*,*) 'inputファイルの読み込み完了'
        
        ! 数値積分法の表示
        write(*,*) '数値積分法: 蛙飛び法'
        
        ! 明示座標ファイルの存在チェック
        ! 優先順位:
        !   1) EXPLICIT_POSITIONS_FILE（入力で指定されていれば最優先）
        !   2) inputs/filled_particles.dat
        !   3) inputs/particle_positions.dat
        use_explicit_positions = .false.
        
        if (len_trim(explicit_positions_file) > 0) then
            positions_file = trim(explicit_positions_file)
            inquire(file=trim(positions_file), exist=use_explicit_positions)
            if (.not. use_explicit_positions) then
                write(*,*) 'エラー: 指定された明示座標ファイルが見つかりません: ', trim(positions_file)
                stop 'read_input_file: EXPLICIT_POSITIONS_FILE not found'
            end if
            write(*,*) '粒子配置: 明示座標ファイルを使用(指定): ', trim(positions_file)
        else
        positions_file = 'inputs/filled_particles.dat'
        inquire(file=trim(positions_file), exist=use_explicit_positions)
        if (use_explicit_positions) then
            write(*,*) '粒子配置: 充填状態ファイルを使用: ', trim(positions_file)
        else
            positions_file = 'inputs/particle_positions.dat'
            inquire(file=trim(positions_file), exist=use_explicit_positions)
            if (use_explicit_positions) then
                write(*,*) '粒子配置: 明示座標ファイルを使用: ', trim(positions_file)
            else
                write(*,*) '粒子配置: 乱数生成（明示座標ファイルなし）'
                end if
            end if
        end if

        call read_walls_file
        call profiler_stop('read_input_file', prof_token)
    end subroutine read_input_file

    !> 斜面壁ファイルを読み込むサブルーチン
    subroutine read_walls_file
        use wall_data_mod
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        integer :: unit_num, ios, line_no
        character(len=256) :: line
        real(8) :: x1, z1, x2, z2
        real(8) :: dx, dz, length_val
        logical :: file_exists
        integer(profiler_tick_kind) :: prof_token
        call profiler_start('read_walls_file', prof_token)
        
        walls_file_exists = .false.
        num_walls = 0
        
        inquire(file=trim(walls_file), exist=file_exists)
        if (.not. file_exists) then
            write(*,*) '斜面壁ファイルなし: ', trim(walls_file)
            call profiler_stop('read_walls_file', prof_token)
            return
        end if
        
        unit_num = 22
        open(unit=unit_num, file=trim(walls_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: 斜面壁ファイルを開けません: ', trim(walls_file)
            call profiler_stop('read_walls_file', prof_token)
            stop 'read_walls_file: open failed'
        end if
        
        line_no = 0
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_no = line_no + 1
            
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            
            read(line, *, iostat=ios) x1, z1, x2, z2
            if (ios /= 0) then
                write(*,*) '警告: 斜面壁ファイル解析失敗 (行', line_no, '): ', trim(line)
                cycle
            end if
            
            dx = x2 - x1
            dz = z2 - z1
            length_val = sqrt(dx*dx + dz*dz)
            if (length_val < 1.0d-10) then
                write(*,*) '警告: 壁長さが短すぎるためスキップ (行', line_no, ')'
                cycle
            end if
            
            if (num_walls >= nw_max) then
                write(*,*) 'エラー: 壁数が上限 (', nw_max, ') を超過しました'
                stop 'read_walls_file: too many walls'
            end if
            
            num_walls = num_walls + 1
            wall_x_start(num_walls) = x1
            wall_z_start(num_walls) = z1
            wall_x_end(num_walls) = x2
            wall_z_end(num_walls) = z2
            wall_length(num_walls) = length_val
            wall_tangent_x(num_walls) = dx / length_val
            wall_tangent_z(num_walls) = dz / length_val
            wall_normal_x(num_walls) = -wall_tangent_z(num_walls)
            wall_normal_z(num_walls) =  wall_tangent_x(num_walls)
            wall_active(num_walls) = .true.
        end do
        
        close(unit_num)
        
        if (num_walls > 0) then
            walls_file_exists = .true.
            write(*,*) '斜面壁を読み込み: ', num_walls, ' 本 (', trim(walls_file), ')'
        else
            write(*,*) '斜面壁ファイルに有効な壁がありません: ', trim(walls_file)
        end if
        call profiler_stop('read_walls_file', prof_token)
    end subroutine read_walls_file

    !> 初期粒子配置と構成を設定するサブルーチン
    subroutine fposit_sub(rmax_out)
        
        use simulation_constants_mod, only: ni_max, PI_VAL
        use simulation_parameters_mod
        use particle_data_mod, only: radius, x_coord, z_coord, rotation_angle, charge
        use cell_system_mod ! モジュールからセルシステム関連変数を取得
        use wall_data_mod, only: num_walls, wall_x_start
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none

        real(8), intent(out) :: rmax_out ! 出力: 最大粒子半径

        integer :: i_layer, j_particle_in_layer, ipx_calc, current_particle_count
        real(8) :: r1_val, r2_val, rn_val, dx_offset, random_uniform_val
        real(8) :: rmin_val ! 宣言をここに移動
        integer :: particles_this_row ! 宣言をここに移動
        ! 明示座標入力用
        integer :: ios
        character(len=256) :: line
        real(8) :: xin, zin, rin, qin
        integer :: read_count
        ! z座標範囲制限用
        real(8) :: z_min_target, z_max_target, z_range, layer_height
        integer :: max_layers_for_range, original_layers
        integer(profiler_tick_kind) :: prof_token
        ! 壁引き抜き用粒子配置範囲
        real(8) :: effective_container_width
        
        ! 分布用変数
        real(8), allocatable :: radius_list(:), charge_list(:)
        integer :: radius_list_count, charge_list_count
        real(8) :: generated_radius, generated_charge
        integer :: radius_list_idx, charge_list_idx
        
        ! 衝突チェック付き配置用変数
        real(8) :: r_avg                          ! 平均半径
        real(8) :: x_candidate, z_candidate       ! 候補位置
        logical :: placement_success              ! 配置成功フラグ
        integer :: retry_count, max_retry         ! リトライカウンタ
        integer :: collision_skip_count           ! 衝突スキップカウント
        logical :: has_collision                  ! 衝突フラグ
        integer :: seed_charge                    ! 電荷生成用の乱数シード（配置用と分離）

        call profiler_start('fposit_sub', prof_token)
        r1_val = particle_radius_large
        r2_val = particle_radius_small

        ! 電荷生成の乱数系列を粒子配置（衝突回避リトライ等）から分離して再現性を上げる。
        ! 同一 RANDOM_SEED のもとで決定的に生成しつつ、配置ロジック変更で電荷が変わりにくくする。
        seed_charge = random_seed + 104729
        if (seed_charge == 0) seed_charge = 1
        
        if (use_explicit_positions) then
            ! 明示座標ファイルから読み込み
            num_particles = 0
            rmax_out = 0.0d0
            rmin_val = 1.0d99
            open(unit=21, file=trim(positions_file), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) 'エラー: 明示座標ファイルを開けません: ', trim(positions_file)
                stop 'fposit_sub: 位置ファイルopen失敗'
            end if
            do
                read(21,'(A)', iostat=ios) line
                if (ios /= 0) exit
                if (len_trim(line) == 0) cycle
                if (line(1:1) == '#' .or. line(1:1) == '!') cycle
                ! 電荷を含めて読み込み (4列目がなければデフォルト値を使用)
                qin = default_charge
                read(line, *, iostat=ios) xin, zin, rin, qin
                if (ios /= 0) then
                    ! 4列目がない場合は3列のみ読み込み
                    read(line, *, iostat=ios) xin, zin, rin
                    if (ios /= 0) cycle
                    qin = default_charge
                end if
                if (rin <= 0.0d0) cycle
                num_particles = num_particles + 1
                if (num_particles > ni_max) then
                    write(*,*) 'エラー: 粒子数がni_maxを超過: ', ni_max
                    stop 'fposit_sub: 粒子が多すぎます'
                end if
                x_coord(num_particles) = xin
                z_coord(num_particles) = zin
                radius(num_particles)  = rin
                charge(num_particles)  = qin
                rotation_angle(num_particles) = 0.0d0
                if (rin > rmax_out) rmax_out = rin
                if (rin < rmin_val) rmin_val = rin
            end do
            close(21)
            if (num_particles <= 0) then
                write(*,*) 'エラー: 明示座標ファイルに有効な粒子がありません: ', trim(positions_file)
                stop 'fposit_sub: 有効粒子なし'
            end if
        else
            ! 乱数配置モード
            
            ! 半径リストファイルの読み込み（ファイル指定の場合）
            radius_list_count = 0
            radius_list_idx = 0
            if (trim(radius_distribution_type) == 'file') then
                allocate(radius_list(ni_max))
                call load_particle_radii_from_file(trim(radius_list_file), radius_list, ni_max, radius_list_count)
                write(*,*) '半径リストファイルから', radius_list_count, '個の半径を読み込みました'
            end if
            
            ! 電荷リストファイルの読み込み（ファイル指定の場合）
            charge_list_count = 0
            charge_list_idx = 0
            if (trim(charge_distribution_type) == 'file') then
                allocate(charge_list(ni_max))
                call load_particle_charges_from_file(trim(charge_list_file), charge_list, ni_max, charge_list_count)
                write(*,*) '電荷リストファイルから', charge_list_count, '個の電荷を読み込みました'
            end if
            
            ! 分布タイプに応じた情報表示
            write(*,*) '半径分布タイプ: ', trim(radius_distribution_type)
            write(*,*) '電荷分布タイプ: ', trim(charge_distribution_type)
            
            ! rmax_out と rmin_val の初期推定（セルサイズ計算用）
            select case (trim(radius_distribution_type))
                case ('fixed')
                    rmax_out = max(r1_val, r2_val)
                    rmin_val = min(r1_val, r2_val)
                case ('file')
                    if (radius_list_count > 0) then
                        rmax_out = maxval(radius_list(1:radius_list_count))
                        rmin_val = minval(radius_list(1:radius_list_count))
                    else
                        rmax_out = r1_val
                        rmin_val = r2_val
                    end if
                case ('uniform')
                    rmax_out = radius_uniform_max
                    rmin_val = radius_uniform_min
                case ('normal')
                    ! 正規分布の場合、3σ範囲を考慮
                    rmax_out = radius_normal_mean + 3.0d0 * radius_normal_std
                    rmin_val = max(1.0d-6, radius_normal_mean - 3.0d0 * radius_normal_std)
                case ('lognormal')
                    ! 対数正規分布の場合、3σ範囲を考慮
                    rmax_out = radius_normal_mean + 3.0d0 * radius_normal_std
                    rmin_val = max(1.0d-6, radius_normal_mean - 3.0d0 * radius_normal_std)
                case ('exponential')
                    ! 指数分布の場合、5倍の平均値を上限とする
                    rmax_out = radius_exponential_mean * 5.0d0
                    rmin_val = 1.0d-6
                case ('bimodal')
                    ! 二峰性分布の場合、両ピークの3σ範囲を考慮
                    rmax_out = max(radius_bimodal_mean1 + 3.0d0 * radius_bimodal_std1, &
                                   radius_bimodal_mean2 + 3.0d0 * radius_bimodal_std2)
                    rmin_val = max(1.0d-6, min(radius_bimodal_mean1 - 3.0d0 * radius_bimodal_std1, &
                                               radius_bimodal_mean2 - 3.0d0 * radius_bimodal_std2))
                case default
                    rmax_out = r1_val
                    rmin_val = r2_val
            end select
            
            ! 平均半径を使用してより密な格子を作成
            r_avg = (rmax_out + rmin_val) / 2.0d0
            rn_val = r_avg + 1.0d-5    ! パッキングのための有効半径（平均半径ベース）
            
            ! 壁引き抜きが有効かつ引き抜き対象の斜面壁が指定されている場合、
            ! 粒子配置のx範囲を壁のx座標までに制限
            if (enable_wall_withdraw .and. withdraw_sloped_wall_id > 0 &
                .and. withdraw_sloped_wall_id <= num_walls) then
                effective_container_width = wall_x_start(withdraw_sloped_wall_id)
                write(*,*) '壁引き抜きモード: 粒子配置x範囲を0〜', effective_container_width, 'に制限'
            else
                effective_container_width = container_width
            end if
            
            ipx_calc = idint(effective_container_width / (2.0d0 * rn_val)) ! 1行あたりの粒子数 (概算)
            
            ! z座標範囲を0.3〜0.6に制限するために層数を調整
            z_min_target = 0.0d0
            z_max_target = 0.5d0
            z_range = z_max_target - z_min_target
            layer_height = 2.0d0 * rn_val
            max_layers_for_range = idint(z_range / layer_height)
            if (max_layers_for_range < 1) max_layers_for_range = 1
            
            ! 元の層数を保存
            original_layers = particle_gen_layers
            
            ! particle_gen_layersを制限範囲に合わせる
            if (particle_gen_layers > max_layers_for_range) then
                write(*,*) '警告: 粒子生成層数を', particle_gen_layers, 'から', max_layers_for_range, &
                          'に制限しました（z=0.3-0.6範囲のため）'
                particle_gen_layers = max_layers_for_range
            end if
            
            write(*,*) '衝突チェック付き粒子配置モード'
            write(*,*) '  平均半径: ', r_avg, ' m'
            write(*,*) '  最大半径: ', rmax_out, ' m'
            write(*,*) '  最小半径: ', rmin_val, ' m'
            
            current_particle_count = 0
            collision_skip_count = 0
            max_retry = 10  ! 位置調整の最大リトライ回数
            
            do i_layer = 1, particle_gen_layers
                if (mod(i_layer, 2) == 0) then  ! 偶数層
                    dx_offset = 2.0d0 * rn_val
                    particles_this_row = ipx_calc - 1
                else                            ! 奇数層
                    dx_offset = rn_val
                    particles_this_row = ipx_calc
                end if
                do j_particle_in_layer = 1, particles_this_row
                    call custom_random(random_seed, random_uniform_val)
                    if (random_uniform_val < 2.0d-1) cycle ! 一部の位置をスキップ
                    
                    ! === 先に半径を決定 ===
                    select case (trim(radius_distribution_type))
                        case ('fixed')
                            call custom_random(random_seed, random_uniform_val)
                            if (random_uniform_val < 0.5d0) then
                                generated_radius = r1_val
                            else
                                generated_radius = r2_val
                            end if
                        case ('file')
                            radius_list_idx = radius_list_idx + 1
                            if (radius_list_idx > radius_list_count) then
                                write(*,*) '警告: 半径リストの粒子数が不足しています。配置を打ち切ります。'
                                write(*,*) '  半径リスト数: ', radius_list_count, ', 配置試行数: ', radius_list_idx
                                goto 100
                            end if
                            generated_radius = radius_list(radius_list_idx)
                        case ('uniform')
                            call generate_uniform_random(random_seed, radius_uniform_min, radius_uniform_max, generated_radius)
                        case ('normal')
                            call generate_normal_random(random_seed, radius_normal_mean, radius_normal_std, generated_radius)
                            if (generated_radius < 1.0d-6) generated_radius = 1.0d-6
                        case ('lognormal')
                            call generate_lognormal_random(random_seed, radius_normal_mean, radius_normal_std, generated_radius)
                            if (generated_radius < 1.0d-6) generated_radius = 1.0d-6
                        case ('exponential')
                            call generate_exponential_random(random_seed, radius_exponential_mean, generated_radius)
                            if (generated_radius < 1.0d-6) generated_radius = 1.0d-6
                        case ('bimodal')
                            call generate_bimodal_random(random_seed, radius_bimodal_mean1, radius_bimodal_std1, &
                                                        radius_bimodal_mean2, radius_bimodal_std2, &
                                                        radius_bimodal_ratio, generated_radius)
                            if (generated_radius < 1.0d-6) generated_radius = 1.0d-6
                        case default
                            generated_radius = r1_val
                    end select
                    
                    ! === 候補位置を計算 ===
                    x_candidate = 2.0d0 * rn_val * (j_particle_in_layer - 1) + dx_offset
                    z_candidate = z_min_target + 2.0d0 * rn_val * (i_layer - 1) + rn_val
                    
                    ! === 衝突チェック付き配置 ===
                    placement_success = .false.
                    do retry_count = 0, max_retry
                        ! 衝突チェック
                        call check_particle_collision_sub(x_candidate, z_candidate, generated_radius, &
                                                         current_particle_count, has_collision)
                        if (.not. has_collision) then
                            placement_success = .true.
                            exit
                        end if
                        
                        ! 衝突があった場合、位置を微調整
                        if (retry_count < max_retry) then
                            call custom_random(random_seed, random_uniform_val)
                            x_candidate = x_candidate + (random_uniform_val - 0.5d0) * generated_radius
                            call custom_random(random_seed, random_uniform_val)
                            z_candidate = z_candidate + (random_uniform_val - 0.5d0) * generated_radius
                            
                            ! 境界チェック
                            if (x_candidate < generated_radius) x_candidate = generated_radius
                            if (x_candidate > effective_container_width - generated_radius) &
                                x_candidate = effective_container_width - generated_radius
                            if (z_candidate < z_min_target + generated_radius) &
                                z_candidate = z_min_target + generated_radius
                            if (z_candidate > z_max_target - generated_radius) &
                                z_candidate = z_max_target - generated_radius
                        end if
                    end do
                    
                    ! 配置失敗の場合はスキップ
                    if (.not. placement_success) then
                        collision_skip_count = collision_skip_count + 1
                        cycle
                    end if
                    
                    ! === 粒子を配置 ===
                    current_particle_count = current_particle_count + 1
                    if (current_particle_count > ni_max) then
                        write(*,*) '粒子数がni_maxを超えました: ', ni_max
                        stop 'fposit_sub: 粒子が多すぎます'
                    end if
                    num_particles = current_particle_count
                    x_coord(num_particles) = x_candidate
                    z_coord(num_particles) = z_candidate
                    rotation_angle(num_particles) = 0.0d0
                    radius(num_particles) = generated_radius

                    ! 電荷の割り当て（分布タイプに応じて）
                    select case (trim(charge_distribution_type))
                        case ('fixed')
                            charge(num_particles) = default_charge
                        case ('file')
                            charge_list_idx = charge_list_idx + 1
                            if (charge_list_idx > charge_list_count) then
                                charge_list_idx = mod(charge_list_idx - 1, charge_list_count) + 1
                            end if
                            charge(num_particles) = charge_list(charge_list_idx)
                        case ('uniform')
                            call generate_uniform_random(seed_charge, charge_uniform_min, charge_uniform_max, generated_charge)
                            charge(num_particles) = generated_charge
                        case ('normal')
                            call generate_normal_random(seed_charge, charge_normal_mean, charge_normal_std, generated_charge)
                            charge(num_particles) = generated_charge
                        case ('lognormal')
                            call generate_lognormal_random(seed_charge, charge_normal_mean, charge_normal_std, generated_charge)
                            charge(num_particles) = generated_charge
                        case ('exponential')
                            call generate_exponential_random(seed_charge, charge_exponential_mean, generated_charge)
                            charge(num_particles) = generated_charge
                        case ('bimodal')
                            call generate_bimodal_random(seed_charge, charge_bimodal_mean1, charge_bimodal_std1, &
                                                        charge_bimodal_mean2, charge_bimodal_std2, &
                                                        charge_bimodal_ratio, generated_charge)
                            charge(num_particles) = generated_charge
                        case default
                            charge(num_particles) = default_charge
                    end select
                end do
            end do
100         continue  ! 半径リスト不足時のジャンプ先
            
            if (collision_skip_count > 0) then
                write(*,*) '衝突により配置をスキップした粒子数: ', collision_skip_count
            end if
            
            ! 半径リスト数 > 配置可能数の警告
            if (trim(radius_distribution_type) == 'file' .and. radius_list_count > 0) then
                if (radius_list_count > num_particles) then
                    write(*,*) '========================================'
                    write(*,*) '警告: 半径リストの粒子数が配置可能数を超えています'
                    write(*,*) '  半径リスト数    : ', radius_list_count
                    write(*,*) '  配置された粒子数: ', num_particles
                    write(*,*) '  未使用の半径数  : ', (radius_list_count - num_particles)
                    write(*,*) '========================================'
                end if
            end if
            
            ! 配列の解放
            if (allocated(radius_list)) deallocate(radius_list)
            if (allocated(charge_list)) deallocate(charge_list)
            
            ! 実際に生成された粒子の半径から rmax_out, rmin_val を再計算
            if (num_particles > 0) then
                rmax_out = maxval(radius(1:num_particles))
                rmin_val = minval(radius(1:num_particles))
            end if
        end if
        write(*,*) '生成された粒子数: ', num_particles

        ! セルサイズ計算 (原文PDF p.35 eq 3.25: C < sqrt(2)*rmin)
        if (disable_cell_algorithm .or. cell_size_override > 0.0d0) then
            ! セル法アルゴリズムを無効化する場合
            if (cell_size_override > 0.0d0) then
                cell_size = cell_size_override
            else
                ! 計算領域全体を1つのセルとして設定
                cell_size = max(container_width, 2.0d0 * rmax_out * particle_gen_layers) * 2.0d0
            end if
            write(*,*) 'セル法アルゴリズムを無効化: cell_size = ', cell_size
        else
            ! 通常のセルサイズ計算
            if (rmin_val > 0.0d0) then
                 cell_size = rmin_val * 1.30d0 ! または入力から。原文ではrmin*1.35d0はコメントアウト
            else
                 cell_size = rmax_out * 1.30d0 ! rminが適切に定義されない場合のフォールバック
            end if
        end if
        
        if (cell_size <= 0.0d0) then
            write(*,*) "エラー: fposit_subでcell_sizeが正ではありません。"
            stop
        endif

        cells_x_dir = idint(container_width / cell_size) + 1
        
        ! cells_z_dirは粒子が到達しうる最大高さをカバーする必要がある
        ! container_heightが指定されている場合はそれを使用、そうでなければ従来の推定方法を使用
        if (container_height > 0.0d0 .and. cell_size > 0.0d0) then
            ! 上壁が指定されている場合は、その高さでセル数を計算
            cells_z_dir = idint(container_height / cell_size) + 1
        else if (num_particles > 0 .and. cell_size > 0.0d0) then
            ! 従来の方法: 生成時の最上部粒子のz座標を使用
            cells_z_dir = idint(z_coord(num_particles) / cell_size) + 10 
        else if (particle_gen_layers > 0 .and. cell_size > 0.0d0 .and. rn_val > 0.0d0) then ! 粒子がない場合でも推定
             if (particle_gen_layers > 0 .and. rn_val > 0 .and. cell_size > 0) then
                cells_z_dir = idint( (2.0d0 * rn_val * (real(particle_gen_layers) -1.0d0) + rn_val) / cell_size) + 10
             else
                cells_z_dir = 20 ! Fallback if values are still problematic
             end if
        else
            cells_z_dir = 20 ! デフォルト値 (粒子も層もない、またはセルサイズが0の場合)
        end if


        if (cells_x_dir * cells_z_dir > nc_max) then
            write(*,*) 'ncl (cell_particle_map)がオーバーフローしました!! 要求セル数: ', cells_x_dir * cells_z_dir
            stop 'fposit_sub: セル配列が小さすぎます'
        end if

        call profiler_stop('fposit_sub', prof_token)
    end subroutine fposit_sub

    !> 材料物性値を初期化し、定数を計算するサブルーチン
    subroutine inmat_sub
        use simulation_constants_mod, only: PI_VAL
        use simulation_parameters_mod, only: time_step, reference_overlap
        use particle_data_mod, only: radius, mass, moment_inertia 
        use cell_system_mod, only: num_particles 
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        integer :: i
        integer(profiler_tick_kind) :: prof_token
        
        call profiler_start('inmat_sub', prof_token)
        ! time_step, particle_densityなどの値はモジュールsimulation_parameters_modで設定されていると仮定
        ! 粒子のポアソン比に基づいてせん断弾性係数と法線方向弾性係数の比(so)を計算
        shear_to_normal_stiffness_ratio = 1.0d0 / (2.0d0 * (1.0d0 + poisson_ratio_particle))

        !$omp parallel do private(i)
        do i = 1, num_particles
            ! 質量: 3D球体 V = 4/3 pi r^3。2Dディスク (面積 pi r^2)の場合、deが面密度ならば。
            ! 元のコードは2Dシミュレーションの文脈でも3D球体の体積で質量を計算しているように見える。
            mass(i) = (4.0d0 / 3.0d0) * PI_VAL * radius(i)**3 * particle_density
            
            ! 慣性モーメント: 3D球体 I = 2/5 m r^2。
            ! 元のコード: pmi(i)=8.d0/15.d0*de*pi*(rr(i)**5)
            ! これは (2/5) * (4/3 pi r^3 de) * r^2 = (2/5) * mass * r^2。球体として正しい。
            moment_inertia(i) = (8.0d0 / 15.0d0) * particle_density * PI_VAL * (radius(i)**5)

            ! 参照食い込み量 δ_ref を平均半径の5%で設定
            if (num_particles > 0) then
                reference_overlap = 0.05d0 * sum(radius(1:num_particles)) / real(num_particles, 8)
            else
                reference_overlap = 0.0d0
            end if  

        end do
        !$omp end parallel do
        call profiler_stop('inmat_sub', prof_token)
    end subroutine inmat_sub

    !> 接触力関連の配列を初期化するサブルーチン
    subroutine init_sub
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod, only: num_particles
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        integer :: i, j
        integer(profiler_tick_kind) :: prof_token

        call profiler_start('init_sub', prof_token)

        ! 実際の粒子数まで繰り返す
        if (num_particles > 0) then
            !$omp parallel do private(i, j)
            do i = 1, num_particles 
                do j = 1, nj_max
                    normal_force_contact(i, j) = 0.0d0
                    shear_force_contact(i, j) = 0.0d0
                    contact_partner_idx(i, j) = 0
                    previous_overlap(i, j) = -1.0d0  ! 負値で「接触なし」を表現
                end do
            end do
            !$omp end parallel do
            ! 増分変位を最初にゼロで初期化
            x_disp_incr(1:num_particles) = 0.0d0
            z_disp_incr(1:num_particles) = 0.0d0
            rot_disp_incr(1:num_particles) = 0.0d0
        end if
        call profiler_stop('init_sub', prof_token)
    end subroutine init_sub

    !> 近傍探索のために粒子をセルに割り当てるサブルーチン
    subroutine ncel_sub
        use simulation_constants_mod, only: nc_max
        use particle_data_mod, only: x_coord, z_coord
        use cell_system_mod
        implicit none
        integer :: i, cell_block_idx
        integer :: ix_cell, iz_cell ! 宣言をここに移動
    
        ! 連結リストをクリア
        if (nc_max > 0) cell_head(1:nc_max) = 0
        if (nc_max > 0) cell_particle_map(1:nc_max) = 0  ! デバッグ用途
        if (num_particles > 0) particle_cell_next(1:num_particles) = 0

        do i = 1, num_particles
            particle_cell_idx(i) = 0 ! 初期化
            if (cell_size <= 0.0d0) then
                write(*,*) "エラー: ncel_subでcell_sizeが0または負です。"
                stop
            endif
            if (cells_x_dir <= 0) then
                 write(*,*) "エラー: ncel_subでcells_x_dirが0または負です。"
                 stop
            endif

            ! 粒子iを含むセルの1次元インデックスを計算
            ix_cell = idint(x_coord(i) / cell_size) ! 座標が負にならないように注意
            iz_cell = idint(z_coord(i) / cell_size)

            ! インデックスが有効範囲 [0, cells_x_dir-1] および [0, cells_z_dir-1] 内にあることを保証
            ix_cell = max(0, min(ix_cell, cells_x_dir - 1))
            iz_cell = max(0, min(iz_cell, cells_z_dir - 1))
            
            cell_block_idx = iz_cell * cells_x_dir + ix_cell + 1 ! Fortranの1ベースインデックス

            if (cell_block_idx > 0 .and. cell_block_idx <= nc_max) then
                 ! 連結リストの先頭に追加
                 particle_cell_next(i) = cell_head(cell_block_idx)
                 cell_head(cell_block_idx) = i

                 ! 旧 single-map も更新（最後に登録された粒子）
                 cell_particle_map(cell_block_idx) = i

                 particle_cell_idx(i) = cell_block_idx
            else
                write(*,*) 'エラー: ncel_subで粒子', i, 'のcell_block_idxが範囲外です。'
                write(*,*) 'x0, z0, c: ', x_coord(i), z_coord(i), cell_size
                write(*,*) 'ix_cell, iz_cell, idx, computed_block_idx: ', ix_cell, iz_cell, cells_x_dir, cell_block_idx
                stop 'ncel_sub: 粒子が無効なセルインデックスにマッピングされました'
            end if
        end do
    end subroutine ncel_sub

    !> 粒子iと壁との接触力を計算するサブルーチン
    subroutine wcont_sub(particle_idx)
        use simulation_parameters_mod, only: container_width, container_height, &
            left_wall_active, bottom_wall_active, right_wall_active, top_wall_active
        use particle_data_mod
        use cell_system_mod, only: num_particles
        use wall_data_mod
        use simulation_constants_mod, only: nj_max
        implicit none
        integer, intent(in) :: particle_idx ! 対象の粒子インデックス
    
        real(8) :: xi, zi, ri_particle ! 粒子iのx座標, z座標, 半径
        real(8) :: wall_angle_sin, wall_angle_cos, overlap_gap ! 壁の法線ベクトル成分, 重なり量
        integer :: wall_contact_slot_idx, wall_partner_id ! 壁の接触スロットインデックス, 壁の相手粒子インデックス
        integer :: wall_idx, slot_idx, first_sloped_slot
        real(8) :: proj_len, closest_x, closest_z
        real(8) :: vec_x, vec_z, dist_sq, dist_val
        real(8) :: dynamic_angle_sin, dynamic_angle_cos
        logical, save :: wall_slot_warning_emitted = .false.
        real(8) :: en_coeff, et_coeff ! 反発係数 (normal, tangential)

        xi = x_coord(particle_idx)
        zi = z_coord(particle_idx)
        ri_particle = radius(particle_idx)
        first_sloped_slot = 14
    
        ! 左壁 (contact_partner_idx = num_particles + 1)
        wall_contact_slot_idx = 11 ! 元のコードでの左壁用の固定スロット
        wall_partner_id = num_particles + 1
        if (left_wall_active .and. xi < ri_particle) then  ! 左壁と接触
            en_coeff = 0.9d0
            et_coeff = 0.9d0
            wall_angle_sin = 0.0d0  ! 法線ベクトル成分 sin(alpha_ij) (粒子中心から壁中心へ向かうベクトル)
            wall_angle_cos = -1.0d0 ! 法線ベクトル成分 cos(alpha_ij)
            overlap_gap = ri_particle - xi ! 元のコードでは dabs(xi)、ここでは重なり量を正とする
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                        ! 接触なし（または壁が無効）
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 下壁 (contact_partner_idx = num_particles + 2)
        wall_contact_slot_idx = 12 ! 元のコードでの下壁用の固定スロット
        wall_partner_id = num_particles + 2
        if (bottom_wall_active .and. zi < ri_particle) then  ! 下壁と接触
            en_coeff = 0.9d0
            et_coeff = 0.9d0
            wall_angle_sin = -1.0d0
            wall_angle_cos = 0.0d0
            overlap_gap = ri_particle - zi 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                        ! 接触なし（または壁が無効）
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 右壁 (contact_partner_idx = num_particles + 3)
        wall_contact_slot_idx = 13 ! 元のコードでの右壁用の固定スロット
        wall_partner_id = num_particles + 3
        if (right_wall_active .and. xi + ri_particle > container_width) then ! 右壁と接触
            en_coeff = 0.9d0
            et_coeff = 0.9d0
            wall_angle_sin = 0.0d0
            wall_angle_cos = 1.0d0
            overlap_gap = (xi + ri_particle) - container_width 
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
            call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
        else                                        ! 接触なし（または壁が無効）
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 上壁 (contact_partner_idx = num_particles + 4)
        ! container_height > 0の場合のみ上壁を有効化
        if (container_height > 0.0d0 .and. top_wall_active) then
            wall_contact_slot_idx = 10 ! 上壁用の固定スロット
            wall_partner_id = num_particles + 4
            if (zi + ri_particle > container_height) then ! 上壁と接触
                en_coeff = 0.9d0
                et_coeff = 0.9d0
                wall_angle_sin = 1.0d0
                wall_angle_cos = 0.0d0
                overlap_gap = (zi + ri_particle) - container_height
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
                call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, wall_angle_sin, wall_angle_cos, overlap_gap, en_coeff, et_coeff)
            else                                            ! 接触なし
                normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
            end if
        else if (container_height > 0.0d0 .and. .not. top_wall_active) then
            ! 上壁が無効化された場合、接触スロットをクリア
            wall_contact_slot_idx = 10
            normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
            contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
        end if
    
        ! 斜面壁 (任意本数、ファイルで定義)
        if (num_walls > 0 .and. first_sloped_slot <= nj_max) then
            do wall_idx = 1, num_walls
                ! 無効な壁はスキップ
                if (.not. wall_active(wall_idx)) then
                    ! 既存の接触スロットがあればクリア
                    wall_partner_id = num_particles + 4 + wall_idx
                    do slot_idx = first_sloped_slot, nj_max
                        if (contact_partner_idx(particle_idx, slot_idx) == wall_partner_id) then
                            normal_force_contact(particle_idx, slot_idx) = 0.0d0
                            shear_force_contact(particle_idx, slot_idx) = 0.0d0
                            contact_partner_idx(particle_idx, slot_idx) = 0
                            exit
                        end if
                    end do
                    cycle
                end if
                
                wall_partner_id = num_particles + 4 + wall_idx
                en_coeff = 0.9d0
                et_coeff = 0.9d0
                proj_len = (xi - wall_x_start(wall_idx)) * wall_tangent_x(wall_idx) + &
                           (zi - wall_z_start(wall_idx)) * wall_tangent_z(wall_idx)
                proj_len = max(0.0d0, min(proj_len, wall_length(wall_idx)))
                closest_x = wall_x_start(wall_idx) + wall_tangent_x(wall_idx) * proj_len
                closest_z = wall_z_start(wall_idx) + wall_tangent_z(wall_idx) * proj_len
    
                vec_x = closest_x - xi
                vec_z = closest_z - zi
                dist_sq = vec_x * vec_x + vec_z * vec_z
    
                if (dist_sq > 1.0d-20) then
                    dist_val = sqrt(dist_sq)
                    dynamic_angle_cos = vec_x / dist_val
                    dynamic_angle_sin = vec_z / dist_val
                else
                    dist_val = 0.0d0
                    dynamic_angle_cos = -wall_normal_x(wall_idx)
                    dynamic_angle_sin = -wall_normal_z(wall_idx)
                end if
    
                overlap_gap = ri_particle - dist_val
    
                wall_contact_slot_idx = 0
                do slot_idx = first_sloped_slot, nj_max
                    if (contact_partner_idx(particle_idx, slot_idx) == wall_partner_id) then
                        wall_contact_slot_idx = slot_idx
                        exit
                    end if
                end do
    
                if (overlap_gap > 0.0d0) then
                    if (wall_contact_slot_idx == 0) then
                        do slot_idx = first_sloped_slot, nj_max
                            if (contact_partner_idx(particle_idx, slot_idx) == 0) then
                                wall_contact_slot_idx = slot_idx
                                exit
                            end if
                        end do
                    end if
    
                    if (wall_contact_slot_idx == 0) then
                        if (.not. wall_slot_warning_emitted) then
                            write(*,*) '警告: 斜面壁との接触を格納するスロットが不足しています (粒子', particle_idx, ')'
                            wall_slot_warning_emitted = .true.
                        end if
                        cycle
                    end if
    
                    contact_partner_idx(particle_idx, wall_contact_slot_idx) = wall_partner_id
                    call actf_sub(particle_idx, wall_partner_id, wall_contact_slot_idx, dynamic_angle_sin, dynamic_angle_cos, &
                                  overlap_gap, en_coeff, et_coeff)
                else
                    if (wall_contact_slot_idx /= 0) then
                        normal_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                        shear_force_contact(particle_idx, wall_contact_slot_idx) = 0.0d0
                        contact_partner_idx(particle_idx, wall_contact_slot_idx) = 0
                    end if
                end if
            end do
        end if
    end subroutine wcont_sub

    !> 粒子iと他の粒子との接触力を計算するサブルーチン
    subroutine pcont_sub(particle_i_idx, rmax_val)
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod
        implicit none

        integer, intent(in) :: particle_i_idx     ! 対象の粒子iのインデックス
        real(8), intent(in) :: rmax_val           ! 最大粒子半径 (探索範囲に使用)

        real(8) :: xi, zi, ri_particle_i
        real(8) :: xj, zj, rj_particle_j
        real(8) :: center_dist, overlap_gap
        real(8) :: contact_angle_sin, contact_angle_cos ! 粒子iから粒子jへの法線ベクトル成分
        integer :: particle_j_idx                 ! 接触相手の粒子jのインデックス
        integer :: contact_slot_for_i_j, contact_slot_for_j_i ! 接触リストのスロット
        integer :: iz_cell_min, iz_cell_max, ix_cell_min, ix_cell_max ! セル探索範囲
        integer :: current_iz_cell, current_ix_cell, cell_block_idx
        integer :: jj, max_particle_contacts_check
        real(8) :: dx, dz               ! 位置差計算用
        integer :: curr_j_idx
        real(8) :: search_extent        ! 近傍探索半径

        max_particle_contacts_check = 10 ! 元のコードのループ(do 11 jj=1,10)から、粒子間接触は最大10個と仮定

        xi = x_coord(particle_i_idx)
        zi = z_coord(particle_i_idx)
        ri_particle_i = radius(particle_i_idx)

        ! セル格子内での探索範囲を決定
        if (cell_size <= 0.0d0) stop "pcont_sub: cell_sizeが正ではありません。"

        ! 以前は ±(2*rmax) で探索していたが、cell_size が大きい場合に隣接セル2つ分を
        ! またぐ衝突を取りこぼすことがあった。そこで探索半径を (2*rmax + cell_size) に拡張する。
        search_extent = 2.0d0 * rmax_val + cell_size

        iz_cell_max = idint((zi + search_extent) / cell_size)
        iz_cell_min = idint((zi - search_extent) / cell_size)
        ix_cell_min = idint((xi - search_extent) / cell_size)
        ix_cell_max = idint((xi + search_extent) / cell_size)

        ! セルインデックスを有効なグリッド境界内に収める
        iz_cell_max = min(iz_cell_max, cells_z_dir - 1)
        iz_cell_min = max(iz_cell_min, 0)
        ix_cell_min = max(ix_cell_min, 0)
        ix_cell_max = min(ix_cell_max, cells_x_dir - 1)
        
        if (iz_cell_max < iz_cell_min .or. ix_cell_max < ix_cell_min) then
             ! 粒子が通常の領域外にある場合や、rmax_valが小さすぎて探索範囲が無効になる場合に発生しうる
             return
        end if

        ! デバッグ出力 (検証モードで粒子1のみ)
        ! if (validation_mode .and. particle_i_idx == 1) then
        !     write(*,'(A,I2,A,I3,A,I3,A,I3,A,I3)') 'Search p',particle_i_idx,': cells x[',ix_cell_min,':',ix_cell_max,'] z[',iz_cell_min,':',iz_cell_max,']'
        !     write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') 'Position: x=',xi,', z=',zi,', search_extent=',search_extent
        ! end if

        do current_iz_cell = iz_cell_min, iz_cell_max      ! z方向のセルループ
            do current_ix_cell = ix_cell_min, ix_cell_max  ! x方向のセルループ
                if (cells_x_dir <=0) stop "pcont_sub: cells_x_dirが正ではありません。"
                cell_block_idx = current_iz_cell * cells_x_dir + current_ix_cell + 1

                if (cell_block_idx <= 0 .or. cell_block_idx > (cells_x_dir * cells_z_dir) ) cycle ! セルが範囲外ならスキップ

                curr_j_idx = cell_head(cell_block_idx) ! セルに属する最初の粒子

                do while (curr_j_idx > 0)
                    particle_j_idx = curr_j_idx

                    if (particle_j_idx == particle_i_idx) then
                        curr_j_idx = particle_cell_next(curr_j_idx)
                        cycle
                    end if

                    xj = x_coord(particle_j_idx)
                    zj = z_coord(particle_j_idx)
                    rj_particle_j = radius(particle_j_idx)

                    dx = xi - xj
                    dz = zi - zj
                    center_dist = sqrt(dx*dx + dz*dz) 
                    overlap_gap = (ri_particle_i + rj_particle_j) - center_dist ! 重なり量 (正なら接触)

                    if (overlap_gap > 0.0d0) then ! 粒子が接触している (重なっている)
                        ! デバッグ出力 (検証モードのみ)
                        ! if (validation_mode) then
                        !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        !         'Contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        ! end if

                        if (center_dist < 1.0d-12) then ! 粒子中心が完全に一致する場合のゼロ除算を回避
                            contact_angle_cos = 1.0d0   ! 暫定的にx軸方向とする
                            contact_angle_sin = 0.0d0
                        else
                            ! 法線ベクトル (i から j へ向かう方向) の成分
                            contact_angle_cos = (xj - xi) / center_dist 
                            contact_angle_sin = (zj - zi) / center_dist
                        end if
                        
                        ! ---- 連絡スロット確保 (i -> j) ----------------
                        contact_slot_for_i_j = 0
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                contact_slot_for_i_j = jj
                                exit
                            end if
                        end do
                        if (contact_slot_for_i_j == 0) then
                            do jj = 1, max_particle_contacts_check
                                if (contact_partner_idx(particle_i_idx, jj) == 0) then
                                    contact_slot_for_i_j = jj
                                    contact_partner_idx(particle_i_idx, jj) = particle_j_idx
                                    exit
                                end if
                            end do
                        end if
                        if (contact_slot_for_i_j == 0) then
                            curr_j_idx = particle_cell_next(curr_j_idx)
                            cycle  ! スロット不足
                        end if

                        ! -------------------------------------------------

                        ! 実際の力計算
                        call actf_sub(particle_i_idx, particle_j_idx, contact_slot_for_i_j, &
                                      contact_angle_sin, contact_angle_cos, overlap_gap, 0.9d0, 0.9d0)
                    else ! 幾何学的な接触なし / 粒子が離れた
                        ! デバッグ出力 (検証モードのみ、最初の数回のみ)
                        ! if (validation_mode .and. particle_i_idx == 1 .and. particle_j_idx == 2) then
                        !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        !         'No contact p',particle_i_idx,' <-> p',particle_j_idx,': dist=',center_dist,', overlap=',overlap_gap,', sumr=',ri_particle_i+rj_particle_j
                        ! end if

                        ! particle_i_idx の particle_j_idx に関する接触情報をクリア
                        do jj = 1, max_particle_contacts_check
                            if (contact_partner_idx(particle_i_idx, jj) == particle_j_idx) then
                                normal_force_contact(particle_i_idx, jj) = 0.0d0
                                shear_force_contact(particle_i_idx, jj) = 0.0d0
                                contact_partner_idx(particle_i_idx, jj) = 0
                                exit
                            end if
                        end do
                    end if

                    curr_j_idx = particle_cell_next(curr_j_idx) ! セル内の次の粒子へ
                end do !! セル内連結リスト走査
            end do ! x方向セルループ
        end do     ! z方向セルループ
    end subroutine pcont_sub

    !> 全粒子間のクーロン力を計算するサブルーチン
    subroutine coulomb_force_sub
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        
        integer :: i, j
        real(8) :: dx, dz, dist, dist_sq, dist_cubed
        real(8) :: force_magnitude, fx, fz
        real(8) :: qi, qj
        
        ! クーロン力が無効化されている場合は何もしない
        if (.not. enable_coulomb_force) then
            return
        end if
        
        ! 全粒子ペアについてクーロン力を計算
        !$omp parallel do schedule(dynamic, 16) private(i, j, qi, qj, dx, dz, dist_sq, dist, dist_cubed, force_magnitude, fx, fz)
        do i = 1, num_particles - 1
            qi = charge(i)
            if (abs(qi) < 1.0d-20) cycle ! 電荷がゼロならスキップ
            
            do j = i + 1, num_particles
                qj = charge(j)
                if (abs(qj) < 1.0d-20) cycle ! 電荷がゼロならスキップ
                
                ! 粒子間の距離ベクトルと距離
                dx = x_coord(j) - x_coord(i)
                dz = z_coord(j) - z_coord(i)
                dist_sq = dx*dx + dz*dz
                
                ! ゼロ除算を回避
                if (dist_sq < 1.0d-20) cycle
                
                dist = sqrt(dist_sq)
                dist_cubed = dist * dist_sq
                
                ! クーロン力の大きさ F = k * q1 * q2 / r²
                ! 力のベクトル成分 F_vec = F * (r_vec / |r|) = k * q1 * q2 * r_vec / r³
                force_magnitude = coulomb_constant * qi * qj / dist_cubed
                
                fx = force_magnitude * dx
                fz = force_magnitude * dz
                
                ! 粒子iに力を加算（反発力なので i<-j 方向、つまり -dx 方向）
                !$omp atomic
                x_force_sum(i) = x_force_sum(i) - fx
                !$omp atomic
                z_force_sum(i) = z_force_sum(i) - fz
                
                ! 粒子jに反作用力を加算（反発力なので j->i 方向、つまり +dx 方向）
                !$omp atomic
                x_force_sum(j) = x_force_sum(j) + fx
                !$omp atomic
                z_force_sum(j) = z_force_sum(j) + fz
            end do
        end do
        !$omp end parallel do
    end subroutine coulomb_force_sub

    !> セル法＋カットオフを用いたクーロン力計算サブルーチン（O(N・近傍)）
    !> 既存のセル構造（cell_head, particle_cell_next）を流用し、
    !> カットオフ半径rc以内の粒子ペアのみを計算する。
    subroutine coulomb_force_cell_cutoff_sub
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force, &
            coulomb_cutoff, coulomb_softening, coulomb_shift_force
        use particle_data_mod
        use cell_system_mod
        implicit none
        
        integer :: i, j
        integer :: ix, iz, ix2, iz2, dxc, dzc
        integer :: cell0, cell1, nspan
        real(8) :: dx, dz, dist_sq, r_eff_sq, r_eff_cubed
        real(8) :: force_factor, fx, fz
        real(8) :: qi, qj
        real(8) :: rc, rc_sq, delta_sq
        real(8) :: rc_eff_sq, rc_eff_cubed  ! softened cutoff用
        real(8), parameter :: eps_charge = 1.0d-20
        
        ! クーロン力が無効化されている場合は何もしない
        if (.not. enable_coulomb_force) return
        
        rc = coulomb_cutoff
        rc_sq = rc * rc
        delta_sq = coulomb_softening * coulomb_softening
        ! shifted-force用: softened cutoff値 (rc^2 + δ^2)^(3/2)
        rc_eff_sq = rc_sq + delta_sq
        rc_eff_cubed = rc_eff_sq * sqrt(rc_eff_sq)
        
        ! カットオフが0以下の場合は旧実装（全対全）にフォールバック
        if (rc <= 0.0d0) then
            call coulomb_force_sub_full_pairs()
            return
        end if
        
        ! 近傍セルの範囲を決定: nspan = ceil(rc / cell_size)
        ! セルサイズが0以下の場合は安全のため1に設定
        if (cell_size > 0.0d0) then
            nspan = ceiling(rc / cell_size)
        else
            nspan = 1
        end if
        
        ! セル走査によるクーロン力計算
        ! 二重カウントを避けるため、同一セル内は j=next(i)、
        ! 異なるセルは「前方向のみ」（dzc>0, または dzc==0 かつ dxc>0）を走査
        !$omp parallel do collapse(2) schedule(dynamic, 16) &
        !$omp private(ix, iz, cell0, i, j, qi, qj, dx, dz, dist_sq, r_eff_sq, r_eff_cubed, force_factor, fx, fz, dxc, dzc, ix2, iz2, cell1)
        do iz = 0, cells_z_dir - 1
            do ix = 0, cells_x_dir - 1
                cell0 = iz * cells_x_dir + ix + 1
                
                ! === 同一セル内ペア ===
                i = cell_head(cell0)
                do while (i > 0)
                    qi = charge(i)
                    if (abs(qi) < eps_charge) then
                        i = particle_cell_next(i)
                        cycle
                    end if
                    
                    ! j = next(i) から開始（i < j を保証）
                    j = particle_cell_next(i)
                    do while (j > 0)
                        qj = charge(j)
                        if (abs(qj) < eps_charge) then
                            j = particle_cell_next(j)
                            cycle
                        end if
                        
                        ! 距離計算
                        dx = x_coord(j) - x_coord(i)
                        dz = z_coord(j) - z_coord(i)
                        dist_sq = dx*dx + dz*dz
                        
                        ! カットオフ判定
                        if (dist_sq < rc_sq) then
                            ! ソフトニング: r_eff^2 = r^2 + δ^2
                            r_eff_sq = dist_sq + delta_sq
                            r_eff_cubed = r_eff_sq * sqrt(r_eff_sq)
                            
                            ! シフトドフォース: F = k*qi*qj*r_vec*(1/r_eff^3 - 1/rc_eff^3)
                            ! ※ softening込みのカットオフ値を使用して r=rc で力が0になるようにする
                            if (coulomb_shift_force) then
                                force_factor = coulomb_constant * qi * qj * (1.0d0/r_eff_cubed - 1.0d0/rc_eff_cubed)
                            else
                                force_factor = coulomb_constant * qi * qj / r_eff_cubed
                            end if
                            
                            fx = force_factor * dx
                            fz = force_factor * dz
                            
                            ! 粒子iに力を加算（反発力なので -dx 方向）
                            !$omp atomic
                            x_force_sum(i) = x_force_sum(i) - fx
                            !$omp atomic
                            z_force_sum(i) = z_force_sum(i) - fz
                            
                            ! 粒子jに反作用力を加算（+dx 方向）
                            !$omp atomic
                            x_force_sum(j) = x_force_sum(j) + fx
                            !$omp atomic
                            z_force_sum(j) = z_force_sum(j) + fz
                        end if
                        
                        j = particle_cell_next(j)
                    end do
                    
                    i = particle_cell_next(i)
                end do
                
                ! === 前方向の近傍セル ===
                ! dzc > 0, または dzc == 0 かつ dxc > 0 のセルのみ走査
                do dzc = 0, nspan
                    do dxc = -nspan, nspan
                        ! 前方向条件: dzc > 0, または (dzc == 0 かつ dxc > 0)
                        if (dzc == 0 .and. dxc <= 0) cycle
                        
                        ix2 = ix + dxc
                        iz2 = iz + dzc
                        
                        ! 境界チェック
                        if (ix2 < 0 .or. ix2 >= cells_x_dir) cycle
                        if (iz2 < 0 .or. iz2 >= cells_z_dir) cycle
                        
                        cell1 = iz2 * cells_x_dir + ix2 + 1
                        
                        ! cell0の全粒子 × cell1の全粒子
                        i = cell_head(cell0)
                        do while (i > 0)
                            qi = charge(i)
                            if (abs(qi) < eps_charge) then
                                i = particle_cell_next(i)
                                cycle
                            end if
                            
                            j = cell_head(cell1)
                            do while (j > 0)
                                qj = charge(j)
                                if (abs(qj) < eps_charge) then
                                    j = particle_cell_next(j)
                                    cycle
                                end if
                                
                                ! 距離計算
                                dx = x_coord(j) - x_coord(i)
                                dz = z_coord(j) - z_coord(i)
                                dist_sq = dx*dx + dz*dz
                                
                                ! カットオフ判定
                                if (dist_sq < rc_sq) then
                                    ! ソフトニング: r_eff^2 = r^2 + δ^2
                                    r_eff_sq = dist_sq + delta_sq
                                    r_eff_cubed = r_eff_sq * sqrt(r_eff_sq)
                                    
                                    ! シフトドフォース（softening込みのカットオフ値を使用）
                                    if (coulomb_shift_force) then
                                        force_factor = coulomb_constant * qi * qj * (1.0d0/r_eff_cubed - 1.0d0/rc_eff_cubed)
                                    else
                                        force_factor = coulomb_constant * qi * qj / r_eff_cubed
                                    end if
                                    
                                    fx = force_factor * dx
                                    fz = force_factor * dz
                                    
                                    ! 粒子iに力を加算
                                    !$omp atomic
                                    x_force_sum(i) = x_force_sum(i) - fx
                                    !$omp atomic
                                    z_force_sum(i) = z_force_sum(i) - fz
                                    
                                    ! 粒子jに反作用力を加算
                                    !$omp atomic
                                    x_force_sum(j) = x_force_sum(j) + fx
                                    !$omp atomic
                                    z_force_sum(j) = z_force_sum(j) + fz
                                end if
                                
                                j = particle_cell_next(j)
                            end do
                            
                            i = particle_cell_next(i)
                        end do
                    end do
                end do
                
            end do
        end do
        !$omp end parallel do
        
    end subroutine coulomb_force_cell_cutoff_sub
    
    !> 全ペア計算版クーロン力（カットオフなし、セル法なし）
    !> coulomb_use_cell=false または coulomb_cutoff<=0 の場合に使用
    subroutine coulomb_force_sub_full_pairs
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force, &
            coulomb_softening, coulomb_shift_force, coulomb_cutoff
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        
        integer :: i, j
        real(8) :: dx, dz, dist_sq, r_eff_sq, r_eff_cubed
        real(8) :: force_factor, fx, fz
        real(8) :: qi, qj
        real(8) :: delta_sq, rc, rc_sq
        real(8) :: rc_eff_sq, rc_eff_cubed  ! softened cutoff用
        real(8), parameter :: eps_charge = 1.0d-20
        real(8), parameter :: eps_dist = 1.0d-20
        
        if (.not. enable_coulomb_force) return
        
        delta_sq = coulomb_softening * coulomb_softening
        rc = coulomb_cutoff
        rc_sq = rc * rc
        ! shifted-force用: softened cutoff値 (rc^2 + δ^2)^(3/2)
        rc_eff_sq = rc_sq + delta_sq
        rc_eff_cubed = rc_eff_sq * sqrt(rc_eff_sq)
        
        ! 全粒子ペアについてクーロン力を計算
        !$omp parallel do schedule(dynamic, 16) private(i, j, qi, qj, dx, dz, dist_sq, r_eff_sq, r_eff_cubed, force_factor, fx, fz)
        do i = 1, num_particles - 1
            qi = charge(i)
            if (abs(qi) < eps_charge) cycle
            
            do j = i + 1, num_particles
                qj = charge(j)
                if (abs(qj) < eps_charge) cycle
                
                dx = x_coord(j) - x_coord(i)
                dz = z_coord(j) - z_coord(i)
                dist_sq = dx*dx + dz*dz
                
                ! カットオフ判定（rc > 0 の場合のみ）
                if (rc > 0.0d0 .and. dist_sq >= rc_sq) cycle
                
                ! ゼロ除算回避
                if (dist_sq < eps_dist) cycle
                
                ! ソフトニング: r_eff^2 = r^2 + δ^2
                r_eff_sq = dist_sq + delta_sq
                r_eff_cubed = r_eff_sq * sqrt(r_eff_sq)
                
                ! シフトドフォース（softening込みのカットオフ値を使用）
                if (coulomb_shift_force .and. rc > 0.0d0) then
                    force_factor = coulomb_constant * qi * qj * (1.0d0/r_eff_cubed - 1.0d0/rc_eff_cubed)
                else
                    force_factor = coulomb_constant * qi * qj / r_eff_cubed
                end if
                
                fx = force_factor * dx
                fz = force_factor * dz
                
                !$omp atomic
                x_force_sum(i) = x_force_sum(i) - fx
                !$omp atomic
                z_force_sum(i) = z_force_sum(i) - fz
                
                !$omp atomic
                x_force_sum(j) = x_force_sum(j) + fx
                !$omp atomic
                z_force_sum(j) = z_force_sum(j) + fz
            end do
        end do
        !$omp end parallel do
        
    end subroutine coulomb_force_sub_full_pairs
    
    !> クーロン力計算の新旧実装を比較する検証サブルーチン
    !> 大きなカットオフで両者の力が一致するか確認する
    subroutine verify_coulomb_implementations
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force, &
            coulomb_cutoff, coulomb_softening, coulomb_shift_force, coulomb_use_cell
        use particle_data_mod
        use cell_system_mod
        implicit none
        
        real(8), allocatable :: fx_cell(:), fz_cell(:)
        real(8), allocatable :: fx_full(:), fz_full(:)
        real(8) :: max_diff_x, max_diff_z, total_diff
        real(8) :: saved_cutoff
        logical :: saved_use_cell, saved_shift
        integer :: i
        
        if (.not. enable_coulomb_force) then
            write(*,*) '[検証] クーロン力が無効のため検証スキップ'
            return
        end if
        
        if (num_particles <= 0) return
        
        allocate(fx_cell(num_particles), fz_cell(num_particles))
        allocate(fx_full(num_particles), fz_full(num_particles))
        
        ! 現在の設定を保存
        saved_cutoff = coulomb_cutoff
        saved_use_cell = coulomb_use_cell
        saved_shift = coulomb_shift_force
        
        ! 力をクリアして保存
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
        end do
        
        ! セル法で計算（大きなカットオフを設定）
        coulomb_cutoff = 1000.0d0  ! 十分大きいカットオフ
        coulomb_shift_force = .false.  ! 比較のためシフトドフォースを無効化
        call coulomb_force_cell_cutoff_sub()
        
        do i = 1, num_particles
            fx_cell(i) = x_force_sum(i)
            fz_cell(i) = z_force_sum(i)
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
        end do
        
        ! 全対全で計算
        coulomb_use_cell = .false.
        call coulomb_force_sub_full_pairs()
        
        do i = 1, num_particles
            fx_full(i) = x_force_sum(i)
            fz_full(i) = z_force_sum(i)
        end do
        
        ! 設定を復元
        coulomb_cutoff = saved_cutoff
        coulomb_use_cell = saved_use_cell
        coulomb_shift_force = saved_shift
        
        ! 差分を計算
        max_diff_x = 0.0d0
        max_diff_z = 0.0d0
        total_diff = 0.0d0
        
        do i = 1, num_particles
            max_diff_x = max(max_diff_x, abs(fx_cell(i) - fx_full(i)))
            max_diff_z = max(max_diff_z, abs(fz_cell(i) - fz_full(i)))
            total_diff = total_diff + abs(fx_cell(i) - fx_full(i)) + abs(fz_cell(i) - fz_full(i))
        end do
        
        write(*,*) '=== クーロン力実装検証 ==='
        write(*,*) '粒子数: ', num_particles
        write(*,*) 'セル法 vs 全対全 最大差(X): ', max_diff_x
        write(*,*) 'セル法 vs 全対全 最大差(Z): ', max_diff_z
        write(*,*) '総差分: ', total_diff
        
        if (max_diff_x < 1.0d-10 .and. max_diff_z < 1.0d-10) then
            write(*,*) '[OK] 新旧実装の力が一致しています'
        else
            write(*,*) '[警告] 新旧実装に差異があります（許容範囲を確認してください）'
        end if
        write(*,*) '=========================='
        
        ! 力をクリア（検証後に再計算されるため）
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
        end do
        
        deallocate(fx_cell, fz_cell, fx_full, fz_full)
        
    end subroutine verify_coulomb_implementations

    !> 蛙飛び法による粒子の位置と速度を更新するサブルーチン
    subroutine nposit_leapfrog_sub(judge_static, phase)
        ! phase: 0=初期化（v(0)→v(Δt/2)）, 1=位置更新, 2=速度更新
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use simulation_parameters_mod, only: time_step
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        
        integer, intent(in) :: phase
        integer, intent(out) :: judge_static
        
        integer :: i
        real(8) :: sum_abs_disp, avg_abs_disp
        real(8) :: grav, dt
        real(8) :: total_kinetic_energy, ke_trans, ke_rot
        
        dt = time_step
        grav = GRAVITY_ACCEL
        
        if (phase == 0) then
            ! 初回のみ: v(0) → v(Δt/2) への変換
            !$omp parallel do private(i)
            do i = 1, num_particles
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * (dt * 0.5d0)
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - grav) * (dt * 0.5d0)
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * (dt * 0.5d0)
            end do
            !$omp end parallel do
            judge_static = 0
            
        else if (phase == 1) then
            ! フェーズ1: 位置更新のみ
            sum_abs_disp = 0.0d0
            !$omp parallel do private(i) reduction(+:sum_abs_disp)
            do i = 1, num_particles
                ! 位置更新: x(t+Δt) = x(t) + v(t+Δt/2) * Δt
                x_disp_incr(i) = x_vel(i) * dt
                z_disp_incr(i) = z_vel(i) * dt
                rot_disp_incr(i) = rotation_vel(i) * dt
                
                x_coord(i) = x_coord(i) + x_disp_incr(i)
                z_coord(i) = z_coord(i) + z_disp_incr(i)
                rotation_angle(i) = rotation_angle(i) + rot_disp_incr(i)
                
                sum_abs_disp = sum_abs_disp + abs(x_disp_incr(i)) + abs(z_disp_incr(i))
            end do
            !$omp end parallel do
            
            ! 静止判定（運動エネルギーベース）
            if (num_particles > 0) then
                total_kinetic_energy = 0.0d0
                !$omp parallel do private(i, ke_trans, ke_rot) reduction(+:total_kinetic_energy)
                do i = 1, num_particles
                    ! 並進運動エネルギー: KE_trans = 0.5 * m * (vx² + vz²)
                    if (mass(i) > 1.0d-20) then
                        ke_trans = 0.5d0 * mass(i) * (x_vel(i)**2 + z_vel(i)**2)
                    else
                        ke_trans = 0.0d0
                    end if
                    
                    ! 回転運動エネルギー: KE_rot = 0.5 * I * ω²
                    if (moment_inertia(i) > 1.0d-20) then
                        ke_rot = 0.5d0 * moment_inertia(i) * rotation_vel(i)**2
                    else
                        ke_rot = 0.0d0
                    end if
                    
                    total_kinetic_energy = total_kinetic_energy + ke_trans + ke_rot
                end do
                !$omp end parallel do
                
                if (total_kinetic_energy < kinetic_energy_threshold) then
                    judge_static = 1
                else
                    judge_static = 0
                end if
            else
                judge_static = 1
            end if
            
        else if (phase == 2) then
            ! フェーズ2: 速度更新のみ（新しい位置での力を使用）
            !$omp parallel do private(i)
            do i = 1, num_particles
                ! 速度更新: v(t+3Δt/2) = v(t+Δt/2) + a(t+Δt) * Δt
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * dt
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - grav) * dt
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * dt
            end do
            !$omp end parallel do
            ! judge_static = 0 ! phase 2では静止判定をリセットしない
        end if
    end subroutine nposit_leapfrog_sub

    !> 粒子iと粒子/壁jとの間の実際の接触力（法線方向およびせん断方向）を計算するサブルーチン
    subroutine actf_sub(p_i, p_j, contact_slot_idx_for_pi, angle_sin, angle_cos, initial_overlap, en_coeff, et_coeff)
        ! p_i: 主となる粒子のインデックス
        ! p_j: 他方の粒子インデックス (<= num_particles の場合) または壁ID (> num_particles の場合)
        ! contact_slot_idx_for_pi: p_i の接触配列における p_j のスロット
        ! angle_sin, angle_cos: p_i から p_j の中心/接触点への法線ベクトル成分
        ! initial_overlap: 幾何学的な重なり量、接触していれば正
        ! en_coeff, et_coeff: 反発係数 (normal, tangential)
        use simulation_constants_mod, only: ni_max
        use simulation_parameters_mod, only: time_step, reference_overlap, rolling_friction_coeff_particle, rolling_friction_coeff_wall
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none

        integer, intent(in) :: p_i, p_j, contact_slot_idx_for_pi
        real(8), intent(in) :: angle_sin, angle_cos, initial_overlap
        real(8), intent(in) :: en_coeff, et_coeff ! 反発係数 (normal, tangential)

        real(8) :: ri_val, rj_val, effective_mass
        real(8) :: kn_normal_stiffness, ks_shear_stiffness ! 法線・せん断バネ定数 Kn, Ks
        real(8) :: damping_coeff_normal, damping_coeff_shear ! 法線・せん断粘性係数 ηn, ηs
        real(8) :: rel_disp_normal_incr, rel_disp_shear_incr ! 法線・せん断方向の相対変位増分 Δun, Δus
        real(8) :: damping_force_normal, damping_force_shear ! 法線・せん断方向の粘性抵抗力 dn, ds
        real(8) :: total_normal_force, total_shear_force     ! 全法線力 fn, 全せん断力 fs
        real(8) :: friction_coeff_current      ! 現在の摩擦係数 μ
        real(8) :: critical_time_step_check    ! 安定性チェック用の時間刻み (ddt)

        real(8) :: x_disp_incr_pi, z_disp_incr_pi, rot_disp_incr_pi ! 粒子iの変位増分
        real(8) :: x_disp_incr_pj, z_disp_incr_pj, rot_disp_incr_pj ! 粒子j(または壁=0)の変位増分
        real(8) :: mass_pi, mass_pj                                 ! 粒子i,jの質量
        real(8) :: r_eff ! 宣言をここに移動
        logical :: is_tracked_pair
        real(8) :: vx_rel, vz_rel, v_rel_n
        real(8) :: tau, alpha, wd, delta_theory
        real(8) :: s1, s2, denom, A_c, B_c
        real(8) :: vx_rel_current, vz_rel_current, v_rel_n_current, v_rel_theory
        ! 転がり摩擦関連
        real(8) :: rolling_coeff_current, rel_angular_vel, rolling_torque_mag
        real(8) :: rolling_torque_i, rolling_torque_j

        ri_val = radius(p_i)
        x_disp_incr_pi = x_disp_incr(p_i)
        z_disp_incr_pi = z_disp_incr(p_i)
        rot_disp_incr_pi = rot_disp_incr(p_i)
        mass_pi = mass(p_i)

        if (p_j <= num_particles) then ! 粒子間
            rj_val = radius(p_j)
            x_disp_incr_pj = x_disp_incr(p_j)
            z_disp_incr_pj = z_disp_incr(p_j)
            rot_disp_incr_pj = rot_disp_incr(p_j)
            mass_pj = mass(p_j)
            if (mass_pi + mass_pj > 1.0d-20) then
                 effective_mass = mass_pi * mass_pj / (mass_pi + mass_pj) ! 等価質量 m_eff = m1*m2/(m1+m2)
            else
                 effective_mass = mass_pi ! フォールバック
            end if
            friction_coeff_current = friction_coeff_particle
            r_eff = ri_val * rj_val / (ri_val + rj_val) ! 等価半径 Reff = ri*rj/(ri+rj)
            
            ! Hertz接触理論に基づく法線剛性(δ_refに基づく) (粒子間)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle / (1.0d0 - poisson_ratio_particle**2) * &
                                 sqrt(reference_overlap)
            ! write(*,*) '粒子間'
            ! write(*,*) 'kn_normal_stiffness=', kn_normal_stiffness

            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        else ! 粒子-壁接触
            rj_val = 0.0d0 ! 壁の半径は実質無限大、または変位にrjは使用しない
            x_disp_incr_pj = 0.0d0 ! 壁は動かないと仮定
            z_disp_incr_pj = 0.0d0
            rot_disp_incr_pj = 0.0d0
            effective_mass = mass_pi   ! 等価質量 m_eff = m1
            friction_coeff_current = friction_coeff_wall
            r_eff = ri_val ! 粒子-壁接触では粒子半径を使用
            
            ! Hertz接触理論に基づく法線剛性(δ_refに基づく) (粒子-壁)
            kn_normal_stiffness = (4.0d0/3.0d0) * sqrt(r_eff) * &
                                 young_modulus_particle * young_modulus_wall / &
                                 ((1.0d0-poisson_ratio_particle**2)*young_modulus_wall + &
                                  (1.0d0-poisson_ratio_wall**2)*young_modulus_particle) * &
                                 sqrt(reference_overlap)
            ! write(*,*) '粒子-壁'
            ! write(*,*) 'kn_normal_stiffness=', kn_normal_stiffness
            
            if (kn_normal_stiffness < 1.0d6) kn_normal_stiffness = 1.0d6 ! 最小値設定
        end if
        ks_shear_stiffness = kn_normal_stiffness * shear_to_normal_stiffness_ratio
        ! write(*,*) 'ks_shear_stiffness=', ks_shear_stiffness

        ! 完全弾性衝突 (e=1) のための粘性係数設定
        if (en_coeff >= 0.99d0) then
            ! 完全弾性衝突の場合、粘性をほぼゼロに設定
            damping_coeff_normal = 0.0d0
            damping_coeff_shear = 0.0d0
        else
            ! 粘性係数 (反発係数に基づいて計算)
            if (effective_mass > 0.0d0 .and. kn_normal_stiffness > 0.0d0 .and. en_coeff > 1.0d-6) then
                damping_coeff_normal = -2.0d0 * log(en_coeff) * sqrt(effective_mass * kn_normal_stiffness / (log(en_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_normal = 0.0d0
            end if
            if (effective_mass > 0.0d0 .and. ks_shear_stiffness > 0.0d0 .and. et_coeff > 1.0d-6) then
                damping_coeff_shear = -2.0d0 * log(et_coeff) * sqrt(effective_mass * ks_shear_stiffness / (log(et_coeff)**2 + PI_VAL**2))
            else
                damping_coeff_shear = 0.0d0
            end if
        end if

        ! 安定性基準のチェック (元のddt、レイリー時間刻みに関連)
        if (kn_normal_stiffness > 1.0d-12) then
           critical_time_step_check = 0.1d0 * sqrt(effective_mass / kn_normal_stiffness)
           if (critical_time_step_check < time_step .and. critical_time_step_check > 1.0d-12) then 
                write(*,*) '警告: 安定性基準違反 - 推奨時間刻み:', critical_time_step_check, '現在:', time_step
                write(*,*) '  Kn=', kn_normal_stiffness, ' M_eff=', effective_mass, ' 粒子:', p_i, p_j
           end if
        end if

        ! 相対変位増分 (現時間ステップ time_step における)
        ! angle成分は粒子p_iからp_jへの法線ベクトルを定義
        rel_disp_normal_incr = (x_disp_incr_pi - x_disp_incr_pj) * angle_cos + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_sin
        rel_disp_shear_incr = -(x_disp_incr_pi - x_disp_incr_pj) * angle_sin + &
                               (z_disp_incr_pi - z_disp_incr_pj) * angle_cos + &
                               (ri_val * rot_disp_incr_pi + rj_val * rot_disp_incr_pj)

        ! 弾性力成分の更新 (式3.7, 3.10)
        if (abs(normal_force_contact(p_i, contact_slot_idx_for_pi)) < 1.0d-8) then ! 新規接触
            ! 新規接触の場合、弾性力を重なり量に基づいて初期化
            normal_force_contact(p_i, contact_slot_idx_for_pi) = kn_normal_stiffness * initial_overlap ! 弾性項抜き出す部分
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0 ! せん断力は初期化時にゼロ
        else
            ! 既存接触の場合、増分で更新
            normal_force_contact(p_i, contact_slot_idx_for_pi) = normal_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                 kn_normal_stiffness * rel_disp_normal_incr 
            shear_force_contact(p_i, contact_slot_idx_for_pi) = shear_force_contact(p_i, contact_slot_idx_for_pi) + &
                                                                ks_shear_stiffness * rel_disp_shear_incr
        end if
        
        ! 粘性抵抗力成分の計算 (式3.6, 3.9)
        if (time_step > 1.0d-20) then
             damping_force_normal = damping_coeff_normal * rel_disp_normal_incr / time_step 
             damping_force_shear = damping_coeff_shear * rel_disp_shear_incr / time_step
        else
             damping_force_normal = 0.0d0
             damping_force_shear  = 0.0d0
        end if

        ! 引張力のチェック (粒子が引き離される場合) - 付着力は考慮しない
        if (normal_force_contact(p_i, contact_slot_idx_for_pi) < 0.0d0) then
            normal_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            shear_force_contact(p_i, contact_slot_idx_for_pi) = 0.0d0
            damping_force_normal = 0.0d0
            damping_force_shear = 0.0d0
            contact_partner_idx(p_i, contact_slot_idx_for_pi) = 0 ! 引張なら接触を切る
            return ! 引き離される場合は力なし
        end if

        ! クーロンの摩擦法則を適用 (式3.11)
        if (abs(shear_force_contact(p_i, contact_slot_idx_for_pi)) > &
            friction_coeff_current * normal_force_contact(p_i, contact_slot_idx_for_pi)) then
            shear_force_contact(p_i, contact_slot_idx_for_pi) = friction_coeff_current * &
                normal_force_contact(p_i, contact_slot_idx_for_pi) * &
                sign(1.0d0, shear_force_contact(p_i, contact_slot_idx_for_pi))
            damping_force_shear = 0.0d0 ! 滑りが発生している場合はせん断粘性なし
        end if

        ! 粘性を含む合計の力 (式3.8, 3.12)
        total_normal_force = normal_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_normal 
        total_shear_force = shear_force_contact(p_i, contact_slot_idx_for_pi) + damping_force_shear

        !---------------------------------------------------------------
        ! 転がり摩擦によるモーメント（非球形粒子の近似）
        !   |M_r| = μ_r * F_n * R_eff
        !   向き  = - sign(ω_rel)   （相対角速度に逆らう向き）
        !   - 粒子間:  μ_r = rolling_friction_coeff_particle
        !   - 粒子-壁: μ_r = rolling_friction_coeff_wall
        !---------------------------------------------------------------
        rolling_coeff_current = 0.0d0
        if (p_j <= num_particles) then
            rolling_coeff_current = rolling_friction_coeff_particle
        else
            rolling_coeff_current = rolling_friction_coeff_wall
        end if

        if (rolling_coeff_current > 0.0d0 .and. total_normal_force > 0.0d0) then
            if (p_j <= num_particles) then
                ! 粒子間: 相対角速度
                rel_angular_vel = rotation_vel(p_i) - rotation_vel(p_j)
            else
                ! 粒子-壁: 壁は回転しないと仮定
                rel_angular_vel = rotation_vel(p_i)
            end if

            if (abs(rel_angular_vel) > 1.0d-12) then
                rolling_torque_mag = rolling_coeff_current * total_normal_force * r_eff

                ! 粒子p_iへの転がり抵抗モーメント
                rolling_torque_i = -rolling_torque_mag * sign(1.0d0, rel_angular_vel)
                !$omp atomic
                moment_sum(p_i) = moment_sum(p_i) + rolling_torque_i

                ! 相手が粒子の場合は反作用モーメントも付与
                if (p_j <= num_particles) then
                    rolling_torque_j = -rolling_torque_i
                    !$omp atomic
                    moment_sum(p_j) = moment_sum(p_j) + rolling_torque_j
                end if
            end if
        end if

        ! 粒子p_iに力を適用 (式3.13)
        ! 法線力は中心を結ぶ線に沿って作用 (angle_cos, angle_sin で定義される iからjへの方向)
        ! せん断力はそれに垂直。
        !$omp atomic
        x_force_sum(p_i) = x_force_sum(p_i) - total_normal_force * angle_cos + total_shear_force * angle_sin
        !$omp atomic
        z_force_sum(p_i) = z_force_sum(p_i) - total_normal_force * angle_sin - total_shear_force * angle_cos
        !$omp atomic
        moment_sum(p_i) = moment_sum(p_i) - ri_val * total_shear_force

        ! 粒子p_jに反作用力を適用 (相手が粒子の場合)
        if (p_j <= num_particles .and. contact_slot_idx_for_pi <= 10) then ! 元の jk < 10 は粒子間接触のチェック
            !$omp atomic
            x_force_sum(p_j) = x_force_sum(p_j) + total_normal_force * angle_cos - total_shear_force * angle_sin
            !$omp atomic
            z_force_sum(p_j) = z_force_sum(p_j) + total_normal_force * angle_sin + total_shear_force * angle_cos
            !$omp atomic
            moment_sum(p_j) = moment_sum(p_j) - rj_val * total_shear_force ! p_jに対するせん断力は同じ大きさ、逆向きの回転効果
            
            ! デバッグ出力 (検証モードのみ)
            ! if (validation_mode .and. abs(total_normal_force) > 1.0d-6) then
            !     write(*,'(A,I2,A,I2,A,ES12.5,A,ES12.5)') 'Force p',p_i,' -> p',p_j,': Fn=',total_normal_force,', Fs=',total_shear_force
            !     write(*,'(A,ES12.5,A,ES12.5)') '  overlap=',initial_overlap,', dist=',sqrt((x_coord(p_i)-x_coord(p_j))**2 + (z_coord(p_i)-z_coord(p_j))**2)
            ! end if
        end if
        
        ! 現在のオーバーラップを記録（次ステップで接触開始検出に使用）
        previous_overlap(p_i, contact_slot_idx_for_pi) = initial_overlap
    end subroutine actf_sub

    !> グラフィック用データを出力するサブルーチン
    !> 出力形式: 標準CSV (analyze_repose_angle.py, animate_pem.py と互換)
    subroutine gfout_sub(iter_step, time_val, rmax_val)
        use simulation_constants_mod, only: nj_max, GRAVITY_ACCEL
        use simulation_parameters_mod, only: container_width, container_height, time_step, output_dir
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: iter_step    ! 現在のイテレーションステップ
        real(8), intent(in) :: time_val, rmax_val ! 現在時刻、最大粒子半径
        integer :: i,j
        real(8) :: dt, grav
        real(8), allocatable :: vx_out(:), vz_out(:), rotation_vel_out(:)
        character(len=512) :: file_path

        if (iter_step == 1) then
            file_path = trim(output_dir) // '/particles.csv'
            open(unit=10, file=trim(file_path), status='replace', action='write')
            ! CSVヘッダー行を最初に一度だけ出力
            write(10,'(A)') 'step,time,id,x,z,radius,vx,vz,omega,angle,mass,charge'
            
            file_path = trim(output_dir) // '/contacts.csv'
            open(unit=11, file=trim(file_path), status='replace', action='write')
        end if

        if (num_particles > 0) then
            ! 出力用補正速度の準備（蛙飛び法のみ 0.5*dt*加速度で補正）
            dt = time_step
            grav = GRAVITY_ACCEL

            allocate(vx_out(num_particles), vz_out(num_particles), rotation_vel_out(num_particles))
            ! 出力用補正速度の計算
            do i = 1, num_particles
                vx_out(i) = x_vel(i) - 0.5d0 * dt * (x_force_sum(i) / mass(i))
                vz_out(i) = z_vel(i) - 0.5d0 * dt * (z_force_sum(i) / mass(i) - grav)
                rotation_vel_out(i) = rotation_vel(i) - 0.5d0 * dt * (moment_sum(i) / moment_inertia(i))
            end do

            ! 標準CSV形式で出力 (step, time, id, x, z, radius, vx, vz, omega, angle, mass, charge)
            do i = 1, num_particles
                write(10,'(I8,",",ES14.7,",",I8,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7,",",ES14.7)') &
                    iter_step, time_val, i, x_coord(i), z_coord(i), radius(i), &
                    vx_out(i), vz_out(i), rotation_vel_out(i), rotation_angle(i), &
                    mass(i), charge(i)
            end do
            
            deallocate(vx_out, vz_out, rotation_vel_out)
        end if
        
        ! 接触力の出力 (オプション、graph21.dより)
        write(11,*) 'Time: ', time_val 
        write(11,*) 'Particle_ID,NumContacts,ShearForces...,NormalForces...,Partners...'
        if (num_particles > 0) then
            do i = 1, num_particles
                 write(11,'(I5,A,I5,A,40(ES10.3,A),40(ES10.3,A),40(I5,A))') &
                     i, ",", count(contact_partner_idx(i,1:nj_max) > 0), ",", &
                     (shear_force_contact(i,j), ",", j=1,nj_max), &
                     (normal_force_contact(i,j), ",", j=1,nj_max), &
                     (contact_partner_idx(i,j), ",", j=1,nj_max)
            end do
        end if
    end subroutine gfout_sub

    !> バックアップデータを出力するサブルーチン
    subroutine bfout_sub
        use simulation_constants_mod
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        integer :: i, j
        real(8) :: rmax_dummy_val ! 元のbfoutはrmaxを必要とするが、メインの呼び出しからは渡されない。
                                  ! リストアに不可欠でないか、粒子半径から取得すると仮定。
        character(len=512) :: file_path
        integer(profiler_tick_kind) :: prof_token

        call profiler_start('bfout_sub', prof_token)
        
        if (num_particles > 0) then
           rmax_dummy_val = maxval(radius(1:num_particles))
        else
           rmax_dummy_val = 0.0d0
        end if

        file_path = trim(output_dir) // '/backl.d'
        open(unit=13, file=trim(file_path), status='replace', action='write')

        write(13,*) num_particles, cells_x_dir, cells_z_dir, particle_gen_layers
        write(13,*) rmax_dummy_val, 0.0d0, container_width, container_height, cell_size, time_step ! current_timeではなく初期t=0を保存すると仮定
        write(13,*) particle_density, friction_coeff_particle, friction_coeff_wall, GRAVITY_ACCEL, PI_VAL
        write(13,*) young_modulus_particle, young_modulus_wall, poisson_ratio_particle, poisson_ratio_wall, shear_to_normal_stiffness_ratio
        
        if (num_particles > 0) then
            write(13,*) (mass(i), moment_inertia(i), i=1,num_particles)
            write(13,*) (x_coord(i), z_coord(i), radius(i), charge(i), i=1,num_particles)
            write(13,*) (x_disp_incr(i), z_disp_incr(i), rot_disp_incr(i), i=1,num_particles) ! u,v,f (dpm)
            write(13,*) (x_vel(i), z_vel(i), rotation_vel(i), i=1,num_particles)            ! u0,v0,f0
            do i = 1, num_particles
                write(13,*) (shear_force_contact(i,j), normal_force_contact(i,j), j=1,nj_max)
                write(13,*) (contact_partner_idx(i,j), j=1,nj_max)
            end do
        end if
        close(13)
        call profiler_stop('bfout_sub', prof_token)
    end subroutine bfout_sub

    !> 充填状態の粒子データを保存するサブルーチン
    subroutine save_filled_particles_sub
        use simulation_parameters_mod, only: output_dir, use_explicit_positions
        use particle_data_mod
        use cell_system_mod, only: num_particles
        use profiling_mod, only: profiler_start, profiler_stop, profiler_tick_kind
        implicit none
        integer :: i, unit_num, ios
        character(len=512) :: file_path
        integer(profiler_tick_kind) :: prof_token
        
        call profiler_start('save_filled_particles_sub', prof_token)
        
        ! 既に filled_particles.dat 等の明示配置から開始している場合は、
        ! sweep中に inputs/filled_particles.dat を上書きして初期条件が混ざるのを避ける。
        if (use_explicit_positions) then
            file_path = trim(output_dir) // '/filled_particles.dat'
        else
            file_path = 'inputs/filled_particles.dat'
        end if
        open(newunit=unit_num, file=trim(file_path), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'エラー: 充填状態ファイルを開けません: ', trim(file_path)
            call profiler_stop('save_filled_particles_sub', prof_token)
            return
        end if
        
        ! ヘッダーコメントを書き込み
        write(unit_num, '(A)') '# 充填・堆積完了時の粒子状態'
        write(unit_num, '(A)') '# フォーマット: x z radius [charge]'
        write(unit_num, '(A)') '# x: 粒子中心のx座標 [m]'
        write(unit_num, '(A)') '# z: 粒子中心のz座標 [m]'
        write(unit_num, '(A)') '# radius: 粒子半径 [m]'
        write(unit_num, '(A)') '# charge: 粒子電荷 [C] (オプション)'
        
        ! 粒子データを書き込み
        if (num_particles > 0) then
            do i = 1, num_particles
                write(unit_num, '(ES14.7,1X,ES14.7,1X,ES14.7,1X,ES14.7)') &
                    x_coord(i), z_coord(i), radius(i), charge(i)
            end do
        end if
        
        close(unit_num)
        write(*,*) '充填状態を保存しました: ', trim(file_path), ' (粒子数: ', num_particles, ')'
        
        call profiler_stop('save_filled_particles_sub', prof_token)
    end subroutine save_filled_particles_sub

    !> 擬似乱数を生成するサブルーチン
    subroutine custom_random(seed_io, random_val_out)
        implicit none
        integer, intent(inout) :: seed_io          ! ジェネレータの現在の状態を保持
        real(8), intent(out)   :: random_val_out   ! [0,1) の一様乱数

        ! 元のコードの定数とロジックを可能な限り再現
        seed_io = seed_io * 65539 
        if (seed_io < 0) then
             seed_io = (seed_io + 2147483647) + 1 
        end if
        random_val_out = dble(seed_io) * 0.4656613d-9 ! 元の正規化定数 (1.0 / 2147483648.0)

    end subroutine custom_random

    !> 一様分布から乱数を生成するサブルーチン
    subroutine generate_uniform_random(seed_io, min_val, max_val, result_out)
        implicit none
        integer, intent(inout) :: seed_io
        real(8), intent(in) :: min_val, max_val
        real(8), intent(out) :: result_out
        real(8) :: u
        
        call custom_random(seed_io, u)
        result_out = min_val + u * (max_val - min_val)
    end subroutine generate_uniform_random

    !> 正規分布から乱数を生成するサブルーチン (Box-Muller法)
    subroutine generate_normal_random(seed_io, mean_val, std_val, result_out)
        use simulation_constants_mod, only: PI_VAL
        implicit none
        integer, intent(inout) :: seed_io
        real(8), intent(in) :: mean_val, std_val
        real(8), intent(out) :: result_out
        real(8) :: u1, u2, z
        
        ! Box-Muller変換
        call custom_random(seed_io, u1)
        call custom_random(seed_io, u2)
        
        ! u1が0に非常に近い場合を回避
        if (u1 < 1.0d-10) u1 = 1.0d-10
        
        z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI_VAL * u2)
        result_out = mean_val + std_val * z
    end subroutine generate_normal_random

    !> 対数正規分布から乱数を生成するサブルーチン
    !! mean_val: 対数正規分布の平均（元スケール）
    !! std_val: 対数正規分布の標準偏差（元スケール）
    subroutine generate_lognormal_random(seed_io, mean_val, std_val, result_out)
        use simulation_constants_mod, only: PI_VAL
        implicit none
        integer, intent(inout) :: seed_io
        real(8), intent(in) :: mean_val, std_val
        real(8), intent(out) :: result_out
        real(8) :: u1, u2, z, mu_ln, sigma_ln, cv_sq
        
        ! 変動係数の二乗
        cv_sq = (std_val / mean_val) ** 2
        
        ! 対数正規分布のパラメータに変換
        ! μ_ln = ln(mean) - 0.5 * ln(1 + (std/mean)^2)
        ! σ_ln = sqrt(ln(1 + (std/mean)^2))
        sigma_ln = sqrt(log(1.0d0 + cv_sq))
        mu_ln = log(mean_val) - 0.5d0 * sigma_ln * sigma_ln
        
        ! Box-Muller変換で標準正規分布を生成
        call custom_random(seed_io, u1)
        call custom_random(seed_io, u2)
        if (u1 < 1.0d-10) u1 = 1.0d-10
        
        z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI_VAL * u2)
        
        ! 対数正規分布に変換
        result_out = exp(mu_ln + sigma_ln * z)
    end subroutine generate_lognormal_random

    !> 指数分布から乱数を生成するサブルーチン（逆関数法）
    !! mean_val: 指数分布の平均（= 1/λ）
    subroutine generate_exponential_random(seed_io, mean_val, result_out)
        implicit none
        integer, intent(inout) :: seed_io
        real(8), intent(in) :: mean_val
        real(8), intent(out) :: result_out
        real(8) :: u
        
        call custom_random(seed_io, u)
        
        ! u が 0 に非常に近い場合を回避
        if (u < 1.0d-10) u = 1.0d-10
        
        ! 逆関数法: X = -mean * ln(U)
        result_out = -mean_val * log(u)
    end subroutine generate_exponential_random

    !> 二峰性分布から乱数を生成するサブルーチン
    !! 確率ratioで第1ピーク(mean1, std1)、それ以外で第2ピーク(mean2, std2)
    subroutine generate_bimodal_random(seed_io, mean1, std1, mean2, std2, ratio, result_out)
        use simulation_constants_mod, only: PI_VAL
        implicit none
        integer, intent(inout) :: seed_io
        real(8), intent(in) :: mean1, std1, mean2, std2, ratio
        real(8), intent(out) :: result_out
        real(8) :: u, u1, u2, z, selected_mean, selected_std
        
        ! どちらのピークを選ぶか決定
        call custom_random(seed_io, u)
        
        if (u < ratio) then
            selected_mean = mean1
            selected_std = std1
        else
            selected_mean = mean2
            selected_std = std2
        end if
        
        ! Box-Muller変換で正規分布を生成
        call custom_random(seed_io, u1)
        call custom_random(seed_io, u2)
        if (u1 < 1.0d-10) u1 = 1.0d-10
        
        z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI_VAL * u2)
        result_out = selected_mean + selected_std * z
    end subroutine generate_bimodal_random

    !> 粒子の衝突チェックサブルーチン
    !! 新しい粒子が既存粒子と重なっているかをチェック
    !! x_new, z_new: 新粒子の座標
    !! r_new: 新粒子の半径
    !! n_existing: 既存粒子数
    !! has_collision: 出力 .true. = 衝突あり, .false. = 衝突なし
    subroutine check_particle_collision_sub(x_new, z_new, r_new, n_existing, has_collision)
        use particle_data_mod, only: x_coord, z_coord, radius
        implicit none
        real(8), intent(in) :: x_new, z_new, r_new
        integer, intent(in) :: n_existing
        logical, intent(out) :: has_collision
        integer :: i
        real(8) :: dx, dz, dist_sq, min_dist_sq
        
        has_collision = .false.
        
        do i = 1, n_existing
            dx = x_new - x_coord(i)
            dz = z_new - z_coord(i)
            dist_sq = dx * dx + dz * dz
            min_dist_sq = (r_new + radius(i)) ** 2
            
            if (dist_sq < min_dist_sq) then
                has_collision = .true.
                return
            end if
        end do
    end subroutine check_particle_collision_sub

end program two_dimensional_dem