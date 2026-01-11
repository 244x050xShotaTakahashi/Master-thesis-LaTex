program particle_1dcollision
    implicit none

    ! 固定パラメータ（与えられた値）
    real(8), parameter :: stiffness_k = 1.0d0
    real(8), parameter :: restitution_coeff_e = 0.25d0
    real(8), parameter :: PI_VAL = 3.14159265358979323846d0
    
    ! 粒子パラメータ（配列化）
    integer, parameter :: num_particles = 2
    real(8), dimension(num_particles) :: x, v, mass_p, radius_p
    real(8), dimension(num_particles) :: v_previous
    
    ! 等価パラメータ
    real(8) :: effective_mass, effective_stiffness, effective_damping

    ! 実行時設定（既定値）
    real(8) :: x0, v0, dt, t_end
    integer :: output_stride

    ! 状態変数
    real(8) :: t
    integer :: n_steps, i, ios

    ! 接触力関連変数
    real(8) :: overlap, contact_force, previous_contact_force

    ! 周波数・減衰比（理論解用）
    real(8) :: omega_n, zeta, omega_d

    ! 理論解
    real(8) :: overlap_theory, v_theory

    ! オーバーラップ理論比較用変数
    logical :: ov_active
    real(8) :: ov_delta0, ov_deltadot0, ov_taccum
    real(8) :: ov_omega_n, ov_zeta, ov_omega_d
    real(8) :: ov_m, ov_k, ov_c
    real(8) :: ov_rmse_accum, ov_max_abs_err
    integer :: ov_num_samples

    ! 出力
    integer :: unit_csv

    ! 既定値（ユーザ指定がなければこれを使用）
    dt = 1.0d-6
    t_end = 5.0d0  
    output_stride = 1
    
    ! 粒子パラメータの設定
    mass_p(1) = 1.0d0      ! 粒子1の質量（移動粒子）
    mass_p(2) = 1.0d0      ! 粒子2の質量（静止粒子）
    radius_p(1) = 1.0d0    ! 粒子1の半径
    radius_p(2) = 1.0d0    ! 粒子2の半径

    ! 初期条件設定
    ! 粒子間距離をほぼ接触直前（gap = 0.001）に設定
    x(1) = 3.0d0                              ! 粒子1：左側から右に移動
    v(1) = 1.0d0                              ! 粒子1の速度
    x(2) = x(1) + radius_p(1) + radius_p(2) + 0.001d0  ! 粒子2：接触直前の位置で静止
    v(2) = 0.0d0                              ! 粒子2の速度（静止）

    ! コマンドライン引数の読み込み（x0, v0は使用しないが互換性のため残す）
    call read_cli_args(x0, v0, dt, t_end)

    ! 等価質量・等価剛性の計算
    effective_mass = (mass_p(1) * mass_p(2)) / (mass_p(1) + mass_p(2))
    effective_stiffness = stiffness_k
    
    ! 等価減衰係数をe, m_eff, kから計算
    effective_damping = -2.0d0 * log(restitution_coeff_e) * sqrt(effective_mass * effective_stiffness / &
                        (log(restitution_coeff_e)**2 + PI_VAL**2))

    ! 物理量の導出（等価パラメータを使用）
    omega_n = sqrt(effective_stiffness / effective_mass)
    if (effective_mass > 0.0d0 .and. effective_stiffness > 0.0d0) then
        zeta = effective_damping / (2.0d0 * sqrt(effective_stiffness * effective_mass))
    else
        zeta = 0.0d0
    end if
    if (zeta < 1.0d0) then
        omega_d = omega_n * sqrt(max(0.0d0, 1.0d0 - zeta*zeta))
    else
        omega_d = 0.0d0
    end if

    ! 見やすいパラメータの標準出力
    write(*,*) '================================='
    write(*,*) '粒子間衝突検証（1次元DEM）'
    write(*,*) '数値積分法: 蛙飛び法 (Leapfrog)'
    write(*,'(A,ES14.6)') 'm1 = ', mass_p(1)
    write(*,'(A,ES14.6)') 'm2 = ', mass_p(2)
    write(*,'(A,ES14.6)') 'm_eff = ', effective_mass
    write(*,'(A,ES14.6)') 'k = ', effective_stiffness
    write(*,'(A,ES14.6)') 'c = ', effective_damping
    write(*,'(A,ES14.6)') 'e = ', restitution_coeff_e
    write(*,'(A,ES14.6)') 'zeta = ', zeta
    write(*,'(A,ES14.6)') 'omega_n = ', omega_n
    write(*,'(A,ES14.6)') 'omega_d = ', omega_d
    write(*,'(A,ES14.6)') 'radius1 = ', radius_p(1)
    write(*,'(A,ES14.6)') 'radius2 = ', radius_p(2)
    write(*,'(A,ES14.6)') 'dt = ', dt
    write(*,'(A,ES14.6)') 't_end = ', t_end
    write(*,*) '================================='

    ! 初期化
    t = 0.0d0
    contact_force = 0.0d0
    previous_contact_force = 0.0d0
    
    ! オーバーラップ理論比較の初期化
    ov_active = .false.
    ov_rmse_accum = 0.0d0
    ov_max_abs_err = 0.0d0
    ov_num_samples = 0

    ! ステップ数
    if (dt > 0.0d0) then
        n_steps = int(t_end / dt)
    else
        n_steps = 0
    end if

    ! 出力ディレクトリ作成（存在しない場合）
    call ensure_data_dir()

    unit_csv = 21
    open(unit=unit_csv, file='data/particle_collision_validation.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'CSVファイルを開けません: data/particle_collision_validation.csv'
        stop 1
    end if
    write(unit_csv,'(A)') 'time,overlap_numeric,overlap_theory,v_rel_numeric,v_rel_theory,mode'

    ! 蛙飛び法の初期化: v(0) → v(dt/2)
    ! 初期状態では接触していないので、接触力はゼロ
    call pcont_sub(v(1), v(2), overlap, contact_force, 0)
    v(1) = v(1) + 0.5d0 * dt * (-contact_force / mass_p(1))
    v(2) = v(2) + 0.5d0 * dt * (contact_force / mass_p(2))  ! 反作用
    previous_contact_force = contact_force

    do i = 0, n_steps
        t = dble(i) * dt

        ! 蛙飛び法 (Leapfrog method)
        ! 位置更新: x(t+dt) = x(t) + v(t+dt/2) * dt
        x(1) = x(1) + v(1) * dt
        x(2) = x(2) + v(2) * dt
        
        ! 時刻tにおける速度を計算（前のステップの力を使用して半ステップ戻す）
        v_previous(1) = v(1) + 0.5d0 * dt * (-previous_contact_force / mass_p(1))
        v_previous(2) = v(2) + 0.5d0 * dt * (previous_contact_force / mass_p(2))
        
        ! 接触判定・力計算（時刻tの速度を使用）
        call pcont_sub(v_previous(1), v_previous(2), overlap, contact_force, i)
        
        ! 速度更新: v(t+3dt/2) = v(t+dt/2) + a(t+dt) * dt
        v(1) = v(1) + dt * (-contact_force / mass_p(1))
        v(2) = v(2) + dt * (contact_force / mass_p(2))  ! 反作用
        previous_contact_force = contact_force
    end do

    close(unit_csv)

    ! 最終統計の出力
    if (ov_num_samples > 0) then
        write(*,*) ''
        write(*,*) '================================='
        write(*,*) 'オーバーラップ理論比較結果'
        write(*,*) '================================='
        write(*,'(A,ES14.6)') 'RMSE: ', sqrt(ov_rmse_accum / dble(ov_num_samples))
        write(*,'(A,ES14.6)') 'MaxAbsErr: ', ov_max_abs_err
        write(*,'(A,I8)') 'サンプル数: ', ov_num_samples
        write(*,*) '================================='
    end if
    write(*,*) 'CSV: data/particle_collision_validation.csv に出力しました。'

contains

    subroutine read_cli_args(x0, v0, dt, t_end)
        implicit none
        real(8), intent(inout) :: x0, v0, dt, t_end
        integer :: argc
        character(len=64) :: arg

        argc = command_argument_count()
        if (argc >= 1) then
            call get_command_argument(1, arg); read(arg,*,err=10) x0
        end if
        if (argc >= 2) then
            call get_command_argument(2, arg); read(arg,*,err=10) v0
        end if
        if (argc >= 3) then
            call get_command_argument(3, arg); read(arg,*,err=10) dt
        end if
        if (argc >= 4) then
            call get_command_argument(4, arg); read(arg,*,err=10) t_end
        end if
10      continue
    end subroutine read_cli_args
  
    !> 粒子-壁間の接触判定とオーバーラップ計算、接触力計算を行うサブルーチン(現状は右方向のみ)
    subroutine wcont_sub(x_particle, v_particle, overlap_out, force_out)
        implicit none
        real(8), intent(in) :: x_particle, v_particle
        real(8), intent(out) :: overlap_out, force_out
        
        ! 壁パラメータ（このサブルーチンでのみ使用、実際の計算では使わない）
        real(8), parameter :: radius = 1.0d0
        real(8), parameter :: wall_x = 10.0d0
        
        ! オーバーラップ計算: δ = (x_particle + radius) - wall_x
        overlap_out = (x_particle + radius) - wall_x
        
        if (overlap_out > 0.0d0) then
            ! 接触している場合、接触力を計算
            call actf_sub(overlap_out, v_particle, force_out, 0)
        else
            ! 接触していない場合
            force_out = 0.0d0
            overlap_out = 0.0d0
            ! 接触終了検出
            if (ov_active) then
                ov_active = .false.
                write(*,*) '接触終了検出: t = ', t
            end if
        end if
    end subroutine wcont_sub
    
    !> 粒子間の接触判定とオーバーラップ計算、接触力計算を行うサブルーチン
    subroutine pcont_sub(v1_at_t, v2_at_t, overlap_out, force_out, step_num)
        implicit none
        real(8), intent(in) :: v1_at_t, v2_at_t  ! 時刻tにおける両粒子の速度
        integer, intent(in) :: step_num           ! ステップ番号（出力制御用）
        real(8), intent(out) :: overlap_out, force_out
        real(8) :: distance, sum_radius, v_rel
        
        ! 粒子間距離とオーバーラップ計算
        distance = abs(x(2) - x(1))
        sum_radius = radius_p(1) + radius_p(2)
        overlap_out = sum_radius - distance
        
        if (overlap_out > 0.0d0) then
            ! 接触している場合、時刻tでの相対速度を計算
            v_rel = v1_at_t - v2_at_t  ! 粒子2から見た粒子1の相対速度
            
            ! 接触力を計算
            call actf_sub(overlap_out, v_rel, force_out, step_num)
        else
            ! 接触していない場合
            force_out = 0.0d0
            overlap_out = 0.0d0
            ! 接触終了検出
            if (ov_active) then
                ov_active = .false.
                write(*,*) '接触終了検出: t = ', t
            end if
        end if
    end subroutine pcont_sub
    
    !> 接触力を計算するサブルーチン
    subroutine actf_sub(overlap_val, v_rel, force_out, step_num)
        implicit none
        real(8), intent(in) :: overlap_val, v_rel  ! v_relを相対速度として受け取る
        integer, intent(in) :: step_num           ! ステップ番号（出力制御用）
        real(8), intent(out) :: force_out
        real(8) :: elastic_force, damping_force
        real(8) :: tau, alpha, wd, delta_theory, v_rel_theory
        real(8) :: A_c, B_c, s1, s2, denom
        character(len=20) :: damping_mode
        
        ! 接触開始検出
        if (.not. ov_active) then
            ! 接触開始時の初期化
            ov_active = .true.
            ov_delta0 = 0.0d0  ! 接触開始時のオーバーラップ
            ov_deltadot0 = v_rel  ! 接触開始時の相対速度
            
            ! 等価パラメータを設定
            ov_m = effective_mass
            ov_k = effective_stiffness
            ov_c = effective_damping
            ov_omega_n = 0.0d0
            ov_zeta = 0.0d0
            ov_omega_d = 0.0d0
            
            if (ov_m > 0.0d0 .and. ov_k > 0.0d0) then
                ov_omega_n = sqrt(ov_k / ov_m)
                if (ov_c > 0.0d0) ov_zeta = ov_c / (2.0d0 * sqrt(ov_k * ov_m))
                if (ov_zeta < 1.0d0) ov_omega_d = ov_omega_n * sqrt(max(0.0d0, 1.0d0 - ov_zeta*ov_zeta))
            end if
            
            ov_taccum = 0.0d0
            ov_rmse_accum = 0.0d0
            ov_max_abs_err = 0.0d0
            ov_num_samples = 0
            
            write(*,*) '接触開始検出: t = ', t
            write(*,'(A,ES14.6)') '  m = ', ov_m
            write(*,'(A,ES14.6)') '  k = ', ov_k
            write(*,'(A,ES14.6)') '  c = ', ov_c
            write(*,'(A,ES14.6)') '  zeta = ', ov_zeta
            write(*,'(A,ES14.6)') '  omega_n = ', ov_omega_n
        end if
        
        ! 接触力の計算
        ! 弾性力（等価剛性を使用）
        elastic_force = effective_stiffness * overlap_val
        
        ! 粘性力（等価減衰係数と相対速度を使用）
        damping_force = effective_damping * v_rel
        
        ! 合力（粒子1に作用する力、左向きが正）
        force_out = elastic_force + damping_force
        
        if (ov_active) then
            ! 理論解との比較
            tau = ov_taccum
            
            ! 減衰モードの判定と理論解の計算
            if (ov_c <= 1.0d-16) then
                ! 減衰なし
                if (ov_omega_n > 0.0d0) then
                    delta_theory = ov_delta0 * cos(ov_omega_n * tau) + &
                                  (ov_deltadot0 / ov_omega_n) * sin(ov_omega_n * tau)
                    v_rel_theory = -ov_omega_n * ov_delta0 * sin(ov_omega_n * tau) + &
                                   ov_deltadot0 * cos(ov_omega_n * tau)
                else
                    delta_theory = ov_delta0
                    v_rel_theory = 0.0d0
                end if
                damping_mode = 'undamped'
            else if (ov_zeta < 1.0d0 - 1.0d-12) then
                ! 減衰振動
                alpha = ov_zeta * ov_omega_n
                wd = ov_omega_d
                if (wd > 0.0d0) then
                    delta_theory = exp(-alpha * tau) * &
                                  (ov_delta0 * cos(wd * tau) + &
                                   ((ov_deltadot0 + alpha * ov_delta0) / wd) * sin(wd * tau))
                    v_rel_theory = exp(-alpha * tau) * &
                                  (-alpha * ov_delta0 * cos(wd * tau) - &
                                   ov_delta0 * wd * sin(wd * tau) + &
                                   (-alpha * (ov_deltadot0 + alpha * ov_delta0) / wd) * sin(wd * tau) + &
                                   (ov_deltadot0 + alpha * ov_delta0) * cos(wd * tau))
                else
                    delta_theory = 0.0d0
                    v_rel_theory = 0.0d0
                end if
                damping_mode = 'underdamped'
            else if (abs(ov_zeta - 1.0d0) <= 1.0d-12) then
                ! 臨界減衰
                alpha = ov_omega_n
                delta_theory = exp(-alpha * tau) * &
                              (ov_delta0 + (ov_deltadot0 + alpha * ov_delta0) * tau)
                v_rel_theory = exp(-alpha * tau) * &
                              (-alpha * ov_delta0 - alpha * (ov_deltadot0 + alpha * ov_delta0) * tau + &
                               (ov_deltadot0 + alpha * ov_delta0))
                damping_mode = 'critical'
            else
                ! 過減衰
                s1 = -ov_omega_n * (ov_zeta - sqrt(ov_zeta*ov_zeta - 1.0d0))
                s2 = -ov_omega_n * (ov_zeta + sqrt(ov_zeta*ov_zeta - 1.0d0))
                denom = s1 - s2
                if (abs(denom) > 1.0d-30) then
                    A_c = (ov_deltadot0 - s2 * ov_delta0) / denom
                    B_c = (s1 * ov_delta0 - ov_deltadot0) / denom
                    delta_theory = A_c * exp(s1 * tau) + B_c * exp(s2 * tau)
                    v_rel_theory = A_c * s1 * exp(s1 * tau) + B_c * s2 * exp(s2 * tau)
                else
                    delta_theory = ov_delta0
                    v_rel_theory = 0.0d0
                end if
                damping_mode = 'overdamped'
            end if
            
            ! CSV出力と誤差計算（output_strideごとに間引く）
            if (mod(step_num, output_stride) == 0) then
                write(unit_csv,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,A)') &
                    t, ',', overlap_val, ',', delta_theory, ',', v_rel, ',', v_rel_theory, ',', trim(damping_mode)
                
                ! 誤差の蓄積（出力したデータのみ）
                ov_rmse_accum = ov_rmse_accum + (overlap_val - delta_theory)**2
                ov_max_abs_err = max(ov_max_abs_err, abs(overlap_val - delta_theory))
                ov_num_samples = ov_num_samples + 1
            end if
            
            ov_taccum = ov_taccum + dt
        end if
    end subroutine actf_sub

    subroutine ensure_data_dir()
        implicit none
        integer :: stat
        ! Fortran 2008: 実行コマンドでディレクトリを作成（存在してもOK）
        call execute_command_line('mkdir -p data', wait=.true., exitstat=stat)
    end subroutine ensure_data_dir

end program particle_1dcollision


