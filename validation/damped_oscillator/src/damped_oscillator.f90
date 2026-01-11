program damped_oscillator
    implicit none

    ! 固定パラメータ
    real(8), parameter :: mass_m      = 1.04720d-2
    real(8), parameter :: stiffness_k = 1.0d6
    real(8) :: damping_c   
    real(8), parameter :: restitution_coeff_e = 0.8d0
    real(8), parameter :: PI_VAL = 3.14159265358979323846d0

    ! 斜面検証用パラメータ
    ! real(8), parameter :: GRAVITY_ACCEL = 9.80665d0
    ! real(8), parameter :: ANGLE = 30.0d0
    ! real(8), parameter :: ANGLE_RAD = ANGLE * PI_VAL / 180.0d0
    ! real(8), parameter :: SLOPE_FORCE = cos(ANGLE_RAD) * GRAVITY_ACCEL
    
    ! 実行時設定（既定値）
    real(8) :: x0, v0, dt, t_end
    integer :: output_stride

    ! 状態変数
    real(8) :: x, v, t
    real(8) :: v_out
    integer :: n_steps, i, ios

    ! 周波数・減衰比（理論解用）
    real(8) :: omega_n, zeta, omega_d

    ! 理論解
    real(8) :: x_theory, v_theory

    ! 誤差評価
    real(8) :: max_abs_err_x, max_abs_err_v

    ! 出力
    integer :: unit_csv

    ! 既定値（ユーザ指定がなければこれを使用）
    x0 = 0.0d0
    v0 = -0.1d0
    dt = 1.0d-6
    t_end = 0.02
    output_stride = 1

    damping_c = -2.0d0 * log(restitution_coeff_e) * sqrt(mass_m * stiffness_k / (log(restitution_coeff_e)**2 + PI_VAL**2))
            
    call read_cli_args(x0, v0, dt, t_end)

    ! 物理量の導出
    omega_n = sqrt(stiffness_k / mass_m)
    if (mass_m > 0.0d0 .and. stiffness_k > 0.0d0) then
        zeta = damping_c / (2.0d0 * sqrt(stiffness_k * mass_m))
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
    write(*,*) 'マスバネダンパ系（減衰振動）'
    write(*,*) '数値積分法: 蛙飛び法 (Leapfrog)'
    write(*,'(A,ES14.6)') 'm = ', mass_m
    write(*,'(A,ES14.6)') 'k = ', stiffness_k
    write(*,'(A,ES14.6)') 'c = ', damping_c
    write(*,'(A,ES14.6)') 'zeta = ', zeta
    write(*,'(A,ES14.6)') 'omega_n = ', omega_n
    write(*,'(A,ES14.6)') 'omega_d = ', omega_d
    write(*,'(A,ES14.6)') 'dt = ', dt
    write(*,'(A,ES14.6)') 't_end = ', t_end
    write(*,*) '================================='

    ! 初期化
    x = x0
    v = v0
    t = 0.0d0

    ! ステップ数
    if (dt > 0.0d0) then
        n_steps = int(t_end / dt)
    else
        n_steps = 0
    end if

    ! 出力ディレクトリ作成（存在しない場合）
    call ensure_data_dir()

    unit_csv = 21
    open(unit=unit_csv, file='data/damped_oscillator.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'CSVファイルを開けません: data/damped_oscillator.csv'
        stop 1
    end if
    write(unit_csv,'(A)') 'time,x_num,v_num,x_theory,v_theory'

    max_abs_err_x = 0.0d0
    max_abs_err_v = 0.0d0

    ! 蛙飛び法の初期化: v(0) → v(dt/2)
    v = v + 0.5d0 * dt * (-(damping_c/mass_m) * v - (stiffness_k/mass_m) * x)

    do i = 0, n_steps
        t = dble(i) * dt

        call theory_solution(t, x0, v0, omega_n, zeta, omega_d, x_theory, v_theory)
        
        ! 出力用の速度(蛙飛び法の半周期分のずれを補正)
        v_out = v + 0.5d0 * (-dt) * (-(damping_c/mass_m) * v - (stiffness_k/mass_m) * x)

        if (mod(i, output_stride) == 0) then
            write(unit_csv,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') t, ',', x, ',', v_out, ',', x_theory, ',', v_theory
        end if

        max_abs_err_x = max(max_abs_err_x, abs(x - x_theory))
        max_abs_err_v = max(max_abs_err_v, abs(v - v_theory))

        ! 蛙飛び法 (Leapfrog method)
        ! x(t+dt) = x(t) + v(t+dt/2) * dt
        x = x + v * dt
        ! v(t+3dt/2) = v(t+dt/2) + a(t+dt) * dt
        v = v + dt * (-(damping_c/mass_m) * v - (stiffness_k/mass_m) * x ) !- SLOPE_FORCE) 
    end do

    close(unit_csv)

    write(*,*) '最大絶対誤差: '
    write(*,'(A,ES14.6)') '  |x_num - x_theory|_max = ', max_abs_err_x
    write(*,'(A,ES14.6)') '  |v_num - v_theory|_max = ', max_abs_err_v
    write(*,*) 'CSV: data/damped_oscillator.csv に出力しました。'

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

    subroutine theory_solution(t, x0, v0, omega_n, zeta, omega_d, x_out, v_out)
        implicit none
        real(8), intent(in) :: t, x0, v0, omega_n, zeta, omega_d
        real(8), intent(out) :: x_out, v_out
        real(8) :: expfac, A, B, s, c

        if (zeta < 1.0d0 - 1.0d-12) then
            A = x0
            if (omega_d > 0.0d0) then
                B = (v0 + zeta*omega_n*x0) / omega_d
            else
                B = 0.0d0
            end if
            expfac = exp(-zeta*omega_n*t)
            s = sin(omega_d*t)
            c = cos(omega_d*t)
            x_out = expfac * (A*c + B*s)
            v_out = expfac * ( -zeta*omega_n*(A*c + B*s) + (-A*omega_d*s + B*omega_d*c) )
        else
            ! 臨界/過減衰の簡易式（ここでは臨界近似）
            expfac = exp(-omega_n*t)
            x_out = expfac * (x0 + (v0 + omega_n*x0)*t)
            v_out = expfac * (v0 - omega_n*(v0 + omega_n*x0)*t)
        end if
    end subroutine theory_solution

    subroutine ensure_data_dir()
        implicit none
        integer :: stat
        ! Fortran 2008: 実行コマンドでディレクトリを作成（存在してもOK）
        call execute_command_line('mkdir -p data', wait=.true., exitstat=stat)
    end subroutine ensure_data_dir

end program damped_oscillator


