! 二次元PEM（KV＋定数k）斜面検証プログラム
module simulation_constants_mod
    implicit none
    integer, parameter :: ni_max = 1024
    integer, parameter :: nj_max = 13
    integer, parameter :: nc_max = 20000
    real(8), parameter :: PI_VAL = 3.141592653589793d0
    real(8), parameter :: GRAVITY_ACCEL = 9.80665d0
end module simulation_constants_mod

module simulation_parameters_mod
    use simulation_constants_mod, only: ni_max
    implicit none
    ! 時間・出力
    real(8) :: time_step  ! 時間ステップ
    integer :: max_steps  ! 最大ステップ数
    integer :: output_interval  ! 出力間隔
    ! 斜面
    real(8) :: slope_angle  ! 斜面角度
    real(8) :: friction_coeff  ! 摩擦係数
    ! 接触モデル（定数剛性）
    real(8) :: normal_stiffness_kn  ! 法線剛性
    real(8) :: shear_stiffness_kt  ! 接線剛性
    real(8) :: shear_to_normal_stiffness_ratio  ! 接線剛性と法線剛性の比
    real(8) :: restitution_coeff_normal_en  ! 法線方向の反発係数
    real(8) :: restitution_coeff_tangent_et  ! 接線方向の反発係数
    ! 粒子
    real(8) :: particle_radius  ! 粒子半径
    real(8) :: particle_density  ! 粒子密度
    ! セル法
    real(8) :: container_width  ! コンテナ幅
    logical :: disable_cell_algorithm  ! セル法を無効にするフラグ
    real(8) :: cell_size_override  ! セルサイズのオーバーライド
    ! 初期配置
    real(8) :: initial_x_contact  ! 初期x位置
    save
end module simulation_parameters_mod

module particle_data_mod
    use simulation_constants_mod, only: ni_max, nj_max
    implicit none
    real(8), dimension(ni_max) :: radius, mass, moment_inertia
    real(8), dimension(ni_max) :: x_coord, z_coord, rotation_angle
    real(8), dimension(ni_max) :: x_vel, z_vel, rotation_vel
    real(8), dimension(ni_max) :: x_force_sum, z_force_sum, moment_sum
    real(8), dimension(ni_max, nj_max) :: normal_force_contact, shear_force_contact
    integer, dimension(ni_max, nj_max) :: contact_partner_idx
    real(8), dimension(ni_max, nj_max) :: previous_overlap
    real(8), dimension(ni_max) :: x_disp_incr, z_disp_incr, rot_disp_incr
    save
end module particle_data_mod

module cell_system_mod
    use simulation_constants_mod, only: ni_max, nc_max
    implicit none
    integer :: num_particles
    integer :: cells_x_dir, cells_z_dir
    real(8) :: cell_size
    integer, dimension(nc_max) :: cell_head
    integer, dimension(ni_max) :: particle_cell_next
    integer, dimension(ni_max) :: particle_cell_idx
    integer, dimension(nc_max) :: cell_particle_map
    save
end module cell_system_mod

program pem2d_slope_kv
    use simulation_constants_mod
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    implicit none
    integer :: it_step, i, static_judge_flag
    real(8) :: current_time
    logical :: leapfrog_initialized
    logical :: in_contact
    real(8) :: v_tangent, v_theory, vx_out, vz_out, omega_out
    real(8) :: u_gen_meas, u_gen_theory, vg_theory, ft_theory, ct_value

    call read_input_file
    call fposit_sub
    call inmat_sub
    call init_sub

    current_time = 0.0d0
    leapfrog_initialized = .false.
    in_contact = .false.

    open(unit=10, file='data/pem_slope_trace.csv', status='replace', action='write')
    write(10,'(A)') 'time,x,z,vx,vz,omega,v_tangent,v_theory,error_percent,contact,Fn,Ft,u_generalized,u_generalized_theory,vg_theory,Ft_theory'

    do it_step = 0, max_steps
        current_time = dble(it_step) * time_step
        if (.not. leapfrog_initialized) then ! 初回のみ: 力計算してから v(0) → v(Δt/2) への変換
            call clear_forces
            call wcont_sub(1, in_contact)
            call nposit_leapfrog_sub(0)
            leapfrog_initialized = .true.
        else ! 初回以外: 位置更新してから力計算
            call nposit_leapfrog_sub(1)
            call clear_forces
            call wcont_sub(1, in_contact)
            call nposit_leapfrog_sub(2)
        end if
        ! 理論速度を計算
        call calculate_theoretical_velocity(current_time, v_theory)
        ! 修正速度を計算
        call output_corrected_velocity(1, vx_out, vz_out, omega_out)
        ! 接線速度を計算（補正速度）
        v_tangent = vx_out * cos(slope_angle) + vz_out * sin(slope_angle)

        ! 接線方向の減衰係数 c_t を逐次再計算（壁相手: m_eff = m）
        if (restitution_coeff_tangent_et > 1.0d-6 .and. shear_stiffness_kt > 0.0d0) then
            ct_value = -2.0d0 * log(restitution_coeff_tangent_et) * &
     &                 sqrt(mass(1) * shear_stiffness_kt / (log(restitution_coeff_tangent_et)**2 + PI_VAL**2))
        else
            ct_value = 0.0d0
        end if

        ! 実測の一般化速度 u' = v_t + r*omega（符号は actf_sub の vt_vel に合わせる）
        u_gen_meas = (x_vel(1) * cos(slope_angle) + z_vel(1) * sin(slope_angle)) + radius(1) * rotation_vel(1)

        ! 理論値（一般化変位u, 速度u', 接線力Ft）
        call kv_tangent_theory(current_time, mass(1), shear_stiffness_kt, ct_value, slope_angle, &
     &                        u_gen_theory, vg_theory, ft_theory)

        ! 出力
        call gfout_sub(current_time, v_tangent, v_theory, in_contact, vx_out, vz_out, omega_out, &
     &                u_gen_meas, u_gen_theory, vg_theory, ft_theory)

        if (mod(it_step, 10000) == 0) then
            write(*,'(A,F10.6,A,F12.6)') 'Time=', current_time, ' V_tangent=', v_tangent
        end if
    end do

    close(10)
contains

    !> 入力ファイルを読み込むサブルーチン
    subroutine read_input_file
        implicit none
        character(len=256) :: line, keyword, input_filename
        integer :: ios, unit_num
        real(8) :: value
        logical :: slope_in_deg
        if (command_argument_count() > 0) then
            call get_command_argument(1, input_filename)
            if (input_filename(1:1) /= '/' .and. input_filename(1:5) /= 'input') then
                input_filename = 'input/' // trim(input_filename)
            end if
        else
            input_filename = 'input/slope_input.dat'
        end if
        unit_num = 20
        open(unit=unit_num, file=input_filename, status='old', action='read', iostat=ios)
        if (ios /= 0) stop '入力ファイルを開けません'
        ! デフォルト
        time_step = 1.0d-5
        max_steps = 200000
        output_interval = 200
        slope_angle = 30.0d0 * PI_VAL / 180.0d0
        friction_coeff = 0.2d0
        normal_stiffness_kn = 1.0d5
        shear_stiffness_kt = 5.0d4
        shear_to_normal_stiffness_ratio = 0.5d0
        restitution_coeff_normal_en = 0.8d0
        restitution_coeff_tangent_et = 0.8d0
        particle_radius = 1.0d-2
        particle_density = 2500.0d0
        initial_x_contact = 2.0d-2
        disable_cell_algorithm = .false.
        cell_size_override = 0.0d0
        container_width = 0.2d0
        do
            read(unit_num,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            read(line, *, iostat=ios) keyword, value
            if (ios /= 0) cycle
            select case (trim(keyword))
            case ('TIME_STEP'); time_step = value
            case ('MAX_STEPS'); max_steps = int(value)
            case ('OUTPUT_INTERVAL'); output_interval = int(value)
            case ('SLOPE_ANGLE'); slope_angle = value * PI_VAL / 180.0d0
            case ('FRICTION_COEFF'); friction_coeff = value
            case ('KN'); normal_stiffness_kn = value
            case ('KT'); shear_stiffness_kt = value
            case ('SHEAR_TO_NORMAL_RATIO'); shear_to_normal_stiffness_ratio = value
            case ('E_NORMAL'); restitution_coeff_normal_en = value 
            case ('E_TANGENT'); restitution_coeff_tangent_et = value
            case ('PARTICLE_RADIUS'); particle_radius = value
            case ('PARTICLE_DENSITY'); particle_density = value
            case ('INITIAL_X_CONTACT'); initial_x_contact = value
            case ('DISABLE_CELL_ALGORITHM'); disable_cell_algorithm = (int(value) == 1)
            case ('CELL_SIZE_OVERRIDE'); cell_size_override = value
            case default
            end select
        end do
        close(unit_num)
        if (shear_stiffness_kt <= 0.0d0) shear_stiffness_kt = max(1.0d0, normal_stiffness_kn * shear_to_normal_stiffness_ratio)
    end subroutine read_input_file

    !> 粒子を配置するサブルーチン
    subroutine fposit_sub
        use simulation_constants_mod, only: PI_VAL
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod
        implicit none
        real(8) :: n_x, n_z, d
        num_particles = 1
        radius(1) = particle_radius
        x_coord(1) = initial_x_contact - radius(1) * sin(slope_angle)
        z_coord(1) = tan(slope_angle) * initial_x_contact + radius(1) * cos(slope_angle)
        rotation_angle(1) = 0.0d0
        x_vel(1) = 0.0d0
        z_vel(1) = 0.0d0
        rotation_vel(1) = 0.0d0
        if (cell_size_override > 0.0d0) then
            cell_size = cell_size_override
        else
            cell_size = radius(1) * 2.0d0
        end if
        container_width = max(0.2d0, initial_x_contact * 2.0d0)
        cells_x_dir = idint(container_width / cell_size) + 1
        cells_z_dir = idint((z_coord(1) * 2.0d0) / cell_size) + 10
        if (cells_x_dir * cells_z_dir > nc_max) stop 'セル数が多すぎます'
    end subroutine fposit_sub

    !> 質量と慣性モーメントを計算するサブルーチン
    subroutine inmat_sub
        use simulation_constants_mod, only: PI_VAL
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer :: i
        do i = 1, num_particles
            mass(i) = (4.0d0 / 3.0d0) * PI_VAL * radius(i)**3 * particle_density
            moment_inertia(i) = (2.0d0 / 5.0d0) * mass(i) * radius(i)**2
        end do
    end subroutine inmat_sub

    !> 初期化するサブルーチン
    subroutine init_sub
        use simulation_constants_mod, only: nj_max
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer :: i, j
        do i = 1, num_particles
            do j = 1, nj_max
                normal_force_contact(i, j) = 0.0d0
                shear_force_contact(i, j) = 0.0d0
                contact_partner_idx(i, j) = 0
                previous_overlap(i, j) = -1.0d0
            end do
            x_disp_incr(i) = 0.0d0
            z_disp_incr(i) = 0.0d0
            rot_disp_incr(i) = 0.0d0
        end do
    end subroutine init_sub

    !> 位置の更新を行うサブルーチン
    subroutine nposit_leapfrog_sub(phase)
        use simulation_parameters_mod, only: time_step
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: phase
        integer :: i
        real(8) :: dt
        dt = time_step
        if (phase == 0) then
            do i = 1, num_particles
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * (dt * 0.5d0)
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - GRAVITY_ACCEL) * (dt * 0.5d0)
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * (dt * 0.5d0)
            end do
        else if (phase == 1) then
            do i = 1, num_particles
                x_disp_incr(i) = x_vel(i) * dt
                z_disp_incr(i) = z_vel(i) * dt
                rot_disp_incr(i) = rotation_vel(i) * dt
                x_coord(i) = x_coord(i) + x_disp_incr(i)
                z_coord(i) = z_coord(i) + z_disp_incr(i)
                rotation_angle(i) = rotation_angle(i) + rot_disp_incr(i)
            end do
        else if (phase == 2) then
            do i = 1, num_particles
                x_vel(i) = x_vel(i) + (x_force_sum(i)/mass(i)) * dt
                z_vel(i) = z_vel(i) + (z_force_sum(i)/mass(i) - GRAVITY_ACCEL) * dt
                rotation_vel(i) = rotation_vel(i) + (moment_sum(i)/moment_inertia(i)) * dt
            end do
        end if
    end subroutine nposit_leapfrog_sub

    !> 接触力をクリアするサブルーチン
    subroutine clear_forces
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer :: i
        do i = 1, num_particles
            x_force_sum(i) = 0.0d0
            z_force_sum(i) = 0.0d0
            moment_sum(i) = 0.0d0
        end do
    end subroutine clear_forces

    !> 壁接触を判定するサブルーチン
    subroutine wcont_sub(particle_idx, contact_flag)
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: particle_idx
        logical, intent(out) :: contact_flag
        real(8) :: xi, zi, ri
        real(8) :: n_x, n_z, d, overlap
        real(8) :: angle_sin, angle_cos
        integer :: wall_slot, wall_id
        logical :: use_rigid_model
        ! 剛体分岐用の補助変数（ここで宣言）
        real(8) :: tol_on, tol_off, corr, vn
        logical :: was_in_contact, contact_now
        xi = x_coord(particle_idx)
        ! write(*,*) 'xi = ', xi
        zi = z_coord(particle_idx)
        ! write(*,*) 'zi = ', zi
        ri = radius(particle_idx)
        n_x = -sin(slope_angle) ! 壁の単位法線ベクトル
        n_z =  cos(slope_angle) ! 壁の単位法線ベクトル
        d = n_x * xi + n_z * zi ! 粒子と壁の距離
        ! write(*,*) 'd = ', d
        angle_sin = sin(slope_angle)
        angle_cos = cos(slope_angle)
        overlap = ri - abs(d) ! 粒子と壁の重なり
        ! write(*,*) 'overlap = ', overlap
        wall_slot = 10 ! 壁のスロット
        wall_id = num_particles + 4 ! 壁のID
        use_rigid_model = .false. ! (restitution_coeff_normal_en >= 0.999d0)

        if (use_rigid_model) then
            ! 剛体分岐: 接触ヒステリシスと運動学的拘束（法線射影＋位置投影）
            tol_on = 1.0d-9
            tol_off = 2.0d-9
            was_in_contact = (previous_overlap(particle_idx, wall_slot) >= 0.0d0)
            if (was_in_contact) then
                contact_now = (overlap > -tol_off)
            else
                contact_now = (overlap > -tol_on)
            end if

            if (contact_now) then
                ! 接触している: 先に運動学的拘束を適用
                corr = ri - d
                if (abs(corr) > 0.0d0) then
                    x_coord(particle_idx) = x_coord(particle_idx) + corr * n_x
                    z_coord(particle_idx) = z_coord(particle_idx) + corr * n_z
                end if
                vn = x_vel(particle_idx) * n_x + z_vel(particle_idx) * n_z
                if (abs(vn) > 0.0d0) then
                    x_vel(particle_idx) = x_vel(particle_idx) - vn * n_x
                    z_vel(particle_idx) = z_vel(particle_idx) - vn * n_z
                end if

                contact_partner_idx(particle_idx, wall_slot) = wall_id
                call actf_sub_rigid(particle_idx, wall_id, wall_slot, angle_sin, angle_cos)
                contact_flag = .true.
            else
                ! 接触していない
                normal_force_contact(particle_idx, wall_slot) = 0.0d0
                shear_force_contact(particle_idx, wall_slot) = 0.0d0
                contact_partner_idx(particle_idx, wall_slot) = 0
                previous_overlap(particle_idx, wall_slot) = -1.0d0
                contact_flag = .false.
            end if
        else
            if (overlap > 0.0d0) then
                ! 接触している（弾性分岐）
                contact_partner_idx(particle_idx, wall_slot) = wall_id
                call actf_sub(particle_idx, wall_id, wall_slot, angle_sin, angle_cos, overlap)
                contact_flag = .true.
            else
                ! 接触していない
                normal_force_contact(particle_idx, wall_slot) = 0.0d0
                shear_force_contact(particle_idx, wall_slot) = 0.0d0
                contact_partner_idx(particle_idx, wall_slot) = 0
                previous_overlap(particle_idx, wall_slot) = -1.0d0
                contact_flag = .false.
            end if
        end if
    end subroutine wcont_sub

    !> 接触力を計算するサブルーチン
    subroutine actf_sub(p_i, p_j, slot, angle_sin, angle_cos, initial_overlap)
        use simulation_constants_mod, only: ni_max, PI_VAL
        use simulation_parameters_mod
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none
        integer, intent(in) :: p_i, p_j, slot
        real(8), intent(in) :: angle_sin, angle_cos, initial_overlap
        real(8) :: ri, rj, m_eff
        real(8) :: kn, kt
        real(8) :: cn, ct
        real(8) :: rel_n, rel_t
        real(8) :: dn, dtan
        real(8) :: fn_total, ft_total
        real(8) :: xn_i, zn_i, xn_j, zn_j
        real(8) :: rel_vx, rel_vz, vn_vel, vt_vel
        real(8) :: nx, nz, tx, tz
        ri = radius(p_i)
        if (p_j <= num_particles) then ! 他の粒子であれば
            rj = radius(p_j)
        else ! 壁であれば
            rj = 0.0d0 
        end if
        ! 質量を計算
        if (p_j <= num_particles) then
            m_eff = mass(p_i) * mass(p_j) / max(1.0d-30, mass(p_i) + mass(p_j))
        else
            m_eff = mass(p_i)
        end if
        ! 剛性を計算
        kn = normal_stiffness_kn
        if (shear_stiffness_kt > 0.0d0) then
            kt = shear_stiffness_kt
            ! write(*,*) 'kt = shear_stiffness_kt'
        else
            kt = kn * shear_to_normal_stiffness_ratio
            ! write(*,*) 'kt = kn * shear_to_normal_stiffness_ratio'
        end if
        ! write(*,*) 'kn = ', kn, ' kt = ', kt
        ! 法線方向の減衰係数を計算
        if (restitution_coeff_normal_en > 1.0d-6 .and. kn > 0.0d0) then
            cn = -2.0d0 * log(restitution_coeff_normal_en) * sqrt(m_eff * kn / (log(restitution_coeff_normal_en)**2 + PI_VAL**2))
        else
            cn = 0.0d0
        end if
        ! 接線方向の減衰係数を計算
        if (restitution_coeff_tangent_et > 1.0d-6 .and. kt > 0.0d0) then ! 反発係数が正で剛性が正の場合
            ct = -2.0d0 * log(restitution_coeff_tangent_et) * sqrt(m_eff * kt / (log(restitution_coeff_tangent_et)**2 + PI_VAL**2))
        else
            ct = 0.0d0
        end if
        ! write(*,*) 'restitution_coeff_normal_en = ', restitution_coeff_normal_en, ' restitution_coeff_tangent_et = ', restitution_coeff_tangent_et
        ! 法線・接線ベクトル
        nx = -angle_sin
        nz =  angle_cos
        tx =  angle_cos
        tz =  angle_sin

        ! 相対位置（法線方向は圧縮を正とする）
        rel_n = -(x_disp_incr(p_i) * nx + z_disp_incr(p_i) * nz)
        rel_t =  (x_disp_incr(p_i) * tx + z_disp_incr(p_i) * tz) + (ri * rot_disp_incr(p_i))
        if (p_j <= num_particles) then
            rel_t = rel_t + rj * rot_disp_incr(p_j)
        end if
        ! write(*,*) 'x_disp_incr(p_i) = ', x_disp_incr(p_i), ' z_disp_incr(p_i) = ', z_disp_incr(p_i)
        ! write(*,*) 'rel_n = ', rel_n, ' rel_t = ', rel_t
        ! 弾性力を計算
        if (abs(normal_force_contact(p_i, slot)) < 1.0d-12) then ! 法線力が小さい場合は初期接触力を設定
            normal_force_contact(p_i, slot) = kn * initial_overlap
            shear_force_contact(p_i, slot) = 0.0d0
        else ! 法線力が大きい場合は接触力を更新
            normal_force_contact(p_i, slot) = normal_force_contact(p_i, slot) + kn * rel_n
            shear_force_contact(p_i, slot) = shear_force_contact(p_i, slot) + kt * rel_t
        end if
        ! 相対速度を計算
        if (p_j <= num_particles) then
            rel_vx = x_vel(p_i) - x_vel(p_j)
            rel_vz = z_vel(p_i) - z_vel(p_j)
        else
            rel_vx = x_vel(p_i)
            rel_vz = z_vel(p_i)
        end if
        ! 法線方向の相対速度
        vn_vel = -(rel_vx * nx + rel_vz * nz)
        ! 接線方向の相対速度
        if (p_j <= num_particles) then ! 他の粒子であれば
            vt_vel = rel_vx * tx + rel_vz * tz + (ri * rotation_vel(p_i) + rj * rotation_vel(p_j))
        else ! 壁であれば
            vt_vel = rel_vx * tx + rel_vz * tz + (ri * rotation_vel(p_i))
        end if
        ! 減衰力を計算
        dn = cn * vn_vel
        dtan = ct * vt_vel
        ! 合計力を計算
        fn_total = normal_force_contact(p_i, slot) + dn
        ft_total = shear_force_contact(p_i, slot) + dtan

        ! 法線力が負の場合は接触を解除
        if (fn_total <= 0.0d0) then
            fn_total = 0.0d0
            ft_total = 0.0d0
            contact_partner_idx(p_i, slot) = 0
            normal_force_contact(p_i, slot) = 0.0d0
            shear_force_contact(p_i, slot) = 0.0d0
            previous_overlap(p_i, slot) = -1.0d0
            return
        end if

        ! せん断力が大きすぎる場合は摩擦力による制御を行う
        if (abs(ft_total) > friction_coeff * fn_total) then
            ft_total = friction_coeff * fn_total * sign(1.0d0, ft_total)
            ! write(*,*) 'ft_total = friction_coeff * fn_total * sign(1.0d0, ft_total)'
        end if
        ! 力の合計を計算
        x_force_sum(p_i) = x_force_sum(p_i) + fn_total * nx - ft_total * tx
        z_force_sum(p_i) = z_force_sum(p_i) + fn_total * nz - ft_total * tz
        moment_sum(p_i) = moment_sum(p_i) - ri * ft_total
        previous_overlap(p_i, slot) = initial_overlap
    end subroutine actf_sub

    !> 剛体接触力（法線・接線方向）を計算するサブルーチン
    subroutine actf_sub_rigid(p_i, p_j, slot, angle_sin, angle_cos)
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use simulation_parameters_mod, only: slope_angle, friction_coeff
        use particle_data_mod
        implicit none
        integer, intent(in) :: p_i, p_j, slot
        real(8), intent(in) :: angle_sin, angle_cos
        real(8) :: ri
        real(8) :: fn_total, ft_total
        real(8) :: tangent_x, tangent_z
        real(8) :: v_tangent
        real(8) :: fn_limit
        real(8) :: desired_static
        real(8) :: slip_sign
        real(8) :: force_tangent_gravity
        real(8) :: v_eps
        real(8) :: mu_c

        ri = radius(p_i)

        ! 法線・接線方向単位ベクトル
        tangent_x =  angle_cos
        tangent_z =  angle_sin

        ! 法線反力（斜面上の剛体仮定）
        fn_total = mass(p_i) * GRAVITY_ACCEL * cos(slope_angle)
        fn_total = max(fn_total, 0.0d0)

        ! 重力の接線成分（正味方向の符号決定に使用）
        force_tangent_gravity = mass(p_i) * GRAVITY_ACCEL * sin(slope_angle)

        fn_limit = friction_coeff * fn_total

        ! 接線方向速度（転がりを含む）
        v_tangent = x_vel(p_i) * tangent_x + z_vel(p_i) * tangent_z + ri * rotation_vel(p_i)

        ! 臨界摩擦係数（剛体球, I = 2/5 m r^2）
        mu_c = (2.0d0 / 7.0d0) * tan(slope_angle)

        if (friction_coeff >= mu_c - 1.0d-12) then
            ! 転がり（非すべり）: Ft = -(2/7) m g sinθ（符号は斜面下向きに抗する）
            ft_total = - (2.0d0 / 7.0d0) * mass(p_i) * GRAVITY_ACCEL * sin(slope_angle)
            ! 念のためクーロン限界は超えない（μ >= μ_c なら満たす）
            if (abs(ft_total) > fn_limit) then
                ft_total = fn_limit * sign(1.0d0, ft_total)
            end if
        else
            ! すべり: クーロン限界
            v_eps = 1.0d-9
            if (abs(v_tangent) > v_eps) then
                slip_sign = sign(1.0d0, v_tangent)
            else
                slip_sign = -sign(1.0d0, force_tangent_gravity)
            end if
            ft_total = fn_limit * slip_sign
        end if

        normal_force_contact(p_i, slot) = fn_total
        shear_force_contact(p_i, slot) = ft_total
        previous_overlap(p_i, slot) = 0.0d0

        x_force_sum(p_i) = x_force_sum(p_i) + fn_total * (-angle_sin) - ft_total * tangent_x
        z_force_sum(p_i) = z_force_sum(p_i) + fn_total * angle_cos - ft_total * tangent_z
        moment_sum(p_i) = moment_sum(p_i) - ri * ft_total
    end subroutine actf_sub_rigid

    !> 理論速度を計算するサブルーチン
    subroutine calculate_theoretical_velocity(t, v_theory)
        use simulation_constants_mod, only: GRAVITY_ACCEL, PI_VAL
        use simulation_parameters_mod, only: slope_angle, friction_coeff
        implicit none
        real(8), intent(in) :: t
        real(8), intent(out) :: v_theory
        real(8) :: mu_c, accel
        ! 臨界摩擦係数を計算
        mu_c = (2.0d0 / 7.0d0) * tan(slope_angle)
        ! write(*,*) 'mu_c = ', mu_c
        ! 摩擦係数が臨界摩擦係数より小さい場合は加速度を計算
        if (friction_coeff < mu_c - 1.0d-6) then ! 摩擦係数が臨界摩擦係数より小さい場合は滑りの理論加速度を計算
            ! 滑り条件: a = g(sin(θ) - μ*cos(θ))
            accel = - GRAVITY_ACCEL * (sin(slope_angle) - friction_coeff * cos(slope_angle))
            write(*,*) 'accel = sliding'
            ! write(*,*) 'friction_coeff = ', friction_coeff, ' mu_c = ', mu_c
            ! write(*,*) 'accel = ', accel
        else ! 摩擦係数が臨界摩擦係数より大きい場合は転がりの理論加速度を計算
            ! 転がり条件: a = g*sin(θ) / (1 + I/(m*r²))
            ! I/(m*r²) = 2/5 なので、a = (5/7) * g * sin(θ)
            accel = - (5.0d0 / 7.0d0) * GRAVITY_ACCEL * sin(slope_angle)
            write(*,*) 'accel = rolling'
            ! write(*,*) 'accel = ', accel
        end if
        v_theory = accel * t
    end subroutine calculate_theoretical_velocity

    !> 接線KV系の理論解（一般化変位u, 速度u' と接線力Ft）を返す
    subroutine kv_tangent_theory(t, mass_val, kt, ct, theta, u_th, up_th, ft_th)
        use simulation_constants_mod, only: GRAVITY_ACCEL
        implicit none
        real(8), intent(in) :: t, mass_val, kt, ct, theta
        real(8), intent(out) :: u_th, up_th, ft_th
        real(8) :: m_eq, omega_n, zeta, omega_d
        real(8) :: u_ss
        real(8) :: expfac, s, c, alpha, zeta_sq, root_term
        real(8) :: r1, r2, A, B

        m_eq = (2.0d0 / 7.0d0) * mass_val
        if (kt <= 0.0d0 .or. m_eq <= 0.0d0) then
            u_th = 0.0d0
            up_th = 0.0d0
            ft_th = 0.0d0
            return
        end if

        omega_n = sqrt(kt / m_eq)
        if (ct > 0.0d0) then
            zeta = ct / (2.0d0 * sqrt(kt * m_eq))
        else
            zeta = 0.0d0
        end if
        zeta_sq = zeta * zeta

        u_ss = (GRAVITY_ACCEL * sin(theta)) * m_eq / kt

        if (zeta < 1.0d0 - 1.0d-12) then
            omega_d = omega_n * sqrt(max(0.0d0, 1.0d0 - zeta_sq))
            alpha = zeta / max(1.0d-30, sqrt(max(0.0d0, 1.0d0 - zeta_sq)))
            expfac = exp(-zeta * omega_n * t)
            s = sin(omega_d * t)
            c = cos(omega_d * t)
            u_th = u_ss - u_ss * expfac * (c + alpha * s)
            up_th = -u_ss * expfac * ( -omega_d * s + alpha * omega_d * c - zeta * omega_n * (c + alpha * s) )
        else if (zeta > 1.0d0 + 1.0d-12) then
            root_term = sqrt(zeta_sq - 1.0d0)
            r1 = -omega_n * (zeta - root_term)
            r2 = -omega_n * (zeta + root_term)
            A = u_ss * r2 / (r2 - r1)
            B = -u_ss * r1 / (r2 - r1)
            u_th = u_ss + A * exp(r1 * t) + B * exp(r2 * t)
            up_th = A * r1 * exp(r1 * t) + B * r2 * exp(r2 * t)
        else
            expfac = exp(-omega_n * t)
            u_th = u_ss - u_ss * (1.0d0 + omega_n * t) * expfac
            up_th = u_ss * (omega_n * omega_n * t) * expfac
        end if

        ft_th = kt * u_th + ct * up_th
    end subroutine kv_tangent_theory

    !> 補正速度を出力するサブルーチン
    subroutine output_corrected_velocity(p_i, vx_out, vz_out, omega_out)
        use simulation_parameters_mod, only: time_step
        use simulation_constants_mod, only: GRAVITY_ACCEL
        use particle_data_mod
        implicit none
        integer, intent(in) :: p_i
        real(8), intent(out) :: vx_out, vz_out, omega_out
        vx_out = x_vel(p_i) - 0.5d0 * time_step * (x_force_sum(p_i) / mass(p_i))
        vz_out = z_vel(p_i) - 0.5d0 * time_step * (z_force_sum(p_i) / mass(p_i) - GRAVITY_ACCEL)
        omega_out = -(rotation_vel(p_i) - 0.5d0 * time_step * (moment_sum(p_i) / moment_inertia(p_i)))
    end subroutine output_corrected_velocity

    !> データを出力するサブルーチン
    subroutine gfout_sub(t, v_tangent, v_theory, contact, vx_out, vz_out, omega_out, &
         u_gen_meas, u_gen_theory, vg_theory, ft_theory)
        use simulation_parameters_mod, only: slope_angle
        use particle_data_mod
        implicit none
        real(8), intent(in) :: t, v_tangent, v_theory, vx_out, vz_out, omega_out
        real(8), intent(in) :: u_gen_meas, u_gen_theory, vg_theory, ft_theory
        logical, intent(in) :: contact
        real(8) :: error_pct
        integer :: c_int
        real(8) :: nx, nz, fn, ft
        if (abs(v_theory) > 1.0d-12) then
            error_pct =  abs(v_tangent - v_theory) / v_theory * 100.0d0
        else
            error_pct = 0.0d0
        end if
        c_int = merge(1, 0, contact)
        nx = -sin(slope_angle)
        nz =  cos(slope_angle)
        fn =  (nx * x_force_sum(1) + nz * z_force_sum(1))
        ft =  (nz * x_force_sum(1) - nx * z_force_sum(1))
        write(10,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,F10.4,A,I1,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
            t, ',', x_coord(1), ',', z_coord(1), ',', vx_out, ',', vz_out, ',', omega_out, ',', &
            v_tangent, ',', v_theory, ',', error_pct, ',', c_int, ',', fn, ',', ft, ',', &
            u_gen_meas, ',', u_gen_theory, ',', vg_theory, ',', ft_theory
    end subroutine gfout_sub

end program pem2d_slope_kv
