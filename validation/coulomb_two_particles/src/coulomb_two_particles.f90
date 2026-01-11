module simulation_parameters_mod
    implicit none
    real(8) :: coulomb_constant = 8.9875517923d9
    logical :: enable_coulomb_force = .true.
end module simulation_parameters_mod

module particle_data_mod
    implicit none
    integer, parameter :: ni_max = 2
    real(8), dimension(ni_max) :: x_coord
    real(8), dimension(ni_max) :: z_coord
    real(8), dimension(ni_max) :: charge
    real(8), dimension(ni_max) :: x_force_sum
    real(8), dimension(ni_max) :: z_force_sum
end module particle_data_mod

module cell_system_mod
    implicit none
    integer :: num_particles = 0
end module cell_system_mod

program coulomb_two_particles
    use simulation_parameters_mod
    use particle_data_mod
    use cell_system_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: step, n_steps, output_stride, unit_out
    real(dp) :: dt, time_val
    real(dp), dimension(ni_max) :: mass, vx_half, vz_half
    real(dp) :: separation0, charge_magnitude
    character(len=*), parameter :: output_file = 'data/coulomb_two_particles.csv'

    num_particles = 2
    dt = 1.0d-6
    n_steps = 200000
    output_stride = 100
    mass = 1.0d-6
    charge_magnitude = 5.0d-8
    charge = charge_magnitude
    separation0 = 2.0d-2
    x_coord(1) = -0.5d0 * separation0
    x_coord(2) =  0.5d0 * separation0
    z_coord = 0.0d0
    vx_half = 0.0d0
    vz_half = 0.0d0

    call zero_forces()
    call coulomb_force_sub()
    call initialize_leapfrog_velocities(dt, mass, vx_half, vz_half)

    open(newunit=unit_out, file=output_file, status='replace', action='write')
    write(unit_out,'(A)') '# mass=' // trim(format_real(mass(1))) // &
                          ',charge=' // trim(format_real(charge_magnitude)) // &
                          ',coulomb_constant=' // trim(format_real(coulomb_constant)) // &
                          ',separation0=' // trim(format_real(separation0)) // &
                          ',time_step=' // trim(format_real(dt)) // &
                          ',steps=' // trim(format_real(real(n_steps, dp)))
    write(unit_out,'(A)') 'time,x1,x2,vx1,vx2,separation'
    time_val = 0.0d0
    call write_state(unit_out, time_val, vx_half, mass, dt)

    do step = 1, n_steps
        call leapfrog_step(dt, mass, vx_half, vz_half)
        time_val = time_val + dt
        if (mod(step, output_stride) == 0) then
            call write_state(unit_out, time_val, vx_half, mass, dt)
        end if
    end do

    close(unit_out)

contains

    subroutine initialize_leapfrog_velocities(dt, mass, vx_half, vz_half)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: vx_half(:)
        real(dp), intent(inout) :: vz_half(:)
        integer :: i

        do i = 1, num_particles
            vx_half(i) = vx_half(i) + 0.5d0 * dt * x_force_sum(i) / mass(i)
            vz_half(i) = vz_half(i) + 0.5d0 * dt * z_force_sum(i) / mass(i)
        end do
    end subroutine initialize_leapfrog_velocities

    subroutine leapfrog_step(dt, mass, vx_half, vz_half)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: mass(:)
        real(dp), intent(inout) :: vx_half(:)
        real(dp), intent(inout) :: vz_half(:)
        integer :: i

        do i = 1, num_particles
            x_coord(i) = x_coord(i) + dt * vx_half(i)
            z_coord(i) = z_coord(i) + dt * vz_half(i)
        end do

        call zero_forces()
        call coulomb_force_sub()

        do i = 1, num_particles
            vx_half(i) = vx_half(i) + dt * x_force_sum(i) / mass(i)
            vz_half(i) = vz_half(i) + dt * z_force_sum(i) / mass(i)
        end do
    end subroutine leapfrog_step

    subroutine zero_forces()
        x_force_sum = 0.0d0
        z_force_sum = 0.0d0
    end subroutine zero_forces

    subroutine write_state(unit, time_val, vx_half, mass, dt)
        integer, intent(in) :: unit
        real(dp), intent(in) :: time_val
        real(dp), intent(in) :: vx_half(:)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: dt
        real(dp) :: separation
        real(dp) :: vx_out(ni_max)
        integer :: i

        do i = 1, num_particles
            vx_out(i) = vx_half(i) - 0.5d0 * dt * x_force_sum(i) / mass(i)
        end do

        separation = x_coord(2) - x_coord(1)
        write(unit,'(ES23.15,",",ES23.15,",",ES23.15,",",ES23.15,",",ES23.15,",",ES23.15)') &
             time_val, x_coord(1), x_coord(2), vx_out(1), vx_out(2), separation
    end subroutine write_state

    function format_real(value) result(str)
        real(dp), intent(in) :: value
        character(len=32) :: str
        write(str,'(ES16.8)') value
        str = adjustl(str)
    end function format_real

    subroutine coulomb_force_sub
        use simulation_parameters_mod, only: coulomb_constant, enable_coulomb_force
        use particle_data_mod
        use cell_system_mod, only: num_particles
        implicit none

        integer :: i, j
        real(8) :: dx, dz, dist, dist_sq, dist_cubed
        real(8) :: force_magnitude, fx, fz
        real(8) :: qi, qj

        if (.not. enable_coulomb_force) return

        !$omp parallel do schedule(dynamic) private(i, j, qi, qj, dx, dz, dist_sq, dist, dist_cubed, force_magnitude, fx, fz)
        do i = 1, num_particles - 1
            qi = charge(i)
            if (abs(qi) < 1.0d-20) cycle

            do j = i + 1, num_particles
                qj = charge(j)
                if (abs(qj) < 1.0d-20) cycle

                dx = x_coord(j) - x_coord(i)
                dz = z_coord(j) - z_coord(i)
                dist_sq = dx*dx + dz*dz

                if (dist_sq < 1.0d-20) cycle

                dist = sqrt(dist_sq)
                dist_cubed = dist * dist_sq

                force_magnitude = coulomb_constant * qi * qj / dist_cubed

                fx = force_magnitude * dx
                fz = force_magnitude * dz

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
    end subroutine coulomb_force_sub

end program coulomb_two_particles

