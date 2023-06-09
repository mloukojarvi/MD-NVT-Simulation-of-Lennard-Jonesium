program main
  ! The main program of the simulation
  use const
  use io
  use simulation
  implicit none

  ! target temperature, density, simulation box edge length
  real(rk) :: temp, dens, L
  ! time, instantaneous temperature, kinetic energy, potential energy, total energy, mean squared displacement
  real(rk) :: t, inst_temp, KE, PE, TE, msd
  ! positions, velocities, forces, intial positions
  real(rk), allocatable :: pos(:,:), vel(:,:), force(:,:), init_pos(:,:)
  ! image coordinates, radial distribution function, speed distribution, potential energies of each particle
  real(rk), allocatable :: img(:,:), rdf(:), spd(:), pepp(:)
  ! list of the filenames
  character(len=29) :: files(5)
  ! current simulation step
  integer(ik) :: step
  
  ! Initialize the simulation and the files
  call read_params(N, temp, dens, tot_steps, dt, gamma_)
  call set_params(N, tot_steps, dt, gamma_)
  allocate(pos(N,3), vel(N,3), force(N,3), init_pos(N,3), img(N,3), rdf(rdf_bins), spd(spd_bins), pepp(N))
  call initialize(temp, dens, L, inst_temp, PE, KE, msd, pos, vel, force, init_pos, img, rdf, spd, pepp)
  call init_files(dens, temp, files)
  
  
  ! Run the simulation
  do step=1, tot_steps
    t = (step-1)*dt
    call info(step, t, temp, dens, PE, KE, inst_temp)
    call write_data(files, step, t, temp, PE, KE, inst_temp, pos, spd, rdf, msd, pepp)
    call simulate(L, temp, pos, vel, force, PE, KE, inst_temp, msd, init_pos, img, rdf, spd, pepp)
  end do
  
end program main
