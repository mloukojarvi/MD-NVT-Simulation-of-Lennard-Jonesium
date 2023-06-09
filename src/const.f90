module const
  ! Constants and constant parameters used in the simulation
  implicit none
  
  integer, parameter :: ik=8, rk=8 ! Sets the integer and real kinds
  real(rk), parameter :: PI = 3.1415926535
  
  ! number of bins for the rdf and speed histograms
  integer(ik) :: rdf_bins = 90, spd_bins = 60 
  ! maximum distance of the interaction, minimum distance for 
  ! the rdf histogram, maximum distance for the speed histogram
  real(rk) :: rc = 3.0, r_min = 0.8, sc = 6.0
  ! number of particles, total number of simulation steps
  integer(ik) :: N, tot_steps
  ! length of time step, thermostat coefficient
  real(rk) :: dt, gamma_
  
contains
  
  subroutine set_params(set_N, set_tot_steps, set_dt, set_gamma_)
  	! Sets the given constant input parameters
  	integer(ik) :: set_N, set_tot_steps
  	real(rk) :: set_dt, set_gamma_
  	
  	N = set_N
  	tot_steps = set_tot_steps
  	dt = set_dt
  	gamma_ = set_gamma_
  end subroutine set_params

end module const
