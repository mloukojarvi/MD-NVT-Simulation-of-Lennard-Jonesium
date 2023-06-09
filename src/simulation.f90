module simulation
  ! Methods for running the simulation
  use const
  implicit none
  
contains

  function box_muller(mu, sigma)
    ! Returns a triplet of random numbers taken from a Gaussian distribution
    ! where the mean and standard deviation are given as arguments
    real(rk), intent(in) :: sigma, mu
    real(rk) :: rn(4), box_muller(3)
    
    call random_number(rn)
    box_muller(1) = mu + sigma*sqrt(-2*log(rn(1)))*cos(2*PI*rn(2))
    box_muller(2) = mu + sigma*sqrt(-2*log(rn(1)))*sin(2*PI*rn(2))
    box_muller(3) = mu + sigma*sqrt(-2*log(rn(3)))*sin(2*PI*rn(4))
  end function box_muller


  subroutine initialize(temp, dens, L, inst_temp, PE, KE, msd, pos, vel, force, init_pos, img, rdf, spd, pepp)
    ! Initializes the simulation by setting the particles into a FCC-lattice  and gives the particle velocities
    ! according to the Maxwell-Boltzmann distribution where the temperature is given as an argument
    real(rk), intent(in) :: temp, dens
    real(rk), intent(inout) :: L, inst_temp, PE, KE, msd
    real(rk), intent(inout) :: pos(N,3), vel(N,3), force(N, 3), init_pos(N,3), img(N,3), rdf(rdf_bins), spd(spd_bins), pepp(N)
    real(rk) :: d, start, pos1(3), pos2(3), avg_velocity(3), r, dr
    integer(ik) :: M, i, j, k, idx
    
    L = (N/dens)**(1.0/3.0)
    M = (N*2)**(1.0/3.0) ! Number of particles in a row
    d = 2*L/M ! Distance between two adjacent particles
    idx = 1
    do i=1, M/2
    	do j=1, M
    	  do k=1, M
    	  	start = 0
    	    if ((mod(k,2) == 0 .and. mod(j,2) == 0) .or. (mod(k,2) .ne. 0 .and.  mod(j,2) .ne. 0)) then
    	      start = d/2.0
    	    end if
    	    pos(idx,:) = (L/2.0 - d/4.0)*[1.0, 1.0, 1.0]
    	    pos(idx,:) = pos(idx,:) - [start + (i-1)*d, (j-1)*d/2.0, (k-1)*d/2.0]
    	    idx = idx + 1 ! Unique index for each particle
    	  end do
    	end do
    end do
    
    img = 0.0
    init_pos = pos
    
    call init_step(L, temp, inst_temp, vel, force, msd, PE, KE, rdf, spd, pepp)
    
    avg_velocity = 0.0
    do i=1, N
      vel(i,:) = box_muller(0.0_rk, sqrt(temp))
      avg_velocity = avg_velocity + 1.0_rk/N*vel(i,:)
    end do
    ! Assert zero initial momentum
    do i=1, N
      vel(i,:) = vel(i,:) - avg_velocity
    end do
    call update(L, pos, vel, force, PE, KE, inst_temp, msd, init_pos, img, rdf, spd, pepp)
  end subroutine initialize
  
  
  function dist(pos1, pos2)
    ! Returns the euclidean distance between the positions given as arguments
    real(rk), intent(in) :: pos1(3), pos2(3)
    real(rk) :: dist
    
    dist = sqrt((sum((pos2-pos1)**2)))
  end function dist
  
  
  function eff_pos(pos1, pos2, L)
    ! Returns the effective position of particle 2 with respect to particle 1 
    ! taking into account the periodic boundary conditions
    real(rk), intent(in) :: pos1(3), pos2(3), L
    real(rk) :: rel_pos(3), eff_pos(3)
    
    rel_pos = pos2-pos1
    eff_pos = pos2-L*nint(rel_pos/L)
  end function eff_pos
  
  
  subroutine init_step(L, temp, inst_temp, vel, force, msd, PE, KE, rdf, spd, pepp)
    ! Initializes the variables before each simulation step
    real(rk), intent(in) :: L, temp, inst_temp, vel(N,3)
    real(rk), intent(inout) :: PE, KE, force(N,3), msd, rdf(rdf_bins), spd(spd_bins), pepp(N)
    real(rk) :: sigma
    integer(ik) :: i
    
    PE = 0.0
    KE = 0.0
    msd = 0.0
    rdf = 0.0
    spd = 0.0
    pepp = 0.0
    sigma = sqrt(2.0_rk*temp*gamma_/dt)
    ! Thermostating
    do i=1, N
      force(i,:) = -gamma_*vel(i,:) + box_muller(0.0_rk, sigma)
    end do
  end subroutine init_step
  
  
  subroutine add_force(i, j, r, pos1, pos2, force)
    ! Adds the force between two particles at the
    ! positions pos1 and pos2 to their total forces
    integer(ik), intent(in) :: i, j
    real(rk), intent(in) :: r, pos1(3), pos2(3)
    real(rk), intent(inout) :: force(N,3)
    real(rk) :: f(3)
    
    f = 24*(2*(1/r**13)-(1/r**7))*(pos1-pos2)/r
    force(i,:) = force(i,:) + f
    force(j,:) = force(j,:) - f
  end subroutine add_force
  
  
  subroutine add_rdf(L, r, rdf)
  	! Adds the distance r to the radial distribution function
    real(rk), intent(in) :: L, r
    real(rk), intent(inout) :: rdf(rdf_bins)
    real(rk) :: dr
    integer(ik) :: idx
    
    dr = (rc-r_min)/rdf_bins
    idx = (r-r_min)/dr+1
    rdf(idx) = rdf(idx) + 2.0_rk/N/(4*PI*r**2*dr)/(N/L**3)
  end subroutine add_rdf
  
  
  subroutine add_spd(speed, spd)
  	! Adds the speed to the speed distribution
    real(rk), intent(in) :: speed
    real(rk), intent(inout) :: spd(spd_bins)
    real(rk) :: dv
    integer(ik) :: idx
    
    dv = sc/spd_bins
    if (speed < sc) then
      idx = speed/dv+1
      spd(idx) = spd(idx) + 1.0_rk/(dv*N)
    end if
  end subroutine add_spd
  
  
  function add_msd(L, init_pos, pos, img_idx)
  	! Adds the displacement from the initial position init_pos to the 
  	! position pos to the mean squared displacement of the whole system
    real(rk), intent(in) :: L, init_pos(3), pos(3), img_idx(3)
    real(rk) :: add_msd
    
    add_msd = 1.0_rk/N*dist((pos + L*img_idx), init_pos)**2
  end function add_msd
  
  
  function add_potential(r)
  	! Adds the potential energy of two particles separated by 
  	! the distance r to the total potential energy of the system
    real(rk), intent(in) :: r
    real(rk) :: add_potential
    
    add_potential = 4*((1/r)**12-(1/r)**6)
  end function add_potential
  
  
  function add_kinetic(speed)
  	! Adds the kinetic energy corresponding to the given 
  	! speed to the total kinetic energy of the system
    real(rk), intent(in) :: speed
    real(rk) :: add_kinetic
    
    add_kinetic = 0.5*(speed**2)
  end function add_kinetic
    
    
  subroutine update_pos(L, pos, vel, force, img)
  	! Updates the positions of the particles taking 
  	! into account the periodic boundary conditions
    real(rk), intent(in) :: L
    real(rk), intent(inout) :: pos(N,3), vel(N,3), force(N,3), img(N,3)
    integer(ik) :: i, j
    
    do i=1, N
      pos(i,:) = pos(i,:) + vel(i,:)*dt + 0.5*force(i,:)*dt**2
      do j=1, 3
        if (2*pos(i,j) .ge. L) then
          pos(i,j) = pos(i,j) - L
          img(i,j) = img(i,j) + 1
        else if (2*pos(i,j) .le. -L) then
          pos(i,j) = pos(i,j) + L
          img(i,j) = img(i,j) - 1
        end if
      end do
    end do   
  end subroutine update_pos
  
  
  subroutine update(L, pos, vel, force, PE, KE, inst_temp, msd, init_pos, img, rdf, spd, pepp)
  	! Updates the rest of the variables
    real(rk), intent(in) :: L
    real(rk), intent(inout) :: PE, KE, inst_temp, msd
    real(rk), intent(inout) :: pos(N,3), vel(N,3), force(N,3), init_pos(N,3), img(N,3), rdf(rdf_bins), spd(spd_bins), pepp(N)
    real(rk) :: pos1(3), pos2(3), speed, r
    integer(ik) :: i, j, k
    
    do i=1, N
      speed = sqrt(sum(vel(i,:)**2))
      call add_spd(speed, spd)
      msd = msd + add_msd(L, init_pos(i,:), pos(i,:), img(i,:))
      KE = KE + add_kinetic(speed)
      do j=1, N
        if (i<j) then ! counts each particle only once
          pos1 = pos(i,:)
          pos2 = eff_pos(pos(i,:), pos(j,:), L)
          r = dist(pos1, pos2)
          if (r < rc) then
            PE = PE + add_potential(r)
            pepp(i) = pepp(i) + add_potential(r)
            call add_force(i, j, r, pos1, pos2, force)
            call add_rdf(L, r, rdf)
          end if
        end if
      end do
    end do
    inst_temp = 2*KE/(3*N-3)
  end subroutine update
  
  
  subroutine simulate(L, temp, pos, vel, force, PE, KE, inst_temp, msd, init_pos, img, rdf, spd, pepp)
  	! Simulates the system forward by one step
    real(rk), intent(in) :: L, temp
    real(rk), intent(inout) :: PE, KE, inst_temp, msd
    real(rk), intent(inout) :: pos(N,3), vel(N,3), force(N,3), init_pos(N,3), img(N,3), rdf(rdf_bins), spd(spd_bins), pepp(N)
    real(rk) :: prev_force(N,3)
    
    prev_force = force
    call update_pos(L, pos, vel, force, img)
    call init_step(L, temp, inst_temp, vel, force, msd, PE, KE, rdf, spd, pepp)
    call update(L, pos, vel, force, PE, KE, inst_temp, msd, init_pos, img, rdf, spd,  pepp)
    vel = vel + 0.5*(force + prev_force)*dt
  end subroutine simulate
  
end module simulation
