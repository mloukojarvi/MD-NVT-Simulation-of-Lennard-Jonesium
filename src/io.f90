module io
  ! Methods for handling input and output
  use const
  implicit none
    
contains 

  subroutine read_params(N, temp, dens, tot_steps, dt, gamma_)
    ! Reads the simulation parameters from the command line
    integer(ik), intent(out) :: N, tot_steps
    real(rk), intent(out) :: temp, dens, dt, gamma_
    integer(ik) :: iarg
    character(len=80) :: arg
    
    iarg=command_argument_count()
    if (iarg /= 6) then
      call get_command_argument(0,arg)
      write(0,'(a,a,a)') 'usage: ',trim(arg), &
      ' number of particles, target temperature, density, total simualtion steps, time step, gamma coefficient'
      stop
    end if
    
    call get_command_argument(1,arg)
    read(arg,*) N
    call get_command_argument(2,arg)
    read(arg,*) temp
    call get_command_argument(3,arg)
    read(arg,*) dens
    call get_command_argument(4,arg)
    read(arg,*) tot_steps
    call get_command_argument(5,arg)
    read(arg,*) dt
    call get_command_argument(6,arg)
    read(arg,*) gamma_ 
  end subroutine read_params
  
  
  subroutine init_files(dens, temp, files)
    ! Initializes the files that are used to store the simulation data
    real(rk), intent(in) :: dens, temp
    character(len=29), intent(inout) :: files(5)
    integer(ik) :: i, ios
    character(len=5) :: str_dens, str_temp
    character(len=26) :: id
    logical :: exists
    
    ! Creates a unique filename for each run
    write(str_dens , '(f5.3)') dens
    write(str_temp , '(f5.3)') temp
    id = '_dens=' // str_dens // '_temp=' // str_temp // '.dat'
    files = ['mea' // id, 'pos' // id, 'spd' // id, 'rdf' // id, 'ppe' // id]
    
    do i=1, 5
      inquire(file=files(i), exist=exists)
      if (exists) then
        open(unit=1, file=files(i), iostat=ios, status='replace')
      else
        open(unit=1, file=files(i), iostat=ios, status='new')
      end if
      if (ios/=0) then
        print '(a,a)', 'Error in creating file ', trim(files(i))
        stop
      end if
      close(1)
    end do
    
    do i=1, 5
      call to_file(files(i), 9_ik, [dfloat(N), temp, dens, dfloat(tot_steps), dt, dfloat(spd_bins), dfloat(rdf_bins), sc, rc])
    end do
  end subroutine init_files
  
  
  subroutine to_file(filename, datapoints, data_to_file)
    ! Writes the given data into the given file
    integer(ik), intent(in) :: datapoints
    real(rk), intent(in) :: data_to_file(datapoints)
    character(len=*), intent(in) :: filename
    integer(ik) :: ios
    
    open(unit=1, file=filename, iostat=ios, status="old", position="append", action="write")
    if (ios/=0) then
      print '(a,a)', 'Error in opening file ', trim(filename)
      stop
    end if
    write(1,*,iostat=ios) data_to_file
    if (ios/=0) then
      print '(a,a)', 'Error in writing data to file ', trim(filename)
      stop
    end if
    close(1)
  end subroutine to_file
  
  
  subroutine write_data(files, step, time, temp, PE, KE, inst_temp, pos, spd, rdf, msd, pepp)
    ! Writes the simulation data into files
    character(len=29), intent(in) :: files(5)
    integer(ik), intent(in) :: step
    real(rk), intent(in) :: time, temp, PE, KE, inst_temp, msd, pos(N,3), spd(spd_bins), rdf(rdf_bins), pepp(N)
    
    call to_file(files(1), 5_ik, [time, inst_temp, PE, KE, msd])
    call to_file(files(4), rdf_bins, rdf)
    call to_file(files(3), spd_bins, spd)
    
    if (mod(step,50) == 0 .or. step == 1) then
      call to_file(files(2), 3*N+1, [time, transpose(pos)])
      call to_file(files(5), N, pepp)
    end if
  end subroutine write_data
  
  
  subroutine info(step, time, temp, dens, PE, KE, inst_temp)
    ! Prints data about the current execution
    integer(ik), intent(in) :: step
    real(rk), intent(in) :: time, temp, dens, PE, KE, inst_temp
    
    if (step == 1) then
      print '(a)', '--------------------------------------------------------------------------------------------------------------'
      print '(a,i3,a,a,i5,a,a,f5.3,a,a,f4.2,a,a,f4.2,a)', '|Number of Particles: ', N, '| ', '|Total Steps: ', tot_steps, '| ', &
                                                          '|Time Step: ', dt, '| ', '|Target Temperature: ', temp, '| ', &
                                                          '|Density: ', dens, '|'
      print '(a)', '--------------------------------------------------------------------------------------------------------------'
    end if
    if (mod(step,1000) == 0 .or. step == 1) then
      print '(a,f5.1,a,a,f4.2,a,a,f8.2,a,a,f8.2,a)', '|Time: ', time, '| ', '|Current temperature: ', inst_temp, '| ', &
                                                     '|Potential Energy: ', PE, '| ', '|Kinetic Energy: ', KE, '|'
      print *, ' '
    end if
  end subroutine info

end module io
