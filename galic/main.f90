program main
  use gadget_tools

  real(kind=4)::lbox,hubble,omega0,omegaL,aexp,mass_in_sol
  integer::npart,i


  call read_gadget(0,'toto.out',-1,.true.,.false.,lbox,hubble,omega0,omegaL,aexp,mass_in_sol,npart)

  write(*,*)npart
  do i=0,10
     write(*,*)pos(1:3,i)
  end do

end program main
