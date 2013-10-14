program ramses
  implicit none  
  real::t1,t2!!!
  integer::i
  ! Read run parameters
  call read_params

  i=0!
  ! Start time integration
  call cpu_time(t1) !
  call adaptive_loop
  call cpu_time(t2)!
  print *, 'temps ecoule',t2-t1,'temps=',i
  i=i+1


end program ramses

