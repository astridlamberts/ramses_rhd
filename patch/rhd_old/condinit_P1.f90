! Initial conditions for a relativistic shock tube
! Taken from Ryu et al, 2006
subroutine condinit(mype)
  use variables
  use hydro_parameters
  implicit none

  integer ::i,j,k,l
  integer::mype
  real(dp)::vy,vz,vx,d,p,h,lor

  do l=1,ngrid

     do k=1,nz
        do j=1,ny
           do i=-1,nx+2
              if ( (x(i))<0.5) then
                 vx=0.0d0
                 vy=0.0d0
                 vz=0.0d0
                 d =10.0d0
                 p =13.3d0
                 h=1.0d0+gamma/(gamma-1.0d0)*p/d
                 lor=(1.d0-(vx**2+vy**2+vz**2))**(-1.d0/2.d0)
                 uin(l,i,j,k,1) = d*lor
                 uin(l,i,j,k,2) = lor**2*d*h*vx
                 uin(l,i,j,k,3) = lor**2*d*h*vy
                 uin(l,i,j,k,4) = lor**2*d*h*vz
                 uin(l,i,j,k,5) = lor**2*d*h-p!p/(gamma-1.d0) + uin(l,i,j,k,1) 
              else
                 vx=0.0d0
                 vy=0.0d0
                 vz=0.0d0
                 d =1.0d0
                 p =1.0e-6
                 h=1.0d0+gamma/(gamma-1.0d0)*p/d
                 lor=(1.-(vx**2+vy**2+vz**2))**(-1.d0/2.d0)
                 uin(l,i,j,k,1) = d*lor
                 uin(l,i,j,k,2) = lor**2*d*h*vx
                 uin(l,i,j,k,3) = lor**2*d*h*vy
                 uin(l,i,j,k,4) = lor**2*d*h*vz
                 uin(l,i,j,k,5) = lor**2*d*h-p! p/(gamma-1.d0) + uin(l,i,j,k,1) 
              endif
              !bx
              uin(l,i,j,k,6) = 0.d0
              
              !by
              uin(l,i,j,k,7) = 0.d0

              !bz
              uin(l,i,j,k,8) = 0.d0              

           end do
        end do
     end do

  end do

end subroutine condinit
