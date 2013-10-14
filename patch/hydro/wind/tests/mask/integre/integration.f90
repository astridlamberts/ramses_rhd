! Program which computes the mean value of the hydro variables in the mask zone using the midpoint method
! check for 1D,2D,3D!!!!!!
! divide by total volume!!!!!!!!!!!



program integration3D


! just for 1 cell for the moment

integer ::n=10 ! number of slices in a cell, in a direction

integer :: i,j,k
!real(kind=8)::x,y,z
real(kind=8)::density

real(kind=8)::integral=0.0d0
real(kind=8)::dx=0.1d0,dy=0.1d0,dz=0.1d0

real(kind=8)::x0=0.d0,y0=0.d0,z0=0.d0

do i=1,n
   do j=1,n
      do k=1,n
! computing the integral in the middle of the slice
         integral =integral +dx*dy*dz*density3D(x0+(i-5-0.5d0)*dx,y0+(j-5-0.5d0)*dy,z0+(k-5-0.5d0)*dz)
      end do   
    end do
end do  

write(*,*),integral

end program integration3D

subroutine integration2D(xx,yy)


! just for 1 cell for the moment

integer ::n=10 ! number of slices in a cell, in a direction

integer :: i,j,k
!real(kind=8)::x,y,z
real(kind=8)::density

real(kind=8)::integral=0.0d0
real(kind=8)::dx=0.1d0,dy=0.1d0! siez of a slice

real(kind=8)::xx,yy   ! position of the center of the cell 

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice
         integral =integral +dx*dy*density2D(xx+(i-5-0.5d0)*dx,yy+(j-5-0.5d0)*dy) 
    end do
end do  

write(*,*),integral

end subroutine integration2D






function density3D(x,y,z) ! function which returns the density for a given position
real(kind=8)::x,y,z
real(kind=8)::rho_0 != 540.0d0 
real(kind=8)::r0    = 0.1d0
real(kind=8)::r ! distance from center of the mask

!rho_0=540.0d0
rho_0=8.0d0
!density=a*x*y*z

r=sqrt(x**2+y**2+z**2)
density=rho_0*(r0/r)**2
!density=rho_0*r

return
end function density3D


function density2D(x,y) ! function which returns the density for a given position
real(kind=8)::x,y
real(kind=8)::rho_0 != 540.0d0 
real(kind=8)::r0    = 0.1d0
real(kind=8)::r ! distance from center of the mask

!rho_0=540.0d0
rho_0=8.0d0
!density=a*x*y*z

r=sqrt(x**2+y**2)
density=rho_0*(r0/r)**2
!density=rho_0*r

return
end function density2D

