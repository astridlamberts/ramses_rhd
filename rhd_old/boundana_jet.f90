!############################# 
!Boundary conditions for a 2D jet coming from the left
! the box is 0<x<30 by -15<y<15. Initially the jet 
! is restricted to -1< y<1 with inflow boundary conditions.
!The rest of the left boundary is zero gradient

subroutine xinner_ana
  use hydro_parameters
  use variables
  use amr_parameters
  implicit none
!
  integer :: j,k
  real(dp)::r,ncell,rmin,rmax,vx,vy,vz,lor,d,h,p,boxlen
  boxlen=30.
  r=1.

  rmin=boxlen/2.-r
  rmax=boxlen/2.+r
  
  do k=ku1,ku2
     do j=ju1,ju2
!        if ((j >(rmin *ny*1./boxlen)) .and. (j<(rmax*ny*1./boxlen))) then
        if (abs(y(j))<1.) then
           !inflowing jet
           vx=0.99d0
           vy=0.0d0
           vz=0.0d0
           d =.1d0
           p =0.01d0
           h=1.0d0+gamma/(gamma-1.0d0)*p/d
           lor=(1.-(vx**2+vy**2+vz**2))**(-1.d0/2.d0)

           uin(1,iu1:iu1+2,j,k,1) = d*lor
           uin(1,iu1:iu1+2,j,k,2) = lor**2*d*h*vx
           uin(1,iu1:iu1+2,j,k,3) = lor**2*d*h*vy
           uin(1,iu1:iu1+2,j,k,4) = lor**2*d*h*vz
           uin(1,iu1:iu1+2,j,k,5) = lor**2*d*h-p

        else

           ! Zero-gradient BC
           uin(1,iu1+2,j,k,:) = uin(1,1,j,k,:)
           uin(1,iu1+1,j,k,:) = uin(1,1,j,k,:)
           uin(1,iu1  ,j,k,:) = uin(1,1,j,k,:)

        endif
     end do
  end do
!
  return
end subroutine xinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine xouter_ana
  use hydro_parameters
  use variables
  implicit none

  integer ::k,j

  ! Periodic BC
  do k=ku1,ku2
     do j=ju1,ju2
        uin(1,iu2-2,j,k,:) = uin(1,1 ,j,k,:)
        uin(1,iu2-1,j,k,:) = uin(1,2 ,j,k,:)
        uin(1,iu2  ,j,k,:) = uin(1,3 ,j,k,:)
     end do
  end do
!
  return
end subroutine xouter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine yinner_ana
#if NDIM>1
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i

  ! Periodic BC
  do k=ku1,ku2
     do i=iu1,iu2
        uin(1,i,ju1  ,k,:) = uin(1,i,ny-2,k,:)
        uin(1,i,ju1+1,k,:) = uin(1,i,ny-1,k,:)
        uin(1,i,ju1+2,k,:) = uin(1,i,ny  ,k,:)
     end do
  end do
!
#endif
  return
end subroutine yinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine youter_ana
#if NDIM>1
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i

  do k=ku1,ku2
     do i=iu1,iu2
        uin(1,i,ju2-2,k,:) = uin(1,i,1 ,k,:)
        uin(1,i,ju2-1,k,:) = uin(1,i,2 ,k,:)
        uin(1,i,ju2  ,k,:) = uin(1,i,3 ,k,:)
     end do
  end do
!
#endif
  return
end subroutine youter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zinner_ana
#if NDIM > 2
  use hydro_parameters
  use variables
  implicit none

  integer ::j,i

  ! Periodic BC
  do j=ju1,ju2
     do i=iu1,iu2
        uin(1,i,j,ku1  ,:) = uin(1,i,j,nz-2,:)
        uin(1,i,j,ku1+1,:) = uin(1,i,j,nz-1,:)
        uin(1,i,j,ku1+2,:) = uin(1,i,j,nz  ,:)
     end do
  end do
!
#endif 
  return
end subroutine zinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zouter_ana
#if NDIM > 2
  use hydro_parameters
  use variables
  implicit none

  integer ::j,i

  ! Periodic BC
  do j=ju1,ju2
     do i=iu1,iu2
        uin(1,i,j,ku2-2,:) = uin(1,i,j,1 ,:)
        uin(1,i,j,ku2-1,:) = uin(1,i,j,2 ,:)
        uin(1,i,j,ku2  ,:) = uin(1,i,j,3 ,:)
     end do
  end do
!
#endif
  return
end subroutine zouter_ana
