!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i
  real(dp)::T1,H,rho0

!!$  do ivar=1,nvar
!!$     do i=1,ncell
!!$        u(i,ivar)=boundary_var(ibound,ivar)
!!$     end do
!!$  end do
  
  ! Add here, if you wish, some user-defined boudary conditions
  ! ........
  T1=t_region(1)/0.5
  H=T1/abs(gravity_params(2))
 rho0=d_region(1)
  DO i=1,ncell
     u(i,1)=rho0*EXP(-x(i,2)/H)
     u(i,ndim+3)=t_region(1)*u(i,1)
     u(i,2)=0.0
     u(i,3)=0.0
     u(i,ndim+2)=u(i,1)*T1/(gamma-1.0)     


 ENDDO
end subroutine boundana
