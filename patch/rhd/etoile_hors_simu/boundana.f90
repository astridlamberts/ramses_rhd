!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use hydro_parameters, ONLY: nvar,boundary_var,gamma
  use wind_parameters,ONLY: Mdot1,vinf1,mach1,x01,y01,z01
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
  real(dp) :: r,theta,lor,rho_a,P_a,d,vx,vy,vz,P,scal1,scal2,h,a1,pi,r2d,phi
  


!  do ivar=1,nvar
!     do i=1,ncell
!        u(i,ivar)=boundary_var(ibound,ivar)
!     end do
!  end do
 
 
 
  ! Add here, if you wish, some user-defined boudary conditions
  ! ........



  ! calcul de la distance par rapport à l'étoile hors champ (positionnee en x1*boxlen,y1*boxlen)

  a1=1.d0
  pi=acos(-1.0d0)
  lor=(1.-vinf1**2)**(-1./2.)
  rho_a = Mdot1*a1/(2.0d0*pi*vinf1*lor)
  P_a=vinf1**2* rho_a/(gamma*Mach1**2-vinf1**2*gamma/(gamma-1))

  do i=1,ncell
#if NDIM==2    
     r=sqrt((x(i,1)-x01/boxlen)**2+(x(i,2)-y01/boxlen)**2)
     theta=atan((x(i,2)-y01/boxlen)/(x(i,1)-x01/boxlen))
     d=mdot1/(2*pi*r*vinf1*lor)
     vx=vinf1*cos(theta)
     vy=vinf1*sin(theta)
     vz=0
     P=P_a*(a1/r)**gamma
#endif

#if NDIM==3    
     r=sqrt((x(i,1)-x01/boxlen)**2+(x(i,2)-y01/boxlen)**2+(x(i,3)-z01/boxlen)**2)
     r2d=sqrt((x(i,1)-x01/boxlen)**2+(x(i,3)-z01/boxlen)**2)

    if ((x(i,1) -x01/boxlen) <0) then 
        theta=2.*pi-acos((x(i,3)-z01/boxlen)/r2d)
     else        
        theta=acos((x(i,3)-z01/boxlen)/r2d)
     endif
     phi=acos((x(i,2)-y01/boxlen)/r)
     vz=vinf1*cos(theta)*sin (phi)
     vx=vinf1*sin(theta)*sin(phi)
     vy=vinf1*cos(phi)
     d=mdot1/(4*pi*r**2**vinf1*lor)
     P=P_a*(a1/r)**gamma

#endif

     h=1.0d0+gamma/(gamma-1.0d0)*P/d
     ! proper density  -> density in lab frame
     u(i,1)= d*lor
     ! proper 3-velocity -> momentum in lab frame
     u(i,2)=lor**2*d*h*vx
     u(i,3)=lor**2*d*h*vy
     u(i,4)=lor**2*d*h*vz
     ! isotropic gas pressure -> total fluid energy
     u(i,5)=lor**2*d*h-P

     !passive scalars

  if (nvar .gt. 5) then
     scal1 =1.0d0 
     scal2 =0.0d0 
     u(i,6)=d*scal1*lor
     u(i,7)=d*scal2*lor
   endif


enddo
  

end subroutine boundana
