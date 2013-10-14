!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use const
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer:: i,j,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0,rcut
  real(dp):: R_star, beta,rho_amb,P_amb,P0,rho0
  integer::  ivar
  real(dp):: eta,nu, alpha, v_w,v_inf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  IF (nn .eq. 0) RETURN


  id=1; iu=2; iv=3; iw=4; ip=ndim+2

!run parameters
  eta = 0.065
   nu = 0.85
   alpha = 0.69


  x0   = 0.5*boxlen
  y0   = 0.5*boxlen
  r0   = 0.195*alpha
!  write(*,*),r0
  z0   = 0.5*boxlen
  R_star = 0.0539 !(if alpha=0.69)
  beta = 1.0
  rho_amb = 0.1*eta*alpha**(-2)*nu**(-1)
   P_amb = 45.d0*nu*eta*alpha**(-2)
  P0 = 361*nu*eta*alpha**(-2)
  rho0 = 540*eta*alpha**(-2)*nu**(-1)
  v_w=nu
  v_inf = 1.36*v_w

 
! There is an initial domain filled with a densité rho_amb and sound speed cs_amb. Afterwards, we add a mask (for r <2r0) with a given densité, pressure and speed. 



q(:,id) = rho_amb
q(:,iu) = 0.0
q(:,iv) = 0.0
q(:,ip) = P_amb


  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center
#if NDIM == 3
     zz=x(i,3)-z0
     rc3=sqrt(xx**2+yy**2+zz*2)
#endif
    

#if NDIM==3
       IF (rc3 .lt. 2*r0) THEN
        q(i,id) = rho0*(r0/rc)**2
        q(i,iu) = v_inf*(xx/rc3)*(1-R_star/rc3)**beta
        q(i,iv) = v_inf*(yy/rc3)*(1-R_star/rc3)**beta
        q(i,iw) = v_inf*(zz/rc3)*(1-R_star/rc3)**beta
        q(i,ip) = P0*(r0/rc)**(10/3)
        ENDIF
#endif

#if NDIM==2   
        !write(*,*)'on est a 2D'
        IF (rc .lt. 2*r0) THEN
           q(i,id) = min(rho0*(r0/rc)**2,4.d0*rho0)
           q(i,iu) = v_inf*(xx/rc)*max((1-R_star/rc)**beta,0.d0)
           q(i,iv) = v_inf*(yy/rc)*max((1-R_star/rc)**beta,0.d0)
           q(i,ip) = min(P0*(r0/rc)**(10/3),2.d0**(2.d0*gamma)*P0)
        ENDIF
#endif
   ENDDO     

     
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
