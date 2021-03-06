!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
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
 
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  integer:: i,id,iu,iv,iw,ip
  real(dp):: x01,y01,z01,rc,xx,yy,zz,r01
  real(dp):: rho_amb,P_amb
  real(dp):: R_star1,P_01,rho_01
  integer :: ivar
  real(dp):: eta1, nu1, alpha1
  real(dp)::G,M
  real(dp)::beta
 
  ! Add here, if you wish, some user-defined initial conditions
 

  IF (nn .eq. 0) RETURN


  id=1; iu=2; iv=3; iw=4; ip=ndim+2

!run parameters
 
  R_star1 = 0.01

  beta = 1.0


! a voir, pour le milieu ambiant j'utilise la premiere etoile comme reference!!!!!!!!!
 ! rho_amb = 0.1*eta1*alpha1**(-2)*nu1**(-1)
 ! P_amb = 45.d0*nu1*eta1*alpha1**(-2)

  rho_amb= 0.1d0
  P_amb=26.4d0

 ! voir les dimensions!!!!!!!!

  !P_01 = 361*nu1*eta1*alpha1**(-2)
  !rho_01 = 540*eta1*alpha1**(-2)*nu1**(-1)


! a voir!!!!!!!!
  G=1.0d0
  M=1.0d0 

!  write(*,*),P_01,P_02,'P1,2'
 ! There is an initial domain filled with a densité rho_amb and sound speed cs_amb. Afterwards, we add a mask (for r <2r0) with a given densité, pressure and speed. 



q(:,id) = rho_amb
q(:,iu) = 0.0
q(:,iv) = 0.0
q(:,ip) = P_amb



!Disk

  DO i=1,nn ! all the cells
     xx=x(i,1)-x01 ! position cell/center
     yy=x(i,2)-y01
     rc=sqrt(xx**2+yy**2)! 2d distance /center


#if NDIM == 3
     zz=x(i,3)-z01
     rc3=sqrt(xx**2+yy**2+zz*2)
     IF (rc3 .lt. 2*r01) THEN
        q(i,id) = min(rho_01*(r01/rc),4.d0*rho_0)
        q(i,iu) = v_inf*(xx/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,iv) = v_inf*(yy/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,iw) = v_inf*(zz/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,ip) = min(P_01*(r01/rc)**(10/3),2.d0**(2.d0*gamma)*P_0)
      ENDIF
#endif

#if NDIM==2   
        !write(*,*)'on est a 2D'

      !!!! Voir condition min et max!!!!
        IF (rc .lt. 2*r01) THEN
        q(i,id) = max(rho_01*(r01/R_star1)**(-3),rho_01)
! a priori pas besoin de la condition!!!!
        q(i,iu) = min(xx/rc*sqrt(G*M/rc),xx/R_star1*sqrt(G*M/R_star1))
        q(i,iv) = min(yy/rc*sqrt(G*M/rc),yy/R_star1*sqrt(G*M/R_star1))
        q(i,ip) = min(P_01*(r01/rc)**(10/3),2.d0**(2.d0*gamma)*P_01)
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
