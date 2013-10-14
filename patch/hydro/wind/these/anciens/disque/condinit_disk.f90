!very first creation of the disk, without mask


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
  real(dp):: rho_amb,P_amb,rho_0,r_0,p,qu,H_0,H,cs,cs_0,p_0,v
  real(dp)::G,M
  integer::  ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  IF (nn .eq. 0) RETURN


! voir les unites
  id=1; iu=2; iv=3; iw=4; ip=ndim+2
  rho_amb=0.01d0
  P_amb=0.0001d0
  rho_0=1.0d0
  r_0=1.0d0
  p=2.0d0
  qu=1.0d0
  H_0=0.1d0
  cs_0=0.1d0
  p_0=cs_0**2*rho_0
  G=1.0d0
  M=1.0d0
  
!position of the center of the disk (in user units, not necessarily [0;1])
  x0=5.0d0
  y0=5.0d0
  z0=5.039d0
  

  q(:,id) = rho_amb
  q(:,iu) = 0.0
  q(:,iv) = 0.0
  q(:,iw) = 0.0
  q(:,ip) = P_amb

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     zz=x(i,3)-z0

     rc=sqrt(xx**2+yy**2)
     H=H_0*(rc/r_0)**((3.-qu)/2.)
#if NDIM==2 
     zz=0.0d0
#endif

     q(i,id)=rho_0*(r_0/rc)**p*exp(-zz**2/2.d0/H**2)
     cs=cs_0*(r_0/rc)**(qu/2.0d0)
     q(i,ip)=cs**2*q(i,id)
     v=sqrt(G*M/rc)
     q(i,iu)=v*yy/rc
     q(i,iv)=-v*xx/rc
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
