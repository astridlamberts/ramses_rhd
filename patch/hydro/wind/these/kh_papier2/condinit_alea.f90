!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
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
  integer::i,j,id,iu,iv,iw,ip
  real(dp)::y0,p0,v0,v,alpha,x0
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::rho1,rho2,r


  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  alpha=0.99 

  y0=y_center(1)

  rho1=1.+alpha
  rho2=1.-alpha
  v=1.
  v0=0.1
  p0=1.



!initialisation de la graine
CALL init_random_seed()        
            
!appel de la procÃ©dure qui gÃ©nÃ¨re les nombres aleatoires

  do i=1,nn
     if(x(i,2) < y0)then
        q(i,id)=1+alpha
        q(i,iu)=v1
        CALL RANDOM_NUMBER(r)*0.01-5e3
        q(i,iv)=r
        q(i,ip)=p0
        q(i,5)=1.
        q(i,6)=0.
     else
        q(i,id)=1-alpha
        q(i,iu)=-v1
        CALL RANDOM_NUMBER(r)*0.01-5e3
        q(i,iv)=r
        q(i,ip)=p0
        q(i,5)=0.
        q(i,6)=1.
     endif
  end do

! initialisation des perturbations de vy avec un bruit blanc
  

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


     

! exemple d'initialisation du generateur en utilisant l'heure
SUBROUTINE init_random_seed()
implicit none
INTEGER :: i, n, heure
INTEGER, DIMENSION(:), ALLOCATABLE :: graine3

! demande au systeme quelle taille de graine il faut          
CALL RANDOM_SEED(size = n)
! on cree le tableau graine
ALLOCATE(graine(n))

! on regarde l'heure
CALL SYSTEM_CLOCK(COUNT=heure)
          
! on utilise l'heure pour initialiser la graine
graine = heure + 37 * (/ (i - 1, i = 1, n) /)
! (/ (i - 1, i = 1, n) /) veut dire : tableau, de cases 1 a n
! dans la case i, valeur i-1
! apres, on multiplie toutes les valeurs par 37, et on ajoute l'heure

! on met la graine dans le generateur
CALL RANDOM_SEED(PUT = graine)
          
! pour que la graine qu'on a utilisee ne prenne pas de place
DEALLOCATE(graine)
end subroutine
