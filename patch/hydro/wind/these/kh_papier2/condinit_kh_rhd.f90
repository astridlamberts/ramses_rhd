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
  real(dp)::rho1,rho2,u2,r,u1,lor1,lor2
  real(dp)::h,lor,Mach1,Mach2,p1,p2
  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  alpha=0.99 

  y0=y_center(1)

  rho1=d_region(1)
  rho2=d_region(2)
  
  u1=u_region(1)
  u2=u_region(2)

  lor1=(1.-u1**2)**(-1./2.)
  lor2=(1.-u2**2)**(-1./2.)

!!!!
  Mach1=p_region(1)
  Mach2=p_region(2)

  
  p1=u1**2*lor1**2*rho1/(gamma*(mach1**2+u1**2*lor1**2-u1**2*lor1**2/(gamma-1.)))
  p2=u2**2*lor2**2*rho1/(gamma*(mach2**2+u1**2*lor2**2-u2**2*lor2**2/(gamma-1.)))


CALL init_random_seed()        
            
!appel de la procÃ©dure qui gÃ©nÃ¨re les nombres aleatoires

  do i=1,nn
     if(x(i,2) <0.5) then
        q(i,id)=rho1!1+alpha
        q(i,iu)=u_region(1)
        CALL RANDOM_NUMBER(r)!*0.01-5e3
        q(i,iv)=r*0.01-5e-3
        q(i,ip)=p1
        q(i,5)=1.
        q(i,6)=0.
     else
        q(i,id)=rho2!1-alpha
        q(i,iu)=u_region(2)
        CALL RANDOM_NUMBER(r)!*0.01-5e3

        q(i,iv)=r*0.01-5e-3
        q(i,ip)=p2
        q(i,5)=0.
        q(i,6)=1.
     endif
  end do

! initialisation des perturbations de vy avec un bruit blanc
  

  do i=1,nn
     ! Convert primitive to conservative variables
     ! specific enthalpy
 
     h=1.0d0+gamma/(gamma-1.0d0)*q(i,5)/q(i,1)
     ! Lorentz factor
     lor=(1.-(q(i,2)**2+q(i,3)**2+q(i,4)**2))**(-1.d0/2.d0)

     ! proper density  -> density in lab frame
     u(i,1)= q(i,1)*lor
     ! proper 3-velocity -> momentum in lab frame
     u(i,2)=lor**2*q(i,1)*h*q(i,2)
     u(i,3)=lor**2*q(i,1)*h*q(i,3)
     u(i,4)=lor**2*q(i,1)*h*q(i,4)
     ! isotropic gas pressure -> total fluid energy
     u(i,5)=lor**2*q(i,1)*h-q(i,5)
    !passive scalars
     do ivar=6,nvar
        u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)*lor
     end do

  end do

end subroutine condinit


     

! exemple d'initialisation du generateur en utilisant l'heure
SUBROUTINE init_random_seed()
implicit none
INTEGER :: i, n, heure
INTEGER, DIMENSION(:), ALLOCATABLE :: graine

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
