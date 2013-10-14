!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters  
  use wind_parameters
  use amr_commons,only:dtnew
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rrpuls,rrstar,ff,rx,ry,rz,rr
  real(dp)::mstar,rpuls,qpuls,epuls,omega,phip,xc,yc,zc
  real(dp)::xpuls,ypuls,zpuls,xstar,ystar,zstar,rxp,ryp,rzp,rxs,rys,rzs
  real(dp)::mtot,estar

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     gmass=Gstar*M1       !gravity_params(1) ! GM
     emass=r01 ! Softening length
     xmass=x01!gravity_params(3) ! Point mass coordinates
     ymass=y01!gravity_params(4)
     zmass=z01!gravity_params(5)

     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)

        f(i,1)=-gmass*rx/rr**3
#if NDIM>1
        f(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rr**3
#endif
    end do
  end if





  ! Pulsar orbiting around a star. Inspired from patch/planets

  if(gravity_type==3)then 
     
     ! if there is no second object, it has no gravitational impact
     if (type2 .ne. 'puls') then
        M2=0.0d0
     endif

     do i=1,ncell
        rxp=0.0d0; ryp=0.0d0; rzp=0.0d0 ! Pulsar
        rxs=0.0d0; rys=0.0d0; rzs=0.0d0 ! Star
        rxs=x(i,1)-x01
        rxp=x(i,1)-x02
#if NDIM>1
        rys=x(i,2)-y01
        ryp=x(i,2)-y02
#endif
#if NDIM>2
        rzs=x(i,3)-z01
        rzp=x(i,3)-z02
#endif
        rrstar=sqrt(rxs**2+rys**2+rzs**2+r01**2)
        rrpuls=sqrt(rxp**2+ryp**2+rzp**2+r02**2)
        f(i,1)=-Gstar*m1*rxs/rrstar**3
        f(i,1)=f(i,1)-Gstar*m2*rxp/rrpuls**3
#if NDIM>1
        f(i,2)=-Gstar*m1*rys/rrstar**3
        f(i,2)=f(i,2)-Gstar*m2*ryp/rrpuls**3
#endif
#if NDIM>2
        f(i,3)=-Gstar*m1*rzs/rrstar**3
        f(i,3)=f(i,3)-Gstar*m2*rzp/rrpuls**3
#endif
    end do
  end if



end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana


