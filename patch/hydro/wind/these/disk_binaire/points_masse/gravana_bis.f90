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
  real(dp)::mtot,gmtot,estar

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
     gmass=GM       !gravity_params(1) ! GM
     emass=rmask1 ! Softening length
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
 ! for 2d non excentric orbit. Will need a leapfrog
  if(gravity_type==3)then 
     Mstar=12.5
     Mtot=12.5+1.4
     GMtot=484+0.1*484
     

     !mstar=488.                 ! GM Be star
     !rpuls=1.               ! Pulsar distance
    ! rmin=gravity_params(1)   ! Hard sphere radius
     qpuls=1.4/Mstar!10.!gravity_params(4)   ! Pulsar to star mass ratio
     epuls=0.05!gravity_params(5)   ! Pulsar gravity softening
     !omega=1.0/rpuls**(1.5)
!     omega2=1.0/rmin**(1.5)
!     omega0=omega1*omega2/(omega1+omega2)
     !phip=omega*t
     !mass center
     xc=0.5*boxlen
     yc=0.5*boxlen
     zc=0.5*boxlen
     estar=0.05

     !On appele le leapfrog pour deplacer les points masses

     call leapfrog(x01,y01,z01,vx01,vy01,vz01,x02,y02,z02,vx02,vy02,vz02,dtnew(nlevelmax),boxlen)
     

     !xpuls=rpuls*(1.0-qpuls)*cos(phip)+xc
     !ypuls=rpuls*(1.0-qpuls)*sin(phip)+yc
     !zpuls=0.0                      +zc
     !xstar=-rpuls*qpuls*cos(phip)     +xc
     !ystar=-rpuls*qpuls*sin(phip)     +yc
     !zstar= 0.0                     +zc
     do i=1,ncell
        !distance for a given cell to both bodies
        rxp=0.0d0; ryp=0.0d0; rzp=0.0d0 ! Pulsar
        rxs=0.0d0; rys=0.0d0; rzs=0.0d0 ! Star
        if(ndim>0)rxp=x(i,1)-x02
        if(ndim>1)ryp=x(i,2)-y02
        if(ndim>2)rzp=x(i,3)-z02
        if(ndim>0)rxs=x(i,1)-x01
        if(ndim>1)rys=x(i,2)-y01
        if(ndim>2)rzs=x(i,3)-z01
        rrpuls=sqrt(rxp**2+ryp**2+rzp**2+epuls**2)
        rrstar=sqrt(rxs**2+rys**2+rzs**2+estar**2)
        !force due to the star
        ff=sqrt(mstar*(1.0-qpuls)/rrstar**3)
!        ff2=sqrt(msun*(1.0-qpla)/rmin**3)
!        ffr=ff1*ff2/(ff1+ff2)
 !       ffr=ffr**2
        if(ndim>0)f(i,1)=-rxs*ff
        if(ndim>1)f(i,2)=-rys*ff
        if(ndim>2)f(i,3)=-rzs*ff
        !force due to the pulsar
        if(ndim>0)f(i,1)=f(i,1)-mstar*qpuls*rxp/rrpuls**3
        if(ndim>1)f(i,2)=f(i,2)-mstar*qpuls*ryp/rrpuls**3
        if(ndim>2)f(i,3)=f(i,3)-mstar*qpuls*rzp/rrpuls**3
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


!###########################################################
!###########################################################
!###########################################################
!###########################################################

!Subroutine which computes the positions of the bodies at each timestep
subroutine leapfrog(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,dt,boxlen)
implicit none

real(kind=8)::x1,y1,z1,vx1,vy1,vz1,xf1,yf1,zf1,vxf1,vyf1,vzf1
real(kind=8)::x2,y2,z2,vx2,vy2,vz2,xf2,yf2,zf2,vxf2,vyf2,vzf2
real(kind=8)::ax1,ay1,az1,ax2,ay2,az2
real(kind=8)::dt,boxlen
real(kind=8)::r1,r2,r
real(kind=8)::GMstar,GMpuls
real(kind=8)::xc,yc,zc

GMstar=484
GMpuls=48
xc=0.5*boxlen
yc=0.5*boxlen
zc=0.5*boxlen

!write(*,*),x,y,vx,vy,det,'leapfrog'
!moving the star
   xf1=x1+vx1*dt
   yf1=y1+vy1*dt
   zf1=z1+vz1*dt

   xf2=x2+vx2*dt
   yf2=y2+vy2*dt
   zf2=z2+vz2*dt

   r=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

   ax1=-GMpuls/(r**2)*(xf1-xf2)/r
   ay1=-GMpuls/(r**2)*(yf1-yf2)/r
   az1=-GMpuls/(r**2)*(zf1-zf2)/r

   vxf1=vx1+ax1*dt
   vyf1=vy1+ay1*dt
   vzf1=vz1+az1*dt
 
   x1=xf1
   y1=yf1
   z1=zf1

   vx1=vxf1
   vy1=vyf1
   vz1=vzf1

!moving the pulsar


!   r2=sqrt((xf2-xc)**2+(yf2-yc)**2+(zf2-zc)**2)

   ax2=-GMstar/(r**2)*(xf2-xf1)/r
   ay2=-GMstar/(r**2)*(yf2-yf1)/r
   az2=-GMstar/(r**2)*(zf2-zf1)/r

   vxf2=vx2+ax2*dt
   vyf2=vy2+ay2*dt
   vzf2=vz2+az2*dt
 
   x2=xf2
   y2=yf2
   z2=zf2

   vx2=vxf2
   vy2=vyf2
   vz2=vzf2

end subroutine leapfrog
