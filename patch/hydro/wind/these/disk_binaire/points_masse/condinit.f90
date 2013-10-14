 !================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  use wind_parameters
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
  integer :: ilevel
  integer:: i,j,id,iu,iv,iw,ip
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: v_a,P_a,rho_a
  integer :: ivar,ncell
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp)::a1=1.0d0 ! separation of the system in simulation units
  real(dp)::csound,omega,H,r3d
  real(dp)::GMtot,Mstar,Mtot,xc,yc,zc,r1,r2,r,gmpuls,gmstar
  real(dp)::inc,peri,noeud,theta,dthetadt




  ncell=ncoarse+twotondim*ngridmax

  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2


!à ameliorer!!!!!!!!!
  v_a = vinf2*(1.0d0-rstar2)**beta2 
  rho_a = Mdot2/(4.0d0*3.1415*a**2*v_a)

#if NDIM==2
 rho_a = Mdot2/(2.0d0*3.1415*a*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach2**2*gamma)

  q(:,id) = rho_amb
  q(:,iu) = 0.0d0
  q(:,iv) = 0.0d0
#if NDIM == 3
  q(:,iw) = 0.0d0
#endif
  q(:,ip) = P_amb




  xc=0.5d0*boxlen
  yc=0.5d0*boxlen
  zc=0.5d0*boxlen

!  Mstar=12.5
!  Mtot=12.5+1.4
  
!  gmstar=73.0!484.0d0
!  GMtot=80!84.0d0+0.1d0*484.0d0
!  gmpuls=7.!48.0d0

  
  !attention, mettre des radians
  inc=0.d0/180.d0*3.1415d0
  peri=0.d0/180.d0*3.1415d0
  noeud=0.d0
  theta=0.d0/180.d0*3.1415d0

  x1=xc+M2/(M1+M2)*a*(1-exc)*(cos(noeud)*cos(peri+theta)-sin(noeud)*sin(peri+theta)*cos(inc))
  y1=yc+M2/(M1+M2)*a*(1-exc)*(sin(noeud)*cos(peri+theta)+cos(noeud)*sin(peri+theta)*cos(inc))
  z1=zc+M2/(M1+M2)*a*(1-exc)*sin(peri+theta)*sin(inc)

  x2=x1-a*(1-exc)*(cos(noeud)*cos(peri+theta)-sin(noeud)*sin(peri+theta)*cos(inc))
  y2=y1-a*(1-exc)*(sin(noeud)*cos(peri+theta)+cos(noeud)*sin(peri+theta)*cos(inc))
  z2=z1-a*(1-exc)*sin(peri+theta)*sin(inc)
  
  r=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
!  write(*,*),r,xini1,xini2,yini1,yini2,zini1,zini2

  ax1=-G*M2/r**2*(x1-x2)/r
  ay1=-G*M2/r**2*(y1-y2)/r
  az1=-G*M2/r**2*(z1-z2)/r
  
  ax2=-G*M1/r**2*(x2-x1)/r
  ay2=-G*M1/r**2*(y2-y1)/r
  az2=-G*M1/r**2*(z2-z1)/r 
  



  dthetadt=sqrt(a*(1-exc**2)*G*(M1+M2)/(r**4))
  
  
  r1=sqrt((x1-xc)**2.0d0+(y1-yc)**2.0d0+(z1-zc)**2.0d0)
  r2=sqrt((x2-xc)**2.0d0+(y2-yc)**2.0d0+(z2-zc)**2.0d0)
  write(*,*),r,r1,r2
 
  if (direction .eq.'d') then 
  ! direct rotation
     
     vx1=dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r1
     vy1=dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r1
     vz1=dthetadt*cos(peri+theta)*sin(inc)*r1
     
     vx2=-dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r2
     vy2=-dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r2
     vz2=-dthetadt*cos(peri+theta)*sin(inc)*r2
  else if (direction .eq. 'i') then
     vx1=-dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r1
     vy1=-dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r1
     vz1=-dthetadt*cos(peri+theta)*sin(inc)*r1
     
     vx2=dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r2
     vy2=dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r2
     vz2=dthetadt*cos(peri+theta)*sin(inc)*r2
  else then 
     write(*,*), "please specify the rotation direction"
  endif


  !star
!  x01=xc-a1*(1.0d0-exc)*1.4d0/Mtot
!  y01=yc
!  z01=zc
  !positions of the 2nd body :pulsar
!  x02   = x01+a1*(1.0d0-exc)
!  y02   = y01
!  z02   = z01



!initial velocity field (necessary for the position update)
!  r=a1
!  r1=sqrt((x01-xc)**2+(y01-yc)**2+(z01-zc)**2)
!  vx01=0.0d0
!  vy01=-sqrt(Gmtot*(2.0D0/r-1.0d0/a1)) !r
!  vy01=-sqrt(Gmpuls*(2.0D0/r-1.0d0/a1)) !r
!  vy01=-sqrt((Gmtot)*(2.0D0/r-1.0d0/a1)) !r
!  vz01=0.0d0
!  r2= sqrt((x02-xc)**2+(y02-yc)**2+(z02-zc)**2)
!  vx02=0.0d0
!    vy02=sqrt(Gmtot*(2.0D0/r-1.0d0/a1)) !r
!  vy02= sqrt(GMstar*(2.0D0/r-1.0d0/a1))
! vy02= sqrt(GMtot*(2.0D0/r-1.0d0/a1))
!  vz02=0.0d0

!  call init_star(q,x,nn,vinf1,rstar1,Mdot1,Mach1,x01,y01,z01,rmask1,beta1)
  call init_puls(q,x,nn,vinf1,rmask1,Mdot1,Mach1,x01,y01,z01)

!  if (type2 .eq. 'puls')then
   call init_puls(q,x,nn,vinf2,rmask2,Mdot2,Mach2,x02,y02,z02)
 !endif
  

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




subroutine init_star(q,x,nn,vinf,rstar,Mdot,Mach,x0,y0,z0,rmask,beta)
  use wind_parameters
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
   real(dp)::rho_a,P_A,v_a  ! variables for r =a
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables!!!!
  real(dp):: int_rho,int_p,int_vx,int_vy,int_vz,vrotx,vroty ! function which computes the mean value of a variable in the mask
  real(dp):: frac
  real(dp)::pi=3.141459265
  integer:: i,j,id,iu,iv,iw,ip
  real::H,Mach,density,v,R3D,v_mask,beta,rho_0w
  real(dp)::vinf,rstar,Mdot,x0,y0,z0,vwind,vrot,rmask,omega,testx,testy,csound
 
 id=1 ; iu=2 ; iv=3; iw=4; ip=ndim+2

  dx=boxlen/2.0d0**(levelmin)

  !compute the density in the wind, at the edge of the mask
  !a=1.0d0 !(a ameliorer!!!)

  v_mask=vinf*(1.0d0-rstar/rmask)!**beta
  rho_0w=Mdot*a**2/(4.0d0*3.1415*(rmask)**2*v_mask)
!#if NDIM==2
! rho_0w = Mdot*a/(2.0d0*3.1415*(a*rmask)*v_mask)
!#endif

  zz = 0.0d0  

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center
     csound=cs_0*(r_0/rmask)**(qu/2.)

#if NDIM == 3 
     zz=x(i,3)-z0
     r3D=sqrt(xx**2+yy**2+zz**2)
     IF (r3d .lt. 1.2d0*rmask) THEN
        omega=0.0d0
        if (disk_presence .eqv. .true) then
          omega=sqrt(GM/(r3d*(rmask**2+r3d**2)))!! no pressure gradient correction yet
        endif  
     H=H_0*(rc/r_0)**((3.-qu)/2.)
     q(i,id)=   frac(xx,yy,zz,rmask,dx)* rho_0w +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id)    
     q(i,ip)=csound**2*q(i,id)/gamma
     vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
     q(i,iu) =  frac(xx,yy,zz,rmask,dx)* (vwind-omega*yy)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iu)
     vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
     q(i,iv) =  frac(xx,yy,zz,rmask,dx)* (vwind+omega*xx)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iv)
     q(i,iw) =  frac(xx,yy,zz,rmask,dx)* int_vz(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iw)
     ENDIF
#endif
     

#if NDIM==2   
     IF (rc .lt. 1.2d0*rmask) THEN
 ! In 2D there is either a disk, or a wind, there is no need having them both. 
! The density has to be continuous across the edge of the disk.
        if (disk_presence .eqv..true.) then
           q(i,id)= frac(xx,yy,zz,rmask,dx)* rho_0*(r_0/rmask)**pe +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id) 
           q(i,ip)=csound**2*q(i,id)/gamma
           omega=sqrt(GM/(rc*(rmask**2+rc**2))+csound**2*(pe+qu)/(rc**2)) 
           q(i,iu)=-omega*yy
           q(i,iv)=+omega*xx
           if ((q(i,iu)*xx+q(i,iv)*yy).gt. 1.0e-14) then
              write(*,*)q(i,iu)*xx+q(i,iv)*yy,'vitesse radiale star'
           endif
     else   
           q(i,id)=   frac(xx,yy,zz,rmask,dx)* rho_0w+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id) 
           q(i,ip)=csound**2*q(i,id)/gamma
           vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
           q(i,iu) =  frac(xx,yy,zz,rmask,dx)* vwind+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iu)
           vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
           q(i,iv) =  frac(xx,yy,zz,rmask,dx)*vwind+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iv)
       endif    
     ENDIF
#endif
  ENDDO     

end subroutine init_star


subroutine init_puls(q,x,nn,vinf,rmask,Mdot,Mach,x0,y0,z0) 
  use wind_parameters
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx 
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: vinf,rmask,Mdot,Mach,beta
  real(dp)::rho_a,P_A ! variables for r =a
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables
  real(dp):: int_rhoP,int_pP
  real(dp):: frac
  real(dp)::pi=3.14159265
  real(dp)::x0,y0,z0
  integer:: i,j,id,iu,iv,iw,ip  
  dx=boxlen/2.0d0**(levelmin)
 
id=1 ; iu=2 ; iv=3; iw=4; ip=ndim+2

 !compute variables at 'a'
  rho_a = Mdot/(4.0d0*pi*a**2*vinf)

!#if NDIM==2
! rho_a = Mdot/(2.0d0*pi*a*vinf)
!#endif

  P_a = vinf**2* rho_a/(Mach**2*gamma)


!pulsar
zz = 0.0d0  

  DO i=1,nn 
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
!has to be checked
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*rmask) THEN
        q(i,iu) = frac(xx,yy,zz,rmask,dx)*vinf*(xx/rc)+(1.0d0-frac(xx,yy,zz,rmask,dx)*q(i,iu))
        q(i,iv) = frac(xx,yy,zz,rmask,dx)*vinf*(yy/rc)+(1.0d0-frac(xx,yy,zz,rmask,dx)*q(i,iv))
        q(i,iw) = frac(xx,yy,zz,rmask,dx)*vinf*(zz/rc)+(1.0d0-frac(xx,yy,zz,rmask,dx)*q(i,iw))
        q(i,id) = frac(xx,yy,zz,rmask,dx)*int_rhoP(xx,yy,zz,Mdot,rmask,dx,a,vinf)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id)
        q(i,ip) = frac(xx,yy,zz,rmask,dx)*int_pP(xx,yy,zz,P_a,rho_a,Mdot,a,rmask,vinf,gamma,dx)&
             +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,ip)
     
      ENDIF
#endif

#if NDIM==2   

     IF (rc .lt. 1.2d0*rmask) THEN
        q(i,iu) = frac(xx,yy,zz,rmask,dx)*vinf*(xx/rc)+(1.0d0-frac(xx,yy,zz,rmask,dx)*q(i,iu))
        q(i,iv) = frac(xx,yy,zz,rmask,dx)*vinf*(yy/rc)+(1.0d0-frac(xx,yy,zz,rmask,dx)*q(i,iv))
        q(i,id) = frac(xx,yy,zz,rmask,dx)*int_rhoP(xx,yy,zz,Mdot,rmask,dx,a,vinf)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id)
        q(i,ip) = frac(xx,yy,zz,rmask,dx)*int_pP(xx,yy,zz,P_a,rho_a,Mdot,a,rmask,vinf,gamma,dx)&
             +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,ip)

      ENDIF
#endif
   
   ENDDO

end subroutine init_puls


subroutine reset_mask(ilevel) ! Loop over grids by vector sweeps
! routine which reinitializes the wind at each timestep
  use hydro_parameters
  use amr_commons
  use hydro_commons
  use cooling_module
  use wind_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif


! trier les variables inutiles!!!!!!!!!!!!!!!!!
  integer::ilevel
  integer::i,icell,igrid,ncache,iskip,ngrid,ilun,ncell
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:nvar)     ::u_mask
  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane
  real(dp)::rc1,xc1,yc1,zc1,xc2,yc2,zc2,rc2
  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar
  real(dp)::a1=1.0d0  ! separation of the system in the simulation, it's the leght scale
  real(dp)::v1,v2,rho1,rho2
  real(dp)::delta_x
  real(dp)::pex,pey,pez
  real(dp)::ax02,ay02,r,rot

  ncell=ncoarse+twotondim*ngridmax


!  if (rotation .eqv. .false.)then
!     x02 = x01 +a1             
!     y02 = y01 
!     z02 = z01 
!  else
! the softening length (=rmask) has to be taken into acount in the rotation speed
!     rot=sqrt(GM/(a*(a**2+rmask1**2)))
!     x02=x01+a1*cos(rot*t)
!     y02=y01+a1*sin(rot*t)
! endif
!circular orbit only
!     if (nstep == 1) then
!        x02   = x01+a1*(1.0d0-exc)
!        y02   = y01
!        vx02=0.0d0
!        vy02=sqrt(GM*(2.0D0/(a1-exc)-1.0d0/a1))
!     else   
!        r=sqrt((x02-x01)**2+(y02-y01)**2)
!        ax02=-GM/r**2*(x02-x01)/r
!        ay02=-GM/r**2*(y02-y01)/r
!        vx02=vx02+ax02*dtnew(ilevel)/2.0d0
!        vy02=vy02+ay02*dtnew(ilevel)/2.0d0
!        endif
  !   else  
 !    write(*,*),dtnew(ilevel)
!     call leapfrog(x02,y02,vx02,vy02,dtnew(nlevelmax))
!write(*,*),x02,y02,vx02,vy02,dtnew(ilevel),'main'
!endif
!endif

  call leapfrog(x01,y01,z01,vx01,vy01,vz01,x02,y02,z02,vx02,vy02,vz02,dtnew(nlevelmax),boxlen)


  v2  = vinf2*(1.0d0-rstar2)**beta2 
 rho2 = Mdot2/(4.0d0*3.1415*a**2*v2)

!#if NDIM==2
! rho2 = Mdot2/(2.0d0*3.1415*a*v2)
!#endif


 ! Conversion factor from user units to cgs units

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in coarse cell units
   dx=0.5D0**ilevel
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do igrid=1,ncache,nvector
  ngrid=MIN(nvector,ncache-igrid+1)
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  ! Loop over cells
  do ind=1,twotondim
   ! Gather cell indices
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do
   ! Gather cell centre positions
! get size of cell
     do idim=1,ndim
         do i=1,ngrid
            delta_x = boxlen/2.0d0**ilevel
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
     end do
           ! Rescale position from code units to user units
     do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
     end do
     ! dx_loc,utile?????


    call init2_puls(xx,u_mask,dx_loc,ngrid,x01,y01,z01,rmask1,Mdot1,vinf1,Mach1,dx)

!     if (type2 .eq. 'puls') then 
         call init2_puls(xx,u_mask,dx_loc,ngrid,x02,y02,z02,rmask2,Mdot2,vinf2,Mach2,dx)
 !    endif    

     do i=1,ngrid
          xc2=xx(i,1)-x02 
          yc2=xx(i,2)-y02
          rc2=sqrt(xc2**2+yc2**2)
          xc1=xx(i,1)-x01 
          yc1=xx(i,2)-y01
          rc1=sqrt(xc1**2+yc1**2)

! a regrouper avec le reste du 3D?
#if NDIM ==3 
         zc2=xx(i,3)-z02
         rc2=sqrt(xc2**2+yc2**2+zc2**2) 
         zc1=xx(i,3)-z01
         rc1=sqrt(xc1**2+yc1**2+zc1**2) 
#endif
         if (rc1 .lt. 1.2*rmask1) then
!!!! voir pout rho et v!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call  reset_puls(xc1,yc1,zc1,rmask1,dx,u_mask,i,ind_cell(i))
!ambient medium!!! radiation pressure.
         endif  
        if (rc2 .lt. 1.2d0*rmask2)then 
           if (type2 .eq. 'puls') then
               call reset_puls(xc2,yc2,zc2,rmask2,dx,u_mask,i,ind_cell(i))
          endif     
        endif
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask


subroutine reset_star(xc,yc,zc,rmask,delta_x,u_mask,i,ind_cell,ilevel,vinf,a,rstar,rc)
                                       
!routine which changes re-initializes the values of the hydro-variables in the mask zones. It is applied to the case of a star. 
!The imput contains the distance to the stellar centre, the size of the star, the size of a cell
! the indice of the cell and the values of the hydro variables in the mask.
!The routine modifies the global variable uold.
  use hydro_commons
  use amr_parameters ! sur?
  use hydro_parameters
  use amr_commons
  implicit none
  real(dp)::frac,vflu
  real(dp),dimension(1:nvector,1:nvar):: u_mask
  real(dp):: outm!,pex,pey,pez
  integer :: i,ilevel,j
  real(dp)::xc,yc,zc,delta_x,rc,rmask,vinf,a,rstar,rho,v
  integer ::ind_cell,nn

 
    do j=1,nvar
  !  uold(ind_cell,j)=frac(xc,yc,zc,rmask,delta_x)*u_mask(i,j)+(1.d0-frac(xc,yc,zc,rmask,delta_x))*uold(ind_cell,j)
       uold(ind_cell,j)=u_mask(i,j)
    enddo
  
end subroutine reset_star

subroutine reset_puls(xc,yc,zc,r0,delta_x,u_mask,i,ind_cell)              
!routine which changes re-initializes the values of the hydro-variables in the mask zones. It is applied to the case of a pulsar. 
!The imput contains the distance to the stellar centre, the size of the pulsar, the size of a cell
! the indice of the cell and the values of the hydro variables in the mask.
!The routine modifies the global variable uold.

  use hydro_commons
  use amr_parameters
  use hydro_parameters
  implicit none
  real(dp)::frac
  real(dp),dimension(1:nvector,1:nvar):: u_mask
  integer ::ivar,i
  integer ::ind_cell
  real(dp)::xc,yc,zc,r0,delta_x
  
      do ivar=1,nvar
           uold(ind_cell,ivar)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,ivar)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,ivar)
     enddo   

end subroutine reset_puls

subroutine init2_star(x,u,dx,nn,x0,y0,z0,rmask,Mdot,vinf,Mach,rstar,beta)
    
!This routine computes the hydro variables in the mask region, in the case of a star
!it returns the vector u and requires all the parameters of the star as input
  use wind_parameters
  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  implicit none

  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer:: i,j,id,iu,iv,iw,ip
  real(dp)::rc,rs,xx,yy,zz,r3d
   integer::  ivar
   real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables !!! save
  real(dp):: int_rho,int_P,H,density
  real(dp):: int_vx,int_vy,int_vz
  real(dp):: delta_x
  real(dp):: boxlen_0,x0,y0,z0,Mdot,vinf,rstar,beta,vrotx,vroty,omega,testx,testy
  real(dp)::v_a,rho_a,P_a,Mach,v,vrot,vwind,rmask,v_rmask,rho_0w,csound,frac

  
  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2


  v_rmask=vinf*(1.0d0-rstar/rmask)**beta
  rho_0w=Mdot/(4.0d0*3.1415*(a*rmask)**2*v_rmask)
!#if NDIM==2
! rho_0w = Mdot/(2.0d0*3.1415*(a*rmask)*v_rmask)
!#endif


  zz = 0.0d0  

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center
     csound=cs_0*(r_0/rmask)**(qu/2.)

#if NDIM == 3
     zz=x(i,3)-z0
     r3d=sqrt(xx**2+yy**2+zz**2)
     if (disk_presence .eqv. .true.)then
        omega=sqrt(GM/(r3d*(rmask**2+r3d**2))) ! non pressure gradient correction yet!!!!)+csound**2*(pe+qu)/(rc**2)) 
     endif
     IF (r3d .lt. 1.2d0*rmask) THEN
        H=H_0*(rc/r_0)**((3.-qu)/2.)
        q(i,id)=rho_0*(r_0/rmask)**pe
        q(i,ip)=csound**2*q(i,id)/gamma
        vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
        q(i,iu)=-omega*yy+vwind
        vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
        q(i,iv)=omega*xx+vwind
        q(i,iw) =  int_vz(xx,yy,zz,vinf,rstar,beta,dx)
        ! Convert primitive to conservative variables
        ! density -> density
        u(i,id)=q(i,id)
        ! velocity -> momentum
        u(i,iu)=q(i,id)*q(i,iu)
        u(i,iv)=q(i,id)*q(i,iv)
        u(i,iw)=q(i,id)*q(i,iw)
        ! kinetic energy
        u(i,ip)=0.0d0
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iu)**2
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iv)**2
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iw)**2
        ! pressure -> total fluid energy
        u(i,ip)=u(i,ip)+q(i,ip)/(gamma-1.0d0)
        ! passive scalars
           do ivar=ndim+3,nvar
              u(i,ivar)=q(i,1)*q(i,ivar)
           end do

        ENDIF
#endif

#if NDIM==2   

     IF (rc .lt. 1.2d0*rmask) THEN
 ! In 2D there is either a disk, or a wind, there is no need having them both. Yet there has to be as mask. 
! The density has to be continuous across the edge of the disk.
        if (disk_presence .eqv..true.) then
           q(i,id)= rho_0*(r_0/rmask)**pe !frac(xx,yy,zz,rmask,dx)* rho_0*(r_0/rmask)**pe +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id) 
           q(i,ip)=csound**2*q(i,id)/gamma
           omega=sqrt(GM/(rc*(rmask**2+rc**2))+csound**2*(pe+qu)/(rc**2)) 
           q(i,iu)=-omega*yy
           q(i,iv)=omega*xx
           if ((q(i,iu)*xx+q(i,iv)*yy).gt. 1.0e-14) then
              write(*,*)q(i,iu)*xx+q(i,iv)*yy,'vitesse radiale star'
           endif
        else   
           q(i,id)=rho_0w!*(r_0/rmask)**pe
           q(i,ip)=csound**2*q(i,id)/gamma
           vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
           q(i,iu)=vwind
           vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
           q(i,iv)=vwind
        endif   
        ! Convert primitive to conservative variables
        ! density -> density
        u(i,id)=q(i,id)
        ! velocity -> momentum
        u(i,iu)=q(i,id)*q(i,iu)
        u(i,iv)=q(i,id)*q(i,iv)
        ! kinetic energy
        u(i,ndim+2)=0.0d0
        u(i,ndim+2)=u(i,ndim+2)+0.5*q(i,id)*q(i,iu)**2
        u(i,ndim+2)=u(i,ndim+2)+0.5*q(i,id)*q(i,iv)**2
        ! pressure -> total fluid energy
        u(i,ndim+2)=u(i,ndim+2)+q(i,ndim+2)/(gamma-1.0d0)
        !passive scalars
         do ivar=ndim+3,nvar
            u(i,ivar)=q(i,1)*q(i,ivar)
         end do

        ENDIF
#endif
   ENDDO 


end subroutine init2_star

subroutine init2_puls(x,u,dx,nn,x0,y0,z0,r0,Mdot,vinf,Mach,delta_x)
 
  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  use wind_parameters 
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  
  integer:: i,j,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0
  real(dp):: rstar, beta
  integer::  ivar
  real(dp):: vinf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: int_rhoP,int_Pp
  real(dp):: delta_x
  real(dp):: Mdot,Mach
  real(dp)::v_a,rho_a,p_a
  IF (nn .eq. 0) RETURN
!faudra voir si tous les use sont utiles

  rho_a = Mdot/(4.0d0*3.141592*a**2*vinf)

!#if NDIM==2
! rho_a = Mdot/(2.0d0*3.141592*a*vinf)
!#endif

  P_a = vinf**2* rho_a/(Mach**2*gamma)

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  zz=0.0d0 

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,iu) = vinf*(xx/rc)
        q(i,iv) = vinf*(yy/rc)
        q(i,iw) = vinf*(zz/rc)
        q(i,id) = int_rhoP(xx,yy,zz,Mdot,r0,delta_x,a,vinf)
        q(i,ip) = int_pP(xx,yy,zz,P_a,rho_a,Mdot,a,r0,vinf,gamma,delta_x)
        ! Convert primitive to conservative variables

        ! density -> density
        u(i,id)=q(i,id)
        ! velocity -> momentum
        u(i,iu)=q(i,id)*q(i,iu)
        u(i,iv)=q(i,id)*q(i,iv)
        u(i,iw)=q(i,id)*q(i,iw)
        ! kinetic energy
        u(i,ip)=0.0d0
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iu)**2
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iv)**2
        u(i,ip)=u(i,ip)+0.5*q(i,id)*q(i,iw)**2
        ! pressure -> total fluid energy
        u(i,ip)=u(i,ip)+q(i,ip)/(gamma-1.0d0)
        ! passive scalars
           do ivar=ndim+3,nvar
              u(i,ivar)=q(i,1)*q(i,ivar)
           end do

        ENDIF
#endif


#if NDIM==2   
        IF (rc .lt. 1.2d0*r0) THEN
         q(i,iu) = vinf*(xx/rc)
         q(i,iv) = vinf*(yy/rc)
         q(i,id) = int_rhoP(xx,yy,zz,Mdot,r0,delta_x,a,vinf)
         q(i,ip) = int_pP(xx,yy,zz,P_a,rho_a,Mdot,a,r0,vinf,gamma,delta_x)
            ! Convert primitive to conservative variables
           ! density -> density
           u(i,id)=q(i,id)
           ! velocity -> momentum
           u(i,iu)=q(i,id)*q(i,iu)
           u(i,iv)=q(i,id)*q(i,iv)
           ! kinetic energy
           u(i,ndim+2)=0.0d0
           u(i,ndim+2)=u(i,ndim+2)+0.5*q(i,id)*q(i,iu)**2
           u(i,ndim+2)=u(i,ndim+2)+0.5*q(i,id)*q(i,iv)**2
           ! pressure -> total fluid energy
           u(i,ndim+2)=u(i,ndim+2)+q(i,ndim+2)/(gamma-1.0d0)
           !passive scalars
           do ivar=ndim+3,nvar
              u(i,ivar)=q(i,1)*q(i,ivar)
           end do
  
        ENDIF
#endif
   ENDDO 


end subroutine init2_puls


function int_rho(xx,yy,zz,Mdot,r0,dx,a,vinf,rstar,beta)
  use amr_parameters
implicit none


real(kind=8)::int_rho
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::rho
real(kind=8)::int
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::vol
real(kind=8)::Mdot,a,vinf,rstar,beta
real(kind=8)::r0 
real(kind=8)::dx


vol = dx**ndim

int =0.0d0

#if NDIM==2
   do i=1,n
       do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         int = int+(dx/n)**ndim*rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),zz,Mdot,r0,a,vinf,rstar,beta)
 
   
     end do
   end do  
#endif

#if NDIM==3
   do i=1,n
       do j=1,n
          do k=1,n
         int = int+(dx/n)**ndim*rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),(zz+(k-5.5d0)*dx/n),Mdot,r0,a,vinf,rstar,beta)
          end do
       end do
   end do
#endif

int_rho=int/vol

return
end function int_rho

function rho(x,y,z,Mdot,r0,a,vinf,rstar,beta) 
  use amr_parameters
implicit none

real(kind=8) :: rho,v
real(kind=8)::x,y,z
real(kind=8)::r0 
real(kind=8)::r 
real(kind=8)::Mdot,a,vinf,rstar,beta

r=sqrt(x**2+y**2+z**2)
 
#if NDIM==2
      if (r==0.0d0)then
         r=rstar/100d0
      endif
        rho=a**2*Mdot/(4.0d0*3.1415d0*r**2*v(r,vinf,rstar,beta))

 !       rho=a*Mdot/(2.0d0*3.1415d0*r*v(r,vinf,rstar,beta))
!!      rho=Mdot/(2.0d0*3.1415d0*r*a*v(r,vinf,rstar,beta))
#endif
 
#if NDIM==3
     if (r==0.0d0)then
         r=rstar/100d0
     endif
        rho=a**2*Mdot/(4.0d0*3.1415d0*r**2*v(r,vinf,rstar,beta))
!!     rho=Mdot/(4.0d0*3.1415d0*(r*a)**2*v(r,vinf,rstar,beta))
#endif

return
end function rho


function int_P(xx,yy,zz,P_a,rho_a,Mdot,a,r0,vinf,rstar,beta,gamma,dx)
  use amr_parameters
implicit none


real(kind=8)::int_P
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::rho,rho_r
real(kind=8)::int
real(kind=8)::xx,yy,zz   
real(kind=8)::vol
real(kind=8)::P_a,rho_a! values at 'a'
real(kind=8)::gamma,Mdot,vinf,rstar,beta,a
real(kind=8)::dx,r0

vol = dx**ndim

int =0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         rho_r = rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),zz,Mdot,r0,a,vinf,rstar,beta)
         int =int +(dx/n)**ndim*P_a*(rho_a/rho_r)**(-gamma) 
      end do
   end do
#endif

#if NDIM==3
    do k=1,n
       do i=1,n
          do j=1,n
             rho_r=rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),(zz+(k-5.5d0)*dx/n),Mdot,r0,a,vinf,rstar,beta)
             int  =int +(dx/n)**ndim*P_a*(rho_a/rho_r)**(-gamma)
          end do
       end do
    end do
#endif


int_P=int/vol


return
end function int_P

function int_vx(xx,yy,zz,vinf,rstar,beta,dx)
 use amr_parameters
implicit none


real(kind=8)::int_vx
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::v
real(kind=8)::int
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::vol
real(kind=8)::vinf
real(kind=8)::rstar
real(kind=8)::dx
real(kind=8)::beta
real(kind=8)::r

vol = dx**ndim

int = 0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2)
         if (r==0.0d0)then
             r=rstar/100d0
         endif
         int =int +(dx/n)**ndim*(xx/r)*v(r,vinf,rstar,beta)
       end do
   end do
#endif


#if NDIM==3 
   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2+(zz+(k-5.5d0)*dx/n)**2)
            if (r==0.0d0)then
               r=rstar/100d0
            endif
            int =int +(dx/n)**ndim*(xx/r)*v(r,vinf,rstar,beta)
         end do
      end do
   enddo
#endif

int_vx=int/vol
return
end function int_vx



function int_vy(xx,yy,zz,vinf,rstar,beta,dx)
  use amr_parameters
implicit none


real(kind=8)::int_vy
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::v
real(kind=8)::int
real(kind=8)::xx,yy,zz  
real(kind=8)::vol
real(kind=8)::dx ! size of the cell
real(kind=8)::vinf
real(kind=8)::rstar
real(kind=8)::beta
real(kind=8)::r

vol= dx**ndim

int=0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2)
         if (r==0.0d0)then
            r=rstar/100d0
         endif
         int =int +(dx/n)**ndim*(yy/r)*v(r,vinf,rstar,beta)
      end do
    end do
#endif

#if NDIM==3
   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2+(zz+(k-5.5d0)*dx/n)**2)
            if (r==0.0d0)then
               r=rstar/100d0
            endif
            int =int +(dx/n)**ndim*(yy/r)*v(r,vinf,rstar,beta)
         end do
      end do
   end do 
#endif

int_vy=int/vol

return
end function int_vy


function int_vz(xx,yy,zz,vinf,rstar,beta,dx)
  use amr_parameters
implicit none
real(kind=8)::int_vz
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::v
real(kind=8)::int
real(kind=8)::xx,yy,zz  
real(kind=8)::vol
real(kind=8)::dx ! size of the cell
real(kind=8)::vinf
real(kind=8)::rstar
real(kind=8)::beta
real(kind=8)::r

vol = dx**ndim

int = 0.0d0

   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2+(zz+(k-5.5d0)*dx/n)**2)
            if (r==0.0d0)then
               r=rstar/100d0
            endif
            int =int +(dx/n)**ndim*(zz/r)*v(r,vinf,rstar,beta)
         end do
      end do
   end do 
int_vz=int/vol
return
end function int_vz

function v(r,vinf,rstar,beta) 
implicit none
real(kind=8):: v
real(kind=8):: rstar 
real(kind=8):: vinf   
real(kind=8):: r 
real(kind=8):: beta  

v=max((vinf*(1.0d0-rstar/r)**beta),vinf/10)

return
end function v


function int_rhoP(xx,yy,zz,Mdot,r0,dx,a,vinf)
  use amr_parameters
implicit none
real(kind=8)::int_rhoP
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::int
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::vol
real(kind=8)::Mdot,a,vinf
real(kind=8)::r0 
real(kind=8)::dx
real(kind=8)::r

vol= dx**ndim

int=0.0d0

#if NDIM==2
   do i=1,n
       do j=1,n
          r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2)
          if (r==0.0d0)then
             write(*,*),'r=0',r0
             r=r0/100d0
          endif
!          int= int+(dx/n)**ndim*(1.0d0/r)
            int =int+(dx/n)**ndim*1.0d0/(r**2)

          ! int= int+(dx/n)**ndim*min(1.0d0/(r**2),1.0d0/(0.01*r0))
       end do
   end do  
int=int*Mdot*a/(2*3.1415*vinf)
#endif

#if NDIM==3
   do i=1,n
       do j=1,n
          do k=1,n
            r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2+(zz+(k-5.5d0)*dx/n)**2)
            if (r==0.0d0)then
               r=r0/100d0
            endif
            int =int+(dx/n)**ndim*1.0d0/(r**2)
!  int= int+(dx/n)**ndim*min(1.0d0/(r**2),1.0d0/(0.01*r0))
            
      !int= int+(dx/n)**ndim*min(1.0d0/(r**2),1.0d0(0.1*r0))
         end do
       end do
   end do
int=int*Mdot*a**2/(4*3.1415*vinf)
#endif

int_rhoP=int/vol
return
end function int_rhoP


function int_Pp(xx,yy,zz,P_a,rho_a,Mdot,a,r0,vinf,gamma,dx)
  use amr_parameters
implicit none
real(kind=8)::int_Pp
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::rho
real(kind=8)::int
real(kind=8)::xx,yy,zz   
real(kind=8)::vol
real(kind=8)::P_a,rho_a
real(kind=8)::gamma,Mdot,vinf,a
real(kind=8)::dx,r0!,P_r0
real(kind=8)::r!,rho_r

vol= dx**ndim

int=0.0d0
!P_r0=P_a*(r0/a)**gamma

#if NDIM==2
   do i=1,n
      do j=1,n
         r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2)
         if (r==0.0d0)then
            r=r0/100d0
            write(*,*),'r=0 P',r0
         endif
         !rho_r=Mdot/(2*3.1415*vinf*r)
 !        int =int +(dx/n)**ndim*P_a*(min((a/r)**(gamma),10*(a/r0)**gamma))
           int =int +(dx/n)**ndim*P_a*(min((a/r)**(2*gamma),10*(a/r0)**(2*gamma)))
  
    end do
   end do
#endif

#if NDIM==3
    do k=1,n
       do i=1,n
          do j=1,n
            r=sqrt((xx+(i-5.5d0)*dx/n)**2+(yy+(j-5.5d0)*dx/n)**2+ (zz+(k-5.5d0)*dx/n)**2)
            if (r==0.0d0)then
               r=r0/100d0
            endif
          !  rho_r=Mdot/(4*3.1415*vinf*r**2)
            int =int +(dx/n)**ndim*P_a*(min((a/r)**(2*gamma),10*(a/r0)**(2*gamma)))
          end do
       end do
    end do
#endif

int_pP=int/vol


return
end function int_pP

function f(vinf,a,rstar,rho_a,v_a,rc)
implicit none
real(kind=8)::f
real(kind=8)::vinf
real(kind=8)::a
real(kind=8)::rstar
real(kind=8)::rho_a
real(kind=8)::v_a
real(kind=8)::rc

f=(a*vinf*v_a*rho_a*rstar)/(rc**3)

return
end function f


function frac(x,y,z,r0,dx)
 use amr_parameters

!function which computes the fraction of a cell being occupied by the mask
implicit none
real(kind=8):: frac
real(kind=8):: dx
real(kind=8):: x,y,z,r
real(kind=8):: xc,yc,zc! postion on the center of the slice
real(kind=8):: r0
integer :: i,j,k
integer :: n=10 

frac=0.d0

#if NDIM==2
   do i=1,n
      do j=1,n
         xc=x+(i-5.5d0)*dx/n
         yc=y+(j-5.5d0)*dx/n
         r=sqrt(xc**2+yc**2)
         if (r .lt. r0)then
            frac=frac+1
         endif
      end do
   end do
#endif

#if NDIM==3
   do k=1,n
      do i=1,n
         do j=1,n
            xc=x+(i-5.5d0)*dx/n
            yc=y+(j-5.5d0)*dx/n
            zc=z+(k-5.5d0)*dx/n
            r=sqrt(xc**2+yc**2+zc**2)
            if (r .lt. r0)then
               frac=frac+1
            endif
         end do
      end do
   end do   
#endif

frac=frac/(n**ndim)

return
end function frac


subroutine read_wind_params(nml_ok)
!subroutine which reads the additionnal namelist
use wind_parameters
use amr_parameters ! pour boxlen
use hydro_parameters
implicit none

 logical::nml_ok
 real::omega_0

 namelist/disk/disk_presence,rho_0,pe,qu,GM,x01,y01,z01,r_0,rmask1,HR
!!  namelist/first_body/rstar1,type1,Mdot1,Mach1,vinf1,rmask1,beta1,f1
! namelist/second_body/rstar2,type2,Mdot2,Mach2,vinf2,rmask2,beta2,f2
 namelist/first_body/type1,Mdot1,Mach1,vinf1,rmask1,rstar1,beta1,f1
 namelist/second_body/type2,Mdot2,Mach2,vinf2,rmask2,rstar2,beta2,f2
 namelist/ambient_medium/a,rotation,exc

 INQUIRE(file='param.nml',exist=nml_ok)
  if(.not. nml_ok)then
    ! if(myid==1)then
    !    write(*,*)'File '//TRIM(infile)//' does not exist'
    ! endif
     write(*,*),'The param.nml is missing'
     call clean_stop
  end if

  open(unit=10,file='param.nml')
  rewind(10)
  read(10,NML=disk)
!  write(*,*),disk_presence,rho_0,pe,qu,GM,x01,y01,z01,r_0,rmask1
  rewind(10)
!  write(*,*),'test'
  read(10,NML=first_body)
  rewind(10)
  read(10,NML=second_body)
  rewind(10)
  read(10,NML=ambient_medium)
  close(10)

!conversion according to the value of a
!  vinf1=vinf1/a
  vinf2=vinf2/a
  
  r_0  = r_0/a
  !r02  = r02/a

 ! rstar1 = rstar1/a
  rstar2 = rstar2/a
  rmask2  =rmask2/a
  rmask1 = rmask1/a

  GM=GM/a**3
  rho_0=rho_0*a**3

!adaptation to the size of the box
  x01 = x01*boxlen
  y01 = y01*boxlen
  z01 = z01*boxlen


!determination of all the other parameters
!  r_0=2*rmask1
  H_0=HR*r_0
  omega_0=sqrt(GM/r_0**3)
  cs_0=H_0*omega_0
  p_0=cs_0**2*rho_0/gamma

  rho_amb=1.0e-8*rho_0
  P_amb=1.0e-8*P_0 




end subroutine  read_wind_params

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!Subroutine which computes the positions of the bodies at each timestep
subroutine leapfrog_bis(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,dt,boxlen)
implicit none

real(kind=8)::x1,y1,z1,vx1,vy1,vz1,xf1,yf1,zf1,vxf1,vyf1,vzf1
real(kind=8)::x2,y2,z2,vx2,vy2,vz2,xf2,yf2,zf2,vxf2,vyf2,vzf2
real(kind=8)::ax1,ay1,az1,ax2,ay2,az2,ax1f,ay1f,az1f,ax2f,ay2f,az2f
real(kind=8)::dt,boxlen
real(kind=8)::r1,r2,r,rf
real(kind=8)::GMstar,GMpuls
real(kind=8)::xc,yc,zc,a

GMstar=73.d0
GMpuls=7.d0
xc=0.5*boxlen
yc=0.5*boxlen
zc=0.5*boxlen

!write(*,*),x,y,vx,vy,det,'leapfrog'
!moving the star

!write(*,*),x1,y1,z1,x2,y2,z2
!   r1=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   a=1.
   !acceleration at t
   !ax1=-GMpuls/(r1**2)*(x1-x2)/r1
   !ay1=-GMpuls/(r1**2)*(y1-y2)/r1
   !az1=-GMpuls/(r1**2)*(z1-z2)/r1

   !ax2=-GMstar/(r1**2)*(x2-x2)/r1
   !ay2=-GMstar/(r1**2)*(y2-y2)/r1
   !az2=-GMstar/(r1**2)*(z2-z2)/r1
   r=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  
   ax1=-GMpuls*(x1-x2)!/r**3.0d0
   ay1=-GMpuls*(y1-y2)!/r**3.0d0
   az1=-GMpuls*(z1-z2)!/r**3.0d0

   ax2=-GMstar*(x2-x1)!/r**3.0d0
   ay2=-GMstar*(y2-y1)!/r**3.0d0
   az2=-GMstar*(z2-z1)!/r**3.0d0


   xf1=x1+vx1*dt+ax1*dt**2.0d0/2.d0
   yf1=y1+vy1*dt+ay1*dt**2.0d0/2.d0
   zf1=z1+vz1*dt+az1*dt**2.0d0/2.d0

   xf2=x2+vx2*dt+ax2*dt**2.0d0/2.d0
   yf2=y2+vy2*dt+ay2*dt**2.0d0/2.d0
   zf2=z2+vz2*dt+az2*dt**2.0d0/2.d0

   rf=sqrt((xf1-xf2)**2+(yf1-yf2)**2+(zf1-zf2)**2)
   write(*,*),r
! acceleration at r+dt
 !  ax1f=-GMpuls/(rf**2)*(xf1-xf2)/rf
 !  ay1f=-GMpuls/(rf**2)*(yf1-yf2)/rf
 !  az1f=-GMpuls/(rf**2)*(zf1-zf2)/rf

 !  ax2f=-GMstar/(rf**2)*(xf2-xf1)/rf
 !  ay2f=-GMstar/(rf**2)*(yf2-yf1)/rf
 !  az2f=-GMstar/(rf**2)*(zf2-zf1)/rf

   ax1f=-GMpuls*(xf1-xf2)!/rf**3.0d0
   ay1f=-GMpuls*(yf1-yf2)!/rf**3.0d0
   az1f=-GMpuls*(zf1-zf2)!/rf**3.0d0

   ax2f=-GMstar*(xf2-xf1)!/rf**3.0d0
   ay2f=-GMstar*(yf2-yf1)!/rf**3.0d0
   az2f=-GMstar*(zf2-zf1)!/rf**3.0d0

   vxf1=vx1+(ax1+ax1f)/2.d0*dt
   vyf1=vy1+(ay1+ay1f)/2.d0*dt
   vzf1=vz1+(az1+az1f)/2.d0*dt
 
   vxf2=vx2+(ax2+ax2f)/2.d0*dt
   vyf2=vy2+(ay2+ay2f)/2.d0*dt
   vzf2=vz2+(az2+az2f)/2.d0*dt
 
   x2=xf2
   y2=yf2
   z2=zf2

   vx2=vxf2
   vy2=vyf2
   vz2=vzf2


   x1=xf1
   y1=yf1
   z1=zf1

   vx1=vxf1
   vy1=vyf1
   vz1=vzf1



end subroutine leapfrog_bis


subroutine leapfrog(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,G,M1,M2,dt,boxlen)
!use wind_parameters
implicit none!



real(kind=8)::x1,y1,z1,vx1,vy1,vz1,xf1,yf1,zf1,vxf1,vyf1,vzf1
real(kind=8)::x2,y2,z2,vx2,vy2,vz2,xf2,yf2,zf2,vxf2,vyf2,vzf2
real(kind=8)::ax1,ay1,az1,ax2,ay2,az2,ax1f,ay1f,az1f,ax2f,ay2f,az2f
real(kind=8)::dt,boxlen
real(kind=8)::r1,r2,r,rf
real(kind=8)::M2,M1,G
real(kind=8)::xc,yc,zc,a




   r=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)


   ax1=-G*m2*(x1-x2)/r**3.0d0                                                                                                                                                                               
   ay1=-G*m2*(y1-y2)/r**3.0d0   
   az1=-G*m2*(z1-z2)/r**3.0d0                                                                                                                                                                               

   ax2=-G*m1*(x2-x1)/r**3.0d0                                                                                                                                                                               
   ay2=-G*m1*(y2-y1)/r**3.0d0                                                                                                                                                                               
   az2=-G*m1*(z2-z1)/r**3.0d0                                                                                                                                                                               


   xf1=x1+vx1*dt+ax1*dt**2.0d0/2.d0
   yf1=y1+vy1*dt+ay1*dt**2.0d0/2.d0
   zf1=z1+vz1*dt+az1*dt**2.0d0/2.d0

   xf2=x2+vx2*dt+ax2*dt**2.0d0/2.d0
   yf2=y2+vy2*dt+ay2*dt**2.0d0/2.d0
   zf2=z2+vz2*dt+az2*dt**2.0d0/2.d0

   rf=sqrt((xf1-xf2)**2+(yf1-yf2)**2+(zf1-zf2)**2)

   ax1f=-G*m2*(xf1-xf2)/rf**3.0d0  
   ay1f=-G*m2*(yf1-yf2)/rf**3.0d0  
   az1f=-G*m2*(zf1-zf2)/rf**3.0d0      
   ax2f=-G*m1*(xf2-xf1)/rf**3.0d0      
   ay2f=-G*m1*(yf2-yf1)/rf**3.0d0      
   az2f=-G*m1*(zf2-zf1)/rf**3.0d0          

   vxf1=vx1+(ax1+ax1f)/2.d0*dt
   vyf1=vy1+(ay1+ay1f)/2.d0*dt
   vzf1=vz1+(az1+az1f)/2.d0*dt

   vxf2=vx2+(ax2+ax2f)/2.d0*dt
   vyf2=vy2+(ay2+ay2f)/2.d0*dt
   vzf2=vz2+(az2+az2f)/2.d0*dt


   x1=xf1
   y1=yf1
   z1=zf1

   vx1=vxf1
   vy1=vyf1
   vz1=vzf1


   x2=xf2
   y2=yf2
   z2=zf2

   vx2=vxf2
   vy2=vyf2
   vz2=vzf2






!GMstar=73.d0
!GMpuls=7.d0
!xc=0.5*boxlen
!yc=0.5*boxlen
!zc=0.5*boxlen

!   r=sqrt((x1-x2)**2+(y1-y2)**2)
!   write(*,*),r   


!   r1=sqrt((x1-xc)**2+(y1-yc)**2)
!   r2=sqrt((x2-xc)**2+(y2-yc)**2)
!   ax1=-GMpuls/(r**2)*(x1-x2)/r
!   ay1=-GMpuls/(r**2)*(y1-y2)/r!

!   ax2=-GMstar/(r**2)*(x2-x1)/r
!   ay2=-GMstar/(r**2)*(y2-y1)/r

!   vxf1=vx1+ax1*dt
!   vyf1=vy1+ay1*dt

!   vxf2=vx2+ax2*dt
!   vyf2=vy2+ay2*dt

!   xf1=(x1)++vx1*dt
!   yf1=(y1)+vy1*dt
   
!   xf2=(x2)+vx2*dt
!   yf2=(y2)+vy2*dt

!   x1=xf1
!   y1=yf1
  
!   x2=xf2
!   y2=yf2
   
!   vx1=vxf1
!   vy1=vyf1

!   vx2=vxf2
!   vy2=vyf2

end subroutine leapfrog
