!================================================================
!================================================================
!================================================================
!================================================================

!all 'a' have been removed except in the beginning when everything is scaled

subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  use wind_params
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
  real(dp):: rho_amb,P_amb,v_a,P_a,rho_a
  integer :: ivar,ncell
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp)::a1=1.0d0 ! separation of the system in simulation units
  real(dp)::r,ax,ay
  real(dp)::Mstar1=18
  real(dp)::Mstar2=12
  integer::index_pas_scal1=1
  integer::index_pas_scal2=2

  ncell=ncoarse+twotondim*ngridmax

  ! If it is the first step, allocate arrays for centered scheme for the force
  if (allocated(px).neqv. .true.)then 
      allocate (px(ncell))
  endif   
 
  if (allocated(py).neqv. .true.)then 
      allocate (py(ncell))
  endif   

  if (allocated(pz).neqv. .true.)then 
      allocate (pz(ncell))
  endif   

 if (allocated(px).neqv. .true.)then 
      write(*,*)'Px could not be allocated'
 endif   
 
  if (allocated(py).neqv. .true.)then 
      write(*,*)'Py could not be allocated'
 endif   

  if (allocated(pz).neqv. .true.)then 
      write(*,*)'Pz could not be allocated'
  endif


  px = 0.0d0
  py = 0.0d0
  pz = 0.0d0

  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2


!à ameliorer!!!!!!!!!
  v_a = vinf1*(1.0d0-rstar1)**beta1 
  rho_a = Mdot1*a1**2/(4.0d0*3.1415*v_a)

#if NDIM==2
 rho_a = Mdot1*a1/(2.0d0*3.1415*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach1**2*gamma)


  P_amb=P_a/R_P
  rho_amb=rho_a/R_rho

  q(:,id) = rho_amb
  q(:,iu) = 0.0d0
  q(:,iv) = 0.0d0
#if NDIM == 3
  q(:,iw) = 0.0d0
#endif
  q(:,ip) = P_amb
!there are no passive scalars in the initial medium

  if (nvar .gt. ip) then
     q(:,ip+1) =0.0d0 
     q(:,ip+2) =0.0d0 
  endif


!initial position of the bodies
!by default the position of the first body is given in the parameter file
  if (rotation .eqv. .true.) then 
     x01 = 0.5*boxlen-(Mstar2)*a1*(1-exc)/(Mstar1+Mstar2)
     x02   = x01+a1*(1.0d0-exc)
  else
     x02  =  x01+a1
  endif
  y02   = y01
  z02   = z01

!if exc neq 0!!!! to do!!!!
!initial velocity field (necessary for the position update)

!  vx02=0.0d0
 
!  vy02=sqrt(Gstar*M*(2.0D0/(x02-x01)-1.0d0/a1))
!!!!! a verifier!!!
!write(*,*),x02,y02,vx02,vy02,'ini'

  
  if (type1 == 'star') then 
   call init_star(q,x,nn,vinf1,rstar1,Mdot1,Mach1,x01,y01,z01,r01,beta1,index_pas_scal1)
 else if (type1 =='puls') then
   call init_puls(q,x,nn,vinf1,r01,Mdot1,Mach1,x01,y01,z01,index_pas_scal1) 
 endif


 if (type2 == 'star') then 
   call init_star(q,x,nn,vinf2,rstar2,Mdot2,Mach2,x02,y02,z02,r02,beta2,index_pas_scal2)
 else if (type2 =='puls') then
   call init_puls(q,x,nn,vinf2,r02,Mdot2,Mach2,x02,y02,z02,index_pas_scal2)
 endif

    
     
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



subroutine init_star(q,x,nn,vinf,rstar,Mdot,Mach,x0,y0,z0,r0,beta,index_pas_scal)

  use amr_parameters
  use wind_params
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: vinf,rstar,r0,Mdot,Mach,beta
  real(dp)::rho_a,P_A,v_a  ! variables for r =a
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables!!!!
  real(dp):: int_rho,int_p,int_vx,int_vy,int_vz ! function which computes the mean value of a variable in the mask
  real(dp):: frac
  real(dp)::pi=3.141459265
  real(dp)::x0,y0,z0
  integer:: i,j,id,iu,iv,iw,ip,index_pas_scal
  real(dp)::a1=1.0d0 ! value of the separation in a AU
 
 id=1 ; iu=2 ; iv=3; iw=4; ip=ndim+2

  dx=boxlen/2.0d0**(levelmin)
 !compute variables at 'a'
  v_a = vinf*(1.0d0-rstar)**beta 
  rho_a = Mdot*a1**2/(4.0d0*pi*v_a)

#if NDIM==2
 rho_a = Mdot*a1/(2.0d0*pi*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach**2*gamma)

  zz = 0.0d0  

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center


!regrouper les (1 -frac...) dans un tableau de taille nvar

#if NDIM == 3 
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,iu) =  frac(xx,yy,zz,r0,dx)* int_vx(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,iu)
        q(i,iv) =  frac(xx,yy,zz,r0,dx)* int_vy(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,iv)
        q(i,iw) =  frac(xx,yy,zz,r0,dx)* int_vz(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,iw)
        q(i,id) =  frac(xx,yy,zz,r0,dx)*int_rho(xx,yy,zz,Mdot,r0,dx,vinf,rstar,beta)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,id)
    q(i,ip) =  frac(xx,yy,zz,r0,dx)*int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,ip)
        if (nvar  .gt. ip) then
           q(i,ip+index_pas_scal)=frac(xx,yy,zz,r0,dx)
        endif

  ENDIF
#endif


#if NDIM==2   
     IF (rc .lt. 1.2d0*r0) THEN
       q(i,iu) =  frac(xx,yy,zz,r0,dx)* int_vx(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,iu)
       q(i,iv) =  frac(xx,yy,zz,r0,dx)* int_vy(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,iv)
       q(i,id) =  frac(xx,yy,zz,r0,dx)*int_rho(xx,yy,zz,Mdot,r0,dx,vinf,rstar,beta)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,id)
   q(i,ip) =  frac(xx,yy,zz,r0,dx)*int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,ip)

       if (nvar  .gt. ip) then
          q(i,ip+index_pas_scal)=frac(xx,yy,zz,r0,dx)
       endif

     ENDIF
#endif
  

  ENDDO     

end subroutine init_star


subroutine init_puls(q,x,nn,vinf,r0,Mdot,Mach,x0,y0,z0,index_pas_scal) 
  use wind_params
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx 
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: vinf,r0,Mdot,Mach,beta
  real(dp)::rho_a,P_A ! variables for r =a
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables
  real(dp):: int_rhoP,int_pP
  real(dp):: frac
  real(dp)::pi=3.14159265
  real(dp)::x0,y0,z0
  integer:: i,j,id,iu,iv,iw,ip  
  real(dp)::a1=1.0d0
  integer::index_pas_scal

  dx=boxlen/2.0d0**(levelmin)
 
id=1 ; iu=2 ; iv=3; iw=4; ip=ndim+2

 !compute variables at 'a'
  rho_a = Mdot*a1**2/(4.0d0*pi*vinf)

#if NDIM==2
 rho_a = Mdot*a1/(2.0d0*pi*vinf)
#endif

  P_a = vinf**2* rho_a/(Mach**2*gamma)


!pulsar
zz = 0.0d0  

  DO i=1,nn 
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,iu) = frac(xx,yy,zz,r0,dx)*vinf*(xx/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,iu)
        q(i,iv) = frac(xx,yy,zz,r0,dx)*vinf*(yy/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,iv)
        q(i,iw) = frac(xx,yy,zz,r0,dx)*vinf*(zz/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,iw)
        q(i,id) = frac(xx,yy,zz,r0,dx)*int_rhoP(xx,yy,zz,Mdot,r0,dx,vinf)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,id)
        q(i,ip) = frac(xx,yy,zz,r0,dx)*int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,ip)
        if (nvar  .gt. ip) then      
           q(i,ip+index_pas_scal) = frac(xx,yy,zz,r0,dx)
        endif


      ENDIF
#endif

#if NDIM==2   

     IF (rc .lt. 1.2d0*r0) THEN
        q(i,iu) = frac(xx,yy,zz,r0,dx)*vinf*(xx/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,iu)
        q(i,iv) = frac(xx,yy,zz,r0,dx)*vinf*(yy/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,iv)
        q(i,id) = frac(xx,yy,zz,r0,dx)*int_rhoP(xx,yy,zz,Mdot,r0,dx,vinf)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,id)
        q(i,ip) = frac(xx,yy,zz,r0,dx)*int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,ip)
        if (nvar  .gt. ip) then
           q(i,ip+index_pas_scal) =frac(xx,yy,zz,r0,dx)
        endif
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
  use wind_params
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
  real(dp)::ax02,ay02,r
  real(dp)::dt ! local timestep
  real(dp)::xbar,ybar !masse center of both stars
  real(dp)::Mstar1=18
  real(dp)::Mstar2=12
  real(dp)::omega
  integer::index_pas_scal1=1
  integer::index_pas_scal2=2

  ncell=ncoarse+twotondim*ngridmax

  !if it is a restart, reallocate the arrays for centered scheme
  if (allocated(px) .neqv. .true.)then
     allocate(px(ncell))
     px=uold(:,2)
  endif
  if (allocated(py) .neqv. .true.)then
     allocate(py(ncell))
     py=uold(:,3)
  endif
  if (allocated(pz) .neqv. .true.)then
     allocate(pz(ncell))
  endif

  if (allocated(px) .neqv. .true.)then
     write(*,*)'unable to allocate Px'
  endif
  if (allocated(py) .neqv. .true.)then
    write(*,*),'unable to allocate Py'
  endif
  if (allocated(pz) .neqv. .true.)then
    write(*,*),'unable to allocate Pz'
  endif

#if NDIM==3
     pz=uold(:,4)
#endif



 !positions of the bodies

  if (rotation .eqv. .false.)then
     x02 = x01 +a1
     y02 = y01 
     z02 = z01 
  else
  !   if (nstep == 1) then
  !   x02 = x01 +a1*(1-exc)
  !   y02 = y01 
  !   z02 = z01 


!        x02   = x01+a1*(1.0d0-exc)
!        y02   = y01
!        vx02=0.0d0
!        vy02=sqrt(Gstar*M*(2.0D0/(x02-x01)-1.0d0/a1))
        
!        r=sqrt((x02-x01)**2+(y02-y01)**2)
!        ax02=-Gstar*M/r**2*(x02-x01)/r
!        ay02=-Gstar*M/r**2*(y02-y01)/r

!        vx02=vx02+ax02*dtnew(ilevel)/2.0d0
!        vy02=vy02+ay02*dtnew(ilevel)/2.0d0
   !  endif   
 !    write(*,*),dtnew(ilevel)
 !    call leapfrog(x02,y02,vx02,vy02,dtnew(ilevel))
!write(*,*),x02,y02,vx02,vy02,dtnew(ilevel),'main'

  !      z02 = z01
   xbar=0.5*boxlen
   ybar=0.5*boxlen

   omega=sqrt(Gstar*(Mstar1+Mstar2)/a1**3)

!for indirect rotation
!    x01=xbar-(Mstar2/(Mstar1+Mstar2))*a1*cos(omega*t)
!   y01=ybar-(Mstar2/(Mstar1+Mstar2))*a1*sin(omega*t)

!   x02=xbar+(Mstar1/(Mstar1+Mstar2))*a1*cos(omega*t)
!   y02=ybar+(Mstar1/(Mstar1+Mstar2))*a1*sin(omega*t)

!for direct rotation
   x01=xbar-(Mstar2/(Mstar1+Mstar2))*a1*cos(omega*t)
   y01=ybar+(Mstar2/(Mstar1+Mstar2))*a1*sin(omega*t)

   x02=xbar+(Mstar1/(Mstar1+Mstar2))*a1*cos(omega*t)
   y02=ybar-(Mstar1/(Mstar1+Mstar2))*a1*sin(omega*t)



!!!voir le cas exc neq 0

 endif   
 
   zc2 = 0.0d0
   zc1 = 0.0d0 


!values for r=a (mettre dans une fonction)
   v1 = vinf1*(1.0d0-rstar1)**beta1 
 rho1 = Mdot1*a1**2/(4.0d0*3.1415*v1)

  v2  = vinf2*(1.0d0-rstar2)**beta2 
 rho2 = Mdot2*a1**2/(4.0d0*3.1415*v2)

#if NDIM==2
 rho1 = Mdot1*a1/(2.0d0*3.1415*v1)
 rho2 = Mdot2*a1/(2.0d0*3.1415*v2)
#endif


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
 
     if (type1=='star')then
        call init2_star(xx,u_mask,dx_loc,ngrid,x01,y01,z01,r01,Mdot1,vinf1,Mach1,rstar1,beta1,dx,index_pas_scal1)
     else if (type1=='puls')then
        call init2_puls(xx,u_mask,dx_loc,ngrid,x01,y01,z01,r01,Mdot1,vinf1,Mach1,dx,index_pas_scal1)
     endif
   

     if (type2=='star')then
        call init2_star(xx,u_mask,dx_loc,ngrid,x02,y02,z02,r02,Mdot2,vinf2,Mach2,rstar2,beta2,dx,index_pas_scal2)
     else if (type2=='puls')then
        call init2_puls(xx,u_mask,dx_loc,ngrid,x02,y02,z02,r02,Mdot2,vinf2,Mach2,dx,index_pas_scal2)
     endif

     do i=1,ngrid
         xc1=xx(i,1)-x01 
         yc1=xx(i,2)-y01
         xc2=xx(i,1)-x02 
         yc2=xx(i,2)-y02
         rc1=sqrt(xc1**2+yc1**2)
         rc2=sqrt(xc2**2+yc2**2)
! a regrouper avec le reste du 3D?
#if NDIM ==3 
         zc1=xx(i,3)-z01
         rc1=sqrt(xc1**2+yc1**2+zc1**2) 
         zc2=xx(i,3)-z02
         rc2=sqrt(xc2**2+yc2**2+zc2**2) 
#endif
         if (rc1 .lt. 1.2d0*r01) then
            ! First body: a star or a pulsar
            if (type1=='star')then 

               call reset_star(xc1,yc1,zc1,r01,dx,u_mask,i,ind_cell(i),ilevel,f1,vinf1,rstar1,rho1,v1,rc1)!#endif

           else if (type1=='puls')then
                call reset_puls(xc1,yc1,zc1,r01,dx,u_mask,i,ind_cell(i))
           endif              
            
        else if ((type2 .ne. 'star').and.(type2.ne.'puls'))then
           !case where there is only 1 body
             f2=0  
           if (f1.eq.1)then
               call  reset_amb(xc1,yc1,zc1,rc1,xc2,yc2,zc2,rc2,ind_cell(i),ilevel,px(ind_cell(i)),py(ind_cell(i)),pz(ind_cell(i)))
               px(ind_cell(i))=uold(ind_cell(i),2)
               py(ind_cell(i))=uold(ind_cell(i),3)  
#if NDIM==3
               pz(ind_cell(i))=uold(ind_cell(i),4)
#endif
           endif    

        else if (rc2 .lt. 1.2d0*r02)then 
            !2nd body: a star or a pulsar
           if (type2=='star')then
              pex=px(ind_cell(i))
              pey=py(ind_cell(i))
              pey=pz(ind_cell(i))
              
              call reset_star(xc2,yc2,zc2,r02,dx,u_mask,i,ind_cell(i),ilevel,f2,vinf2,rstar2,rho2,v2,rc2,pex,pey,pez)
              px(ind_cell(i))=uold(ind_cell(i),2)
              py(ind_cell(i))=uold(ind_cell(i),3)
#if NDIM==3
              pz(ind_cell(i))=uold(ind_cell(i),4)
#endif

          else if (type2=='puls')then    
              call reset_puls(xc2,yc2,zc2,r02,dx,u_mask,i,ind_cell(i))
          endif    
          
        elseif(( rc1 .gt. 1.2d0*r01) .and. (rc2 .gt. 1.2d0*r02))then!
          !ambient medium
          call reset_amb(xc1,yc1,zc1,rc1,xc2,yc2,zc2,rc2,ind_cell(i),ilevel,px(ind_cell(i)),py(ind_cell(i)),pz(ind_cell(i)))

          px(ind_cell(i))=uold(ind_cell(i),2)
          py(ind_cell(i))=uold(ind_cell(i),3)
#if NDIM==3
          pz(ind_cell(i))=uold(ind_cell(i),4)
#endif

       endif
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask




subroutine reset_star(xc,yc,zc,r0,delta_x,u_mask,i,ind_cell,ilevel,fs,vinf,rstar,rho,v,rc)
                                       
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
  real(dp):: outm
  integer :: i,ilevel,fs,ivar
  real(dp)::xc,yc,zc,r0,delta_x,rc,vinf,a,rstar,rho,v,f
  integer ::ind_cell    

      uold(ind_cell,1)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,1)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,1)

#if NDIM==2
      vflu=sqrt(uold(ind_cell,2)**2+uold(ind_cell,3)**2)/uold(ind_cell,1)
 ! the force works so there is an extra energy density related to the speed of the flow
      outm = uold(ind_cell,4)+fs*f(vinf,rstar,rho,v,rc)*vflu*dtnew(ilevel)
      uold(ind_cell,4)= frac(xc,yc,zc,r0,delta_x)*u_mask(i,4)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm

 ! the wind is accelerated by the force
      outm = uold(ind_cell,2)+xc/rc*dtnew(ilevel)*fs*f(vinf,rstar,rho,v,rc)
      uold(ind_cell,2)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,2)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm

      outm = uold(ind_cell,3)+yc/rc*dtnew(ilevel)*fs*f(vinf,rstar,rho,v,rc)          
      uold(ind_cell,3)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,3)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm
#endif

#if NDIM==3
      vflu=sqrt(uold(ind_cell,2)**2+uold(ind_cell,3)**2+uold(ind_cell,4)**2)/uold(ind_cell,1)
              
 ! the force works so there is an extra energy density related to the speed of the flow
      outm = uold(ind_cell,5)+fs*f(vinf,rstar,rho,v,rc)*vflu*dtnew(ilevel)
      uold(ind_cell,5)= frac(xc,yc,zc,r0,delta_x)*u_mask(i,5)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm

 ! the wind is accelerated by the force
      outm = uold(ind_cell,2)+xc/rc*dtnew(ilevel)*fs*f(vinf,rstar,rho,v,rc)
      uold(ind_cell,2)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,2)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm

      outm = uold(ind_cell,3)+yc/rc*dtnew(ilevel)*fs*f(vinf,rstar,rho,v,rc)          
      uold(ind_cell,3)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,3)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm

      outm = uold(ind_cell,4)+yc/rc*dtnew(ilevel)*fs*f(vinf,rstar,rho,v,rc)          
      uold(ind_cell,4)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,4)+(1.d0-frac(xc,yc,zc,r0,delta_x))*outm
    do ivar=ndim+3,nvar
           uold(ind_cell,ivar)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,ivar)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,ivar)
    enddo  

#endif       

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


subroutine reset_amb(xc1,yc1,zc1,rc1,xc2,yc2,zc2,rc2,ind_cell,ilevel,pex,pey,pez)
  use hydro_commons
  use amr_parameters
  use amr_commons !for dtnew
  use hydro_parameters
  use wind_params
  implicit none
  real(dp)::f,pex,pey,pez,f3d
  real(dp)::vflu,vx,vy,vz
  real(dp)::xc1,yc1,zc1,rc1
  real(dp)::xc2,yc2,zc2,rc2
  real(dp)::v1,v2,rho1,rho2
  integer::ind_cell,ilevel
  real(dp)::a1=1.0d0
  real(dp)::c1,c2,fx1,fy1,fx2,fy2,fz1,fz2

   v1 = vinf1*(1.0d0-rstar1)**beta1 
 rho1 = Mdot1*a1**2/(4.0d0*3.1415*v1)

   v2 = vinf2*(1.0d0-rstar2)**beta2 
 rho2 = Mdot2*a1**2/(4.0d0*3.1415*v2)

#if NDIM==2
 rho1 = Mdot1*a1/(2.0d0*3.1415*v1)
 rho2 = Mdot2*a1/(2.0d0*3.1415*v2)
#endif


#if NDIM==2
 fx1=  f1*f(vinf1,rstar1,rho1,v1,rc1)*xc1/rc1
 fy1=  f1*f(vinf1,rstar1,rho1,v1,rc1)*yc1/rc1
 fz1=  0
 fx2=  f2*f(vinf2,rstar2,rho2,v2,rc2)*xc2/rc2
 fy2=  f2*f(vinf2,rstar2,rho2,v2,rc2)*yc2/rc2
 fz2=  0

#endif
#if NDIM==3
 fx1=  f1*f3d(vinf1,rstar1,rho1,v1,rc1)*xc1/rc1
 fy1=  f1*f3d(vinf1,rstar1,rho1,v1,rc1)*yc1/rc1
 fz1=  f1*f3d(vinf1,rstar1,rho1,v1,rc1)*zc1/rc1
 fx2=  f2*f3d(vinf2,rstar2,rho2,v2,rc2)*xc2/rc2
 fy2=  f2*f3d(vinf2,rstar2,rho2,v2,rc2)*yc2/rc2
 fz2=  f2*f3d(vinf2,rstar2,rho2,v2,rc2)*zc2/rc2
#endif

!centered scheme for velocity
vx= 0.5d0*(uold(ind_cell,2)+pex)/uold(ind_cell,1)
vy= 0.5d0*(uold(ind_cell,3)+pey)/uold(ind_cell,1)
vz= 0.0d0

#if NDIM==3
   vz= 0.5d0*(uold(ind_cell,4)+pez)/uold(ind_cell,1)
#endif

  c1 = (fx1*vx+fy1*vy*fz1*vz)*dtnew(ilevel)
  c2 = (fx2*vx+fy2*vy*fz2*vz)*dtnew(ilevel)

!pressure
  uold(ind_cell,ndim+2)=uold(ind_cell,ndim+2)+c1+c2
  uold(ind_cell,2)=uold(ind_cell,2)+(fx1+fx2)*dtnew(ilevel)
  uold(ind_cell,3)=uold(ind_cell,3)+(fy1+fy2)*dtnew(ilevel)
#if NDIM==3
  uold(ind_cell,4)=uold(ind_cell,4)+(fz1+fz2)*dtnew(ilevel)
#endif


end subroutine reset_amb


subroutine init2_star(x,u,dx,nn,x0,y0,z0,r0,Mdot,vinf,Mach,rstar,beta,delta_x,index_pas_scal)
    
!This routine computes the hydro variables in the mask region, in the case of a star
!it returns the vector u and requires all the parameters of the star as input

  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  use wind_params
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
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables !!! save
  real(dp):: int_rho,int_P
  real(dp):: int_vx,int_vy,int_vz
  real(dp):: delta_x
  real(dp):: boxlen_0
  real(dp):: Mdot,Mach
  real(dp)::v_a,rho_a,P_a
  real(dp)::a1=1.0d0
  integer::index_pas_scal

  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  v_a = vinf*(1.0d0-rstar)**beta 
  rho_a = Mdot*a1**2/(4*3.1415d0*v_a)
  
#if NDIM==2
  rho_a = Mdot*a1/(2*3.1415d0*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach**2*gamma)

  zz=0.0d0 
  
  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0) THEN

        q(i,iu) =  int_vx(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,iv) =  int_vy(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,iw) =  int_vz(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,id) =  int_rho(xx,yy,zz,Mdot,r0,delta_x,vinf,rstar,beta)
        q(i,ip) =  int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,delta_x)
        if (nvar .gt. ip ) then
           q(i,ip+1)=0.0d0
           q(i,ip+2)=0.0d0
           q(i,ip+index_pas_scal)= 1.0d0
        endif
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
        q(i,iu) = int_vx(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,iv) = int_vy(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,id) = int_rho(xx,yy,zz,Mdot,r0,delta_x,vinf,rstar,beta)
        q(i,ip) = int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,delta_x)
        if (nvar .gt. ip) then
           q(i,ip+1)=0.0d0
           q(i,ip+2)=0.0d0
           q(i,ip+index_pas_scal) = 1.0d0
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



subroutine init2_puls(x,u,dx,nn,x0,y0,z0,r0,Mdot,vinf,Mach,delta_x,index_pas_scal)
 
  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  use wind_params 
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
  real(dp)::a1=1.0d0
  integer ::index_pas_scal

 
  IF (nn .eq. 0) RETURN
!faudra voir si tous les use sont utiles

  rho_a = Mdot*a1**2/(4.0d0*3.141592*vinf)

#if NDIM==2
  rho_a = Mdot*a1/(2.0d0*3.141592*vinf)
#endif

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
        q(i,id) = int_rhoP(xx,yy,zz,Mdot,r0,delta_x,vinf)
        q(i,ip) = int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,delta_x)
        if (nvar .gt. ip) then
           q(i,ip+1)=0.0d0
           q(i,ip+2)=0.0d0
           q(i,ip+index_pas_scal)=1.0d0
        endif
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
         q(i,id) = int_rhoP(xx,yy,zz,Mdot,r0,delta_x,vinf)
         q(i,ip) = int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,delta_x)
         if (nvar .gt. ip) then
            q(i,ip+1)=0.0
            q(i,ip+2)=0.0
            q(i,ip+index_pas_scal)= 1.0d0
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


end subroutine init2_puls


function int_rho(xx,yy,zz,Mdot,r0,dx,vinf,rstar,beta)
  use amr_parameters
implicit none


real(kind=8)::int_rho
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::rho
real(kind=8)::int
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::vol
real(kind=8)::Mdot,vinf,rstar,beta
real(kind=8)::r0 
real(kind=8)::dx


vol = dx**ndim

int =0.0d0

#if NDIM==2
   do i=1,n
       do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         int = int+(dx/n)**ndim*rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),zz,Mdot,r0,vinf,rstar,beta)
             end do
   end do  
#endif

#if NDIM==3
   do i=1,n
       do j=1,n
          do k=1,n
         int = int+(dx/n)**ndim*rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),(zz+(k-5.5d0)*dx/n),Mdot,r0,vinf,rstar,beta)
          end do
       end do
   end do
#endif

int_rho=int/vol

return
end function int_rho

function rho(x,y,z,Mdot,r0,vinf,rstar,beta) 
  use amr_parameters
implicit none

real(kind=8) :: rho,v
real(kind=8)::x,y,z
real(kind=8)::r0 
real(kind=8)::r 
real(kind=8)::Mdot,vinf,rstar,beta

r=sqrt(x**2+y**2+z**2)
 
#if NDIM==2
      if (r==0.0d0)then
         r=rstar/100d0
      endif
        rho=Mdot/(2.0d0*3.1415d0*r*v(r,vinf,rstar,beta))
#endif
 
#if NDIM==3
     if (r==0.0d0)then
         r=rstar/100d0
     endif
     rho=Mdot/(4.0d0*3.1415d0*r**2*v(r,vinf,rstar,beta))
#endif

return
end function rho


function int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,dx)
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
real(kind=8)::gamma,Mdot,vinf,rstar,beta
real(kind=8)::dx,r0

vol = dx**ndim

int =0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         rho_r = rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),zz,Mdot,r0,vinf,rstar,beta)
         int =int +(dx/n)**ndim*P_a*(rho_a/rho_r)**(-gamma) 
      end do
   end do
#endif

#if NDIM==3
    do k=1,n
       do i=1,n
          do j=1,n
             rho_r=rho((xx+(i-5.5d0)*dx/n),(yy+(j-5.5d0)*dx/n),(zz+(k-5.5d0)*dx/n),Mdot,r0,vinf,rstar,beta)
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


function int_rhoP(xx,yy,zz,Mdot,r0,dx,vinf)
  use amr_parameters
implicit none
real(kind=8)::int_rhoP
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::int
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::vol
real(kind=8)::Mdot,vinf
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
          int= int+(dx/n)**ndim*(1.0d0/r)
      end do
   end do  
int=int*Mdot/(2*3.1415*vinf)


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
         end do
       end do
   end do
int=int*Mdot/(4*3.1415*vinf)
#endif

int_rhoP=int/vol
return
end function int_rhoP


function int_Pp(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,dx)
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
real(kind=8)::gamma,Mdot,vinf
real(kind=8)::dx,r0!,P_r0
real(kind=8)::r!,rho_r
real(kind=8)::a1=1.0d0

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
              int =int +(dx/n)**ndim*P_a*(min((a1/r)**(gamma),10*(a1/r0)**gamma))
                
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
                 int =int +(dx/n)**ndim*P_a*(min((a1/r)**(2*gamma),10*(a1/r0)**(2*gamma)))
          end do
       end do
    end do
#endif

int_pP=int/vol


return
end function int_pP

function f(vinf,rstar,rho_a,v_a,rc)
implicit none
real(kind=8)::f
real(kind=8)::vinf
real(kind=8)::rstar
real(kind=8)::rho_a
real(kind=8)::v_a
real(kind=8)::rc

f=(vinf*v_a*rho_a*rstar)/(rc**3)

return
end function f

function f3d(vinf,rstar,rho_a,v_a,rc)
implicit none
real(kind=8)::f3d
real(kind=8)::vinf
real(kind=8)::rstar
real(kind=8)::rho_a
real(kind=8)::v_a
real(kind=8)::rc

f3d=(vinf*v_a*rho_a*rstar)/(rc**4)

return
end function f3d

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
use wind_params
use amr_parameters ! pour boxlen

implicit none

 logical::nml_ok

 namelist/first_body/type1,Mdot1,Mach1,vinf1,r01,rstar1,beta1,f1
 namelist/second_body/type2,Mdot2,Mach2,vinf2,r02,rstar2,beta2,f2
 namelist/ambient_medium/a,x01,y01,z01,rotation,R_P,R_rho,exc,M,Gstar

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
  read(10,NML=first_body)
  rewind(10)
  read(10,NML=second_body)
  rewind(10)
  read(10,NML=ambient_medium)
  close(10)


!conversion according to the value of a
  vinf1=vinf1/a
  vinf2=vinf2/a
  
  r01  = r01/a
  r02  = r02/a

  rstar1 = rstar1/a
  rstar2 = rstar2/a
 
!adaptation to the size of the box
  x01 = x01*boxlen
  y01 = y01*boxlen
  z01 = z01*boxlen

 Gstar=Gstar/a**3

end subroutine  read_wind_params


subroutine leapfrog(x,y,vx,vy,det,ilevel)
use wind_params
implicit none

real(kind=8)::x,y,vx,vy,xf,yf,vxf,vyf
real(kind=8)::ax,ay
real(kind=8)::det
real(kind=8)::r
integer::ilevel
!write(*,*),x,y,vx,vy,det,'leapfrog'
write(*,*),det,ilevel

   xf=x+vx*det
   yf=y+vy*det
   
   r=sqrt((xf-x01)**2+(yf-y01)**2)
   
   ax=-Gstar*M/(r**2)*(xf-x01)/r
   ay=-Gstar*M/(r**2)*(yf-y01)/r

   vxf=vx+ax*det
   vyf=vy+ay*det
   
   x=xf
   y=yf
   
   vx=vxf
   vy=vyf

end subroutine leapfrog
