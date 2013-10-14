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

  if (disk_presence .eqv..true.) then
    DO i=1,nn ! all the cells
       xx=x(i,1)-x01 ! position cell/center
       yy=x(i,2)-y01
       zz=0.0d0
#if NDIM==3
       zz=x(i,3)-z01
#endif
       rc=sqrt(xx**2+yy**2)
       r3D=sqrt(xx**2+yy**2+zz**2)
       H=H_0*(rc/r_0)**((3.-qu)/2.)
       q(i,id)=max(rho_0*(r_0/rc)**pe*exp(-zz**2/2.d0/H**2),rho_amb)
       csound=cs_0*(r_0/rc)**(qu/2.0d0)
       q(i,ip)=csound**2*q(i,id)/gamma

!   velocity  not strictly keplerian due to pressure
       omega=sqrt(GM/(r3d*(rmask1**2+r3d**2)))
#if NDIM==2
       omega=sqrt(GM/(rc*(rmask1**2+rc**2))+csound**2*(pe+qu)/(rc**2)) 
#endif
       q(i,iu)=-omega*yy
       q(i,iv)=omega*xx
       if ((q(i,iu)*xx+q(i,iv)*yy).gt. 1.0e-13) then
         write(*,*)q(i,iu)*xx+q(i,iv)*yy,'pbe de vitesse radiale ini'
       endif
       
    ENDDO
  endif  

  call init_star(q,x,nn,vinf1,rstar1,Mdot1,Mach1,x01,y01,z01,rmask1,beta1)


!positions of the 2nd body
  x02   = x01+a1*(1.0d0-exc)
  y02   = y01
  z02   = z01

!initial velocity field (necessary for the position update)
!  vx02=0.0d0
!  vy02=sqrt(GM*(2.0D0/(a1-exc)-1.0d0/a1))

  if (type2 .eq. 'puls')then
   call init_puls(q,x,nn,vinf2,rmask2,Mdot2,Mach2,x02,y02,z02)
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


  if (rotation .eqv. .false.)then
     x02 = x01 +a1             
     y02 = y01 
     z02 = z01 
  else
! the softening length (=rmask) has to be taken into acount in the rotation speed
     rot=sqrt(GM/(a*(a**2+rmask1**2)))
     x02=x01+a1*cos(rot*t)
     y02=y01+a1*sin(rot*t)
 endif
!circular orbit only
!     if (nstep == 1) then
!        x02   = x01+a1*(1.0d0-exc)
!        y02   = y01
!        vx02=0.0d0
!        vy02=sqrt(GM*(2.0D0/(a1-exc)-1.0d0/a1))
!        r=sqrt((x02-x01)**2+(y02-y01)**2)
!        ax02=-GM/r**2*(x02-x01)/r
!        ay02=-GM/r**2*(y02-y01)/r
!        vx02=vx02+ax02*dtnew(ilevel)/2.0d0
!        vy02=vy02+ay02*dtnew(ilevel)/2.0d0
!     else  
 !    write(*,*),dtnew(ilevel)
  !   call leapfrog(x02,y02,vx02,vy02,dtnew(nlevelmax))
!write(*,*),x02,y02,vx02,vy02,dtnew(ilevel),'main'
!endif
!endif
   z02 = z01
 
   zc2 = 0.0d0
   zc1 = 0.0d0 

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


    call init2_star(xx,u_mask,dx_loc,ngrid,x01,y01,z01,rmask1,Mdot1,vinf1,Mach1,rstar1,beta1)

     if (type2 .eq. 'puls') then 
         call init2_puls(xx,u_mask,dx_loc,ngrid,x02,y02,z02,rmask2,Mdot2,vinf2,Mach2,dx)
     endif    

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
           call reset_star(xc1,yc1,zc1,rmask1,dx_loc,u_mask,i,ind_cell(i),ilevel,vinf1,a,rstar1,rc1)
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
           write(*,*),'wind'
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
  !namelist/first_body/type1,Mdot1,Mach1,vinf1,r01,rstar1,beta1,f1
 namelist/second_body/rstar2,type2,Mdot2,Mach2,vinf2,rmask2,beta2,f2
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
  write(*,*),disk_presence,rho_0,pe,qu,GM,x01,y01,z01,r_0,rmask1
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



subroutine leapfrog(x,y,vx,vy,det)
use wind_parameters
implicit none

real(kind=8)::x,y,vx,vy,xf,yf,vxf,vyf
real(kind=8)::ax,ay
real(kind=8)::det
real(kind=8)::r
!write(*,*),x,y,vx,vy,det,'leapfrog'

   xf=x+vx*det
   yf=y+vy*det
   
   r=sqrt((xf-x01)**2+(yf-y01)**2)

   ax=-GM/(r**2)*(xf-x01)/r
   ay=-GM/(r**2)*(yf-y01)/r

   vxf=vx+ax*det
   vyf=vy+ay*det
   
   x=xf
   y=yf
   
   vx=vxf
   vy=vyf

end subroutine leapfrog
