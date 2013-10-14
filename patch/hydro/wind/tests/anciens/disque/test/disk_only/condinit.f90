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
  integer :: ivar,ncell
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real::H,cs,v,Mach,omega,r3d,vrotx,vroty,vrot
  real(dp)::testx,testy,csound


  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  q(:,id) = rho_amb
  q(:,iu) = 0.0d0
  q(:,iv) = 0.0d0
#if NDIM == 3
  q(:,iw) = 0.0d0
#endif
  q(:,ip) = P_amb

!creation of the disk!
 
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
       q(i,iu)=omega*yy
       q(i,iv)=-omega*xx
       if ((q(i,iu)*xx+q(i,iv)*yy).ne. 0) then
         write(*,*)q(i,iu)*xx+q(i,iv)*yy,'pbe de vitesse radiale ini'
       endif
       
    ENDDO
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

 !compute variables at 'a'
!   v_a = vinf*(1.0d0-rstar)
!  rho_a = Mdot/(4.0d0*pi*a**2*v_a)
!#if NDIM==2
! rho_a = Mdot/(2.0d0*pi*a*v_a)
!#endif

  !compute the density in the wind, at the edge of the mask
  !a=1.0d0 !(a ameliorer!!!)


  v_mask=vinf*(1.0d0-rstar/rmask)!**beta
  rho_0w=Mdot*a**2/(4.0d0*3.1415*(rmask)**2*v_mask)
#if NDIM==2
 rho_0w = Mdot*a/(2.0d0*3.1415*(a*rmask)*v_mask)
#endif

!  v_r0=vinf*(1.0d0-rstar/r_0)**beta
!  rho_0w=Mdot/(4.0d0*pi*(a*r_0)**2*v_r0)
!#if NDIM==2
! rho_0w = Mdot/(2.0d0*pi*(a*r_0)*v_r0)
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
        if (disk_presence .eqv. .true.) then
           omega=sqrt(GM/(r3d*(rmask**2+r3d**2)))!! no pressure gradient correction yet
        endif  
     H=H_0*(rc/r_0)**((3.-qu)/2.)
       ! q(i,id)=   frac(xx,yy,zz,rmask,dx)* rho_0w*(r_0/rmask)**pe +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id) 
     q(i,id)=   frac(xx,yy,zz,rmask,dx)* rho_0w +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id)    
     q(i,ip)=csound**2*q(i,id)/gamma
     vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
      !  vrot=vrotx(xx,yy,zz,GM,rmask,dx) !v*yy/rc
      !  q(i,iu) =  frac(xx,yy,zz,rmask,dx)* (vwind+vrot)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iu)
!        vrot=vroty(xx,yy,zz,GM,rmask,dx)! -v*xx/rc
     q(i,iu) =  frac(xx,yy,zz,rmask,dx)* (vwind+omega*yy)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iu)
     vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
     q(i,iv) =  frac(xx,yy,zz,rmask,dx)* (vwind-omega*xx)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iv)
     q(i,iw) =  frac(xx,yy,zz,rmask,dx)* int_vz(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,iw)
     ENDIF
#endif


#if NDIM==2   
     IF (rc .lt. 1.2d0*rmask) THEN
 ! In 2D there is either a disk, or a wind, there is no need having them both. Yet there has to be as mask. 
! The density has to be continuous across the edge of the disk.
        if (disk_presence .eqv..true.) then
           q(i,id)= frac(xx,yy,zz,rmask,dx)* rho_0*(r_0/rmask)**pe +(1.d0-frac(xx,yy,zz,rmask,dx))*q(i,id) 
           q(i,ip)=csound**2*q(i,id)/gamma
           omega=sqrt(GM/(rc*(rmask**2+rc**2))+csound**2*(pe+qu)/(rc**2)) 
           q(i,iu)=omega*yy
           q(i,iv)=-omega*xx
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
  integer::ilevel,j
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
  real(dp)::v,rho
  real(dp)::delta_x
  
  ncell=ncoarse+twotondim*ngridmax
 
!values for r=a (mettre dans une fonction)
   v = vinf1*(1.0d0-rstar1/rmask1)
 rho = Mdot1*a**2/(4.0d0*3.1415*(rmask1)**2*v)


#if NDIM==2
 rho = Mdot1*a/(2.0d0*3.1415**rmask1*v)
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
          !!!  delta_x = boxlen/2.0d0**ilevel  sert à rien?????
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
     end do
           ! Rescale position from code units to user units
     do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
     end do


     call init2_star(xx,u_mask,dx_loc,ngrid,x01,y01,z01,rmask1,Mdot1,vinf1,Mach1,rstar1,beta1)
!verifier delta_x!


     do i=1,ngrid
         xc1=xx(i,1)-x01 
         yc1=xx(i,2)-y01
         rc1=sqrt(xc1**2+yc1**2)
#if NDIM ==3 
         zc1=xx(i,3)-z01
         rc1=sqrt(xc1**2+yc1**2+zc1**2) 
#endif
                
         if (rc1 .lt. 1.2*rmask1) then
            
!!!! voir pout rho et v!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call reset_star(xc1,yc1,zc1,rmask1,dx_loc,u_mask,i,ind_cell(i),ilevel,vinf1,a,rstar1,rho,v,rc1)


!ambient medium!!! radiation pressure.
         endif    
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask

subroutine reset_star(xc,yc,zc,rmask,delta_x,u_mask,i,ind_cell,ilevel,vinf,a,rstar,rho,v,rc)
                                       
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
#if NDIM==2
 rho_0w = Mdot/(2.0d0*3.1415*(a*rmask)*v_rmask)
#endif





  zz = 0.0d0  

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center
     csound=cs_0*(r_0/rmask)**(qu/2.)
 ! v_a = vinf*(1.0d0-rstar)
 ! rho_a = Mdot/(4*3.1415d0*a**2*v_a)
  
!#if NDIM==2
!  rho_a = Mdot/(2*3.1415d0*a*v_a)
!#endif

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
        q(i,iu)=omega*yy+vwind
        vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
        q(i,iv)=-omega*xx+vwind
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
           q(i,iu)=omega*yy
           q(i,iv)=-omega*xx
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


function int_rho(xx,yy,zz,Mdot,r0,dx,a,vinf,rstar,beta)
  use amr_parameters
implicit none


real(kind=8)::int_rho
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::rho
real(kind=8)::int,a
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
      
      rho=Mdot/(2.0d0*3.1415d0*r*a*v(r,vinf,rstar,beta))
#endif
 
#if NDIM==3
     if (r==0.0d0)then
         r=rstar/100d0
     endif
      
     rho=Mdot/(4.0d0*3.1415d0*(r*a)**2*v(r,vinf,rstar,beta))
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

function testx(xx,yy,zz,GM,rmask,dx,csound,pe,qu)
 use amr_parameters
implicit none

real(kind=8)::testx
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::GM,omega,rmask
real(kind=8)::int
real(kind=8)::xx,yy,zz
real(kind=8)::vol,csound,pe,qu
real(kind=8)::dx ! size of the cell
real(kind=8)::r,r3d

vol = dx**ndim

int = 0.0d0
#if NDIM==2
   do i=1,n
      xx=xx+(i-5.5d0)*dx/n
      do j=1,n
         yy=yy+(j-5.5d0)*dx/n
         r=sqrt(xx**2+yy**2)
         if (r==0.0d0)then
         !       r=rstar/100d0   
         write(*,*),'r=0'
         endif
         omega=sqrt(GM/(r*(rmask**2+r**2))+csound**2*(pe+qu)/(r**2))
        ! if (xx/yy .lt.0)then
           int=int +(dx/n)**ndim*omega*yy
        ! else
        !    int =int +(dx/n)**ndim*omega*yy
        ! endif   
      end do
   end do 
#endif

#if NDIM==3
   do i=1,n
      xx=xx+(i-5.5d0)*dx/n
      do j=1,n
         yy=yy+(j-5.5d0)*dx/n
         do k=1,n
            zz=zz+(k-5.5d0)*dx/n
            r=sqrt(xx**2+yy**2)
            r3d=sqrt(xx**2+yy**2+zz**2)
            if (r==0.0d0)then
      !       r=rstar/100d0
             write(*,*),'r=0'
            endif
            omega=sqrt(GM/(r3d**2+rmask**2)**(3./2.))
  
        int =int +(dx/n)**ndim*omega*yy  
      end do
   end do 
   enddo
#endif

testx=int/vol
return

end function testx

function testy(xx,yy,zz,GM,rmask,dx,csound,pe,qu)
 use amr_parameters
implicit none

real(kind=8)::vroty,testy
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::omega,GM,rmask
real(kind=8)::int
real(kind=8)::xx,yy,zz
real(kind=8)::vol,csound,pe,qu
real(kind=8)::dx ! size of the cell
real(kind=8)::r,r3d

vol = dx**ndim

int = 0.0d0

#if NDIM==2
   do i=1,n
      xx=xx+(i-5.5d0)*dx/n
      do j=1,n
         yy=yy+(j-5.5d0)*dx/n
         r=sqrt(xx**2+yy**2)
         if (r==0.0d0)then
         !       r=rstar/100d0   
         write(*,*),'r=0'
         endif
         omega=sqrt(GM/(r*(rmask**2+r**2))+csound**2*(pe+qu)/(r**2))
        ! if (xx/yy .lt.0)then
        !    int =int -(dx/n)**ndim*omega*xx
        ! else
            int =int -(dx/n)**ndim*omega*xx 
        ! endif   
     end do
   end do 
#endif

#if NDIM==3
   do i=1,n
      xx=xx+(i-5.5d0)*dx/n
      do j=1,n
         yy=yy+(j-5.5d0)*dx/n
         do k=1,n
            zz=zz+(k-5.5d0)*dx/n
            r=sqrt(xx**2+yy**2)
            r3d=sqrt(xx**2+yy**2+zz**2)
            if (r==0.0d0)then
      !       r=rstar/100d0
             write(*,*),'r=0'
            endif
            omega=sqrt(GM/(r3d**2+rmask**2)**(3./2.))
            int =int -(dx/n)**ndim*omega*r3d*xx/r
          enddo
       end do
   end do 
#endif

testy=int/vol

return

end function testy


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


function v(r,vinf,rstar,beta) 
implicit none
real(kind=8):: v
real(kind=8):: rstar 
real(kind=8):: vinf   
real(kind=8):: r 
real(kind=8):: beta  

v=max((vinf*(1.0d0-rstar/r)**beta),vinf/1000)

return
end function v



subroutine read_wind_params(nml_ok)
!subroutine which reads the additionnal namelist
!use wind_params
use wind_parameters
use amr_parameters ! pour boxlen
use hydro_parameters

implicit none

 logical::nml_ok
 real(dp)::omega_0
 
namelist/disk/disk_presence,rho_0,pe,qu,GM,x01,y01,z01,rstar1,rmask1,vinf1,Mdot1

namelist/wind/GM,x01,y01,z01,rstar1,rmask1,vinf1,Mdot1,beta1,Mach1


! namelist/first_body/type1,Mdot1,Mach1,vinf1,r01,rstar1,beta1,f1
! namelist/second_body/type2,Mdot2,Mach2,vinf2,r02,rstar2,beta2,f2
! namelist/ambient_medium/a,x01,y01,z01,rotation,R_P,R_rho,exc,M,Gstar

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
  rewind(10)
  read(10,NML=wind)
  rewind(10)

!conversion according to the value of a
  vinf1=vinf1/a
!  vinf2=vinf2/a
  
  rmask1  = rmask1/a
 ! r02  = r02/a

  rstar1 = rstar1/a
 ! rstar2 = rstar2/a
 
!adaptation to the size of the box
  x01 = x01*boxlen
  y01 = y01*boxlen
  z01 = z01*boxlen

  GM=GM/a**3
  rho_0=rho_0*a**3

!determination of all the other parameters
  r_0=2*rmask1
  H_0=0.1*r_0
  omega_0=sqrt(GM/r_0**3)
  cs_0=H_0*omega_0
  p_0=cs_0**2*rho_0/gamma

  rho_amb=1.0e-8*rho_0
  P_amb=1.0e-8*P_0 
  

end subroutine  read_wind_params

