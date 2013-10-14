!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  use disk_parameters
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
  real::H,cs,v,Mach,omega,r3d


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

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     zz=0.0d0
#if NDIM==3
     zz=x(i,3)-z0
#endif
     rc=sqrt(xx**2+yy**2)
     r3D=sqrt(xx**2+yy**2+zz**2)
     H=H_0*(rc/r_0)**((3.-qu)/2.)
     q(i,id)=max(rho_0*(r_0/rc)**pe*exp(-zz**2/2.d0/H**2),rho_amb)
     cs=cs_0*(r_0/rc)**(qu/2.0d0)
     q(i,ip)=cs**2*q(i,id)
!the velocity is not strictly keplerian due to pressure
!     omega=sqrt(G*M/(rc)**3+cs**2*pe/(r_0*rc))
     omega=sqrt(G*M/(rmask**2+r3d**2)**(3./2.))
     v=omega*r3d 
     q(i,iu)=v*yy/rc
     q(i,iv)=-v*xx/rc
 ENDDO     

   call init_star(q,x,nn,vinf,rstar,Mdot,Mach,x0,y0,z0,rmask)
     
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



subroutine init_star(q,x,nn,vinf1,rstar1,Mdot1,Mach1,x01,y01,z01,r01)
  use disk_parameters
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
   real(dp)::rho_a,P_A,v_a  ! variables for r =a
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables!!!!
  real(dp):: int_rho,int_p,int_vx,int_vy,int_vz ! function which computes the mean value of a variable in the mask
  real(dp):: frac
  real(dp)::pi=3.141459265
   integer:: i,j,id,iu,iv,iw,ip
  real::H,cs,Mach,density,v,R3D
  real(dp)::vinf1,rstar1,Mdot1,Mach1,x01,y01,z01,r01,vwind,vrot
 
 id=1 ; iu=2 ; iv=3; iw=4; ip=ndim+2

  dx=boxlen/2.0d0**(levelmin)
 !compute variables at 'a'
!   v_a = vinf*(1.0d0-rstar)
!  rho_a = Mdot/(4.0d0*pi*a**2*v_a)
!#if NDIM==2
! rho_a = Mdot/(2.0d0*pi*a*v_a)
!#endif
  zz = 0.0d0  

  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/center
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /center

#if NDIM == 3 
     zz=x(i,3)-z0
     r3D=sqrt(xx**2+yy**2+zz**2)
     cs=cs_0*(r_0/r01)**(qu/2.)
     IF (r3d .lt. 1.2d0*r01) THEN
        H=H_0*(rc/r_0)**((3.-qu)/2.)
        q(i,id)=rho_0/100*(r_0/r01)**pe
        q(i,ip)=cs**2*q(i,id)
        v=r3D*sqrt(G*M/(r3D**2+r01**2)**(3./2.))
        vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
        vrot=v*yy/rc
        q(i,iu) =  frac(xx,yy,zz,r01,dx)* (vwind+vrot)+(1.d0-frac(xx,yy,zz,r01,dx))*q(i,iu)
        vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
        vrot=-v*xx/rc
        q(i,iv) =  frac(xx,yy,zz,r01,dx)* (vwind+vrot)+(1.d0-frac(xx,yy,zz,r01,dx))*q(i,iv)
        q(i,iw) =  frac(xx,yy,zz,r01,dx)* int_vz(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r01,dx))*q(i,iw)

     !   q(i,id) =  frac(xx,yy,zz,r0,dx)*int_rho(xx,yy,zz,Mdot,r0,dx,a,vinf,rstar,beta)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,id)
     !   q(i,ip) =  frac(xx,yy,zz,r0,dx)*int_P(xx,yy,zz,P_a,rho_a,Mdot,a,r0,vinf,rstar,beta,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,ip
     ENDIF
#endif


#if NDIM==2   
     cs=cs_0*(r_0/r01)**(qu/2.)
     IF (rc .lt. 1.2d0*r01) THEN
        q(i,id)=rho_0*(r_0/r01)**pe
        q(i,ip)=cs**2*q(i,id)
        v=rc*sqrt(G*M/(rc**2+r01**2)**(3./2.))
        q(i,iu)=v*yy/rc
        q(i,iv)=-v*xx/rc
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
  use disk_parameters
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
  real(dp)::v,rho
  real(dp)::delta_x
  
  ncell=ncoarse+twotondim*ngridmax
 
!values for r=a (mettre dans une fonction)
   v = vinf*(1.0d0-rstar)
 rho = Mdot/(4.0d0*3.1415*a**2*v)


#if NDIM==2
 rho = Mdot/(2.0d0*3.1415*a*v)
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

     call init2_star(xx,u_mask,dx_loc,ngrid)

     do i=1,ngrid
         xc1=xx(i,1)-x0 
         yc1=xx(i,2)-y0
         rc1=sqrt(xc1**2+yc1**2)
#if NDIM ==3 
         zc1=xx(i,3)-z0
         rc1=sqrt(xc1**2+yc1**2+zc1**2) 
#endif
         if (rc1 .lt. 1.2*rmask) then
             call reset_star(xc1,yc1,zc1,rmask,dx,u_mask,i,ind_cell(i),ilevel,vinf,rstar,rho,v,rc1)
                       
!ambient medium!!! radiation pressure.
         endif    
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask

subroutine reset_star(xc,yc,zc,r01,delta_x,u_mask,i,ind_cell,ilevel,vinf1,a1,rstar1,rho,v,rc)
                                       
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
  real(dp)::xc,yc,zc,delta_x,rc,r01,vinf1,a1,rstar1,rho,v
  integer ::ind_cell    

  do j=1,nvar
      uold(ind_cell,j)=frac(xc,yc,zc,r01,delta_x)*u_mask(i,j)+(1.d0-frac(xc,yc,zc,r01,delta_x))*uold(ind_cell,j)
  enddo

end subroutine reset_star

subroutine init2_star(x,u,dx,nn)
    
!This routine computes the hydro variables in the mask region, in the case of a star
!it returns the vector u and requires all the parameters of the star as input
  use disk_parameters
  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  implicit none

  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer:: i,j,id,iu,iv,iw,ip
  real(dp)::rc,rs,xx,yy,zz
   integer::  ivar
   real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables !!! save
  real(dp):: int_rho,int_P,H,density,cs
  real(dp):: int_vx,int_vy,int_vz
  real(dp):: delta_x
  real(dp):: boxlen_0
  real(dp)::v_a,rho_a,P_a,Mach,v,r3d,vrot,vwind

  
  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

 ! v_a = vinf*(1.0d0-rstar)
 ! rho_a = Mdot/(4*3.1415d0*a**2*v_a)
  
!#if NDIM==2
!  rho_a = Mdot/(2*3.1415d0*a*v_a)
!#endif

  zz=0.0d0 
  
  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     r3d=sqrt(xx**2+yy**2+zz**2)
     rc=sqrt(xx**2+yy**2)
     cs=cs_0*(r_0/rmask)**(qu/2.)
     IF (r3d .lt. 1.2d0*rmask) THEN
          H=H_0*(rc/r_0)**((3.-qu)/2.)
          q(i,id)=rho_0/100.*(r_0/rmask)**pe
          q(i,ip)=cs**2*q(i,id)
          v=r3d*sqrt(G*M/(r3d**2+rmask**2)**(3./2.))
          vwind=int_vx(xx,yy,zz,vinf,rstar,beta,dx)
          vrot=v*yy/rc
          q(i,iu) = vwind+vrot
          vwind=int_vy(xx,yy,zz,vinf,rstar,beta,dx)
          vrot=-v*xx/rc
          q(i,iv) =  vwind+vrot
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

           cs=cs_0*(r_0/rmask)**(qu/2.)
#if NDIM==2   
        IF (rc .lt. 1.2d0*rmask) THEN
           q(i,id)=rho_0*(r_0/rmask)**pe
           q(i,ip)=cs**2*q(i,id)
           v=rc*sqrt(G*M/(rc**2+rmask**2)**(3./2.))
           q(i,iu)=v*yy/rc
           q(i,iv)=-v*xx/rc
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

v=max((vinf*(1.0d0-rstar/r)**beta),vinf/10)

return
end function v

