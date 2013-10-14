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
  integer :: ilevel
  integer:: i,j,id,iu,iv,iw,ip
  real(dp):: x01,y01,z01,x02,y02,z02,rc,rs,xx,yy,zz,r0,rcut
  real(dp):: rho_amb,P_amb
  real(dp):: R_star1,P_01,rho_01,v_w1,v_inf1
  real(dp):: R_star2,P_02,rho_02,v_w2,v_inf2
  integer :: ivar
  real(dp):: eta1, nu1, alpha
  real(dp):: eta2, nu2 
  real(dp):: eta_w,beta
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: integration2D_rho,integration2D_p,integration2D_vx,integration2D_vy ! function which computes the mean value of a variable
  real(dp):: fract
  real(dp):: Mdot_1,Mdot_2 ! mass loss rate
  real(dp):: boxlen_0 ! size of a small box
  real(dp):: delta_x

  IF (nn .eq. 0) RETURN


  id=1; iu=2; iv=3; iw=4; ip=ndim+2


  delta_x=boxlen/2.0d0**(levelmin)

!run parameters
   
  eta1   = 0.065d0
  nu1    = 0.85d0
  alpha =  0.69d0

  boxlen_0 = alpha*2.5d0
  r0    = 0.078d0*boxlen_0/boxlen

  x01   = 0.5d0*boxlen!-0.4d0*boxlen_0/boxlen)
  y01   = 0.5d0*boxlen
  z01   = 0.5d0*boxlen
  
 

  R_star1 = 0.0539*boxlen_0/boxlen


  beta = 1.0

  P_amb=26.4d0
  rho_amb=0.1d0

  P_01=0.358d0
  rho_01=540.d0

  v_w1=nu1
  v_inf1 = 1.36*v_w1



  Mdot_1=1.0d0/(4.d0*3.14d0*rho_01*v_w1*r0**2)
 



!computing eta_w=v1*Mdot1/v2*Mdot2. The wind should only depend on this parametre

!(eta1*nu1)/(eta2*nu2)
  
!  write(*,*),'eta wind vaut',eta_w


 ! There is an initial domain filled with a densité rho_amb and sound speed cs_amb. Afterwards, we add a mask (for r <2r0) with a given densité, pressure and speed. 



q(:,id) = rho_amb
q(:,iu) = 0.0
q(:,iv) = 0.0
q(:,ip) = P_amb



! First star!!!!!!!!!!!


  DO i=1,nn ! all the cells
     xx=x(i,1)-x01 ! position cell/center
     yy=x(i,2)-y01
     rc=sqrt(xx**2+yy**2)! 2d distance /center


#if NDIM == 3
     zz=x(i,3)-z01
     rc3=sqrt(xx**2+yy**2+zz*2)
     IF (rc3 .lt. 2*r01) THEN
        q(i,id) = min(rho_0*(r01/rc)**2,4.d0*rho_0)
        q(i,iu) = v_inf*(xx/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,iv) = v_inf*(yy/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,iw) = v_inf*(zz/rc3)*max((1-R_star1/rc3)**beta,0d0)
        q(i,ip) = min(P_0*(r01/rc)**(10/3),2.d0**(2.d0*gamma)*P_0)
      ENDIF
#endif

#if NDIM==2   
        !write(*,*)'on est a 2D'
        IF (rc .lt. 2.0d0*r0) THEN

          ! q(i,id) =  fract(xx,yy,r0,delta_x)*integration2D_rho(xx,yy,rho_01,r0,delta_x)+(1.d0-fract(xx,yy,r0,delta_x))*q(i,id)
      !     write(*,*),fract(xx,yy,r0,delta_x),'frac rho'
          ! q(i,iu) =  fract(xx,yy,R_star1,delta_x)*integration2D_vx(xx,yy,v_inf1,r0,delta_x)+(1.d0-fract(xx,yy,r0,delta_x))*q(i,iu)
  !  write(*,*),fract(xx,yy,r0,delta_x),'frac vx'
          ! q(i,iv) =  fract(xx,yy,R_star1,delta_x)*integration2D_vy(xx,yy,v_inf1,r0,delta_x)+(1.d0-fract(xx,yy,r0,delta_x))*q(i,iv)
! write(*,*),fract(xx,yy,r0,delta_x),'frac vy'
          ! q(i,ip) = fract(xx,yy,r0,delta_x)*integration2D_P(xx,yy,P_01,r0,delta_x)+(1.d0-fract(xx,yy,r0,delta_x))*q(i,ip)
!write(*,*),fract(xx,yy,r0,delta_x),'frac vz'
!
           q(i,id) =  fract(xx,yy,r0,delta_x)*integration2D_rho(xx,yy,rho_01,r0,delta_x)&
                +(1.d0-fract(xx,yy,r0,delta_x))*q(i,id)
           q(i,iu) =  fract(xx,yy,r0,delta_x)*integration2D_vx(xx,yy,v_inf1,R_star1,delta_x)&
                +(1.d0-fract(xx,yy,r0,delta_x))*q(i,iu)
           q(i,iv) =  fract(xx,yy,r0,delta_x)*integration2D_vy(xx,yy,v_inf1,R_star1,delta_x)&
                +(1.d0-fract(xx,yy,r0,delta_x))*q(i,iv)
           q(i,ip) = fract(xx,yy,r0,delta_x)*integration2D_P(xx,yy,P_01,r0,delta_x)&
                +(1.d0-fract(xx,yy,r0,delta_x))*q(i,ip)


 !          q(i,id) = min(rho_01*(r0/rc)**2,4.d0*rho_01)
 !          q(i,iu) = v_inf1*(xx/rc)*max((1-R_star1/rc)**beta,0.d0)
 !          q(i,iv) = v_inf1*(yy/rc)*max((1-R_star1/rc)**beta,0.d0)
 !          q(i,ip) = min(P_01*(r0/rc)**(10/3),2.d0**(2.d0*gamma)*P_01)
!           q(i,id)=integration2D_rho(xx,yy,rho_01,r0,delta_x)
!           q(i,iu)=integration2D_vx(xx,yy,v_inf1,R_star1,delta_x)
!           q(i,iv)=integration2D_vy(xx,yy,v_inf1,R_star1,delta_x)
!           q(i,ip)=integration2D_P(xx,yy,P_01,r0,delta_x)

        ENDIF
#endif
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



subroutine reset_mask(ilevel) ! Loop over grids by vector sweeps
! routine which reinitializes the wind at each timestep

  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif


! trier les variables inutiles!!!!!!!!!!!!!!!!!
  integer::ilevel
  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
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
  real(dp)::x01,y01,z01,rc,xcen,ycen,r0,x02,y02,z02
  real(dp)::fract
  real(dp)::P_01,P_02,rho_01,rho_02,v_w1,v_w2,v_inf1,v_inf2
  real(dp):: eta1, nu1
  real(dp):: eta2, nu2
  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar

  real(dp)::boxlen_0 !size of small box
  real(dp)::delta_x!
  
  delta_x=boxlen/2.0d0**ilevel ! size of a cell     

  !write(*,*),delta_x,'deltax'

!a ameliorer!!!!!!!!!
   eta1   = 0.0650
   nu1    = 0.8500
   alpha = 0.69
 
   boxlen_0 = alpha*2.5d0

   v_w1=nu1
   v_inf1 = 1.36*v_w1
   P_01=0.358d0
   rho_01=540.d0

  r0    = 0.078*boxlen_0/boxlen   

  x01   = 0.5d0*boxlen!-0.4d0*boxlen_0/boxlen)
  y01   = 0.5d0*boxlen
  z01   = 0.5d0*boxlen
    
!    write(*,*),x01,x02,'x01,x02'

! !!! copie de trucs dans init_flow_fine, c'est peut etre inutile mais faut voir
 ! Conversion factor from user units to cgs units

! à enlever????????????
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
     do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
     end do
           ! Rescale position from code units to user units
     do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
     end do
    ! Call reset of the mask routine, changes uold in the mask
     call condinit2(xx,u_mask,dx_loc,ngrid,x01,y01,z01,r0,rho_01,P_01,v_inf1,alpha,delta_x)
  !   call condinit2(xx,u_mask,dx_loc,ngrid,x02,y02,z02,r0,rho_02,P_02,v_inf2,alpha,delta_x)
     


! mettre 3D!
     do ivar=1,nvar
         do i=1,ngrid
            xcen=xx(i,1)-x01 ! position cell/centre          !!!
            ycen=xx(i,2)-y01
            rc=sqrt(xcen**2+ycen**2)! 2d distance /centre          !!!
            if (rc .lt. 2*r0) then
            !   uold(ind_cell(i),ivar)=u_mask(i,ivar)
               uold(ind_cell(i),ivar)=fract(xcen,ycen,r0,delta_x)*u_mask(i,ivar)&
                    +(1.d0-fract(xcen,ycen,r0,delta_x))*uold(ind_cell(i),ivar)

            endif
         end do
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask


subroutine condinit2(x,u,dx,nn,x0,y0,z0,r0,rho_0,P_0,v_inf,alpha,delta_x)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use hydro_commons  !!!sur?
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
  real(dp):: R_star, beta,rho_amb,P_amb,P_0,rho_0
  integer::  ivar
  real(dp):: eta,nu, alpha, v_w,v_inf
  !real(dp):: coef_x0,coef_y0,coef_z0,coef_r0
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: integration2D_rho,integration2D_P
  real(dp):: integration2D_vx,integration2D_vy
  real(dp):: delta_x
  real(dp):: boxlen_0

  IF (nn .eq. 0) RETURN


! faire le cas 1D!3D! Passer les variables hydro en argument

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

!run parameters

  boxlen_0=alpha*2.5d0

  R_star = 0.054d0*boxlen_0/boxlen

  !x0    = 0.5d0*(1.0d0-coef_x0*boxlen_0/boxlen)! position of the centre of the mask
  !y0    = 0.5d0*(1.0d0-coef_y0*boxlen_0/boxlen)
  !z0    = 0.5d0*(1.0d0-coef_z0*boxlen_0/boxlen)
  !r0    = coef_r0*boxlen_0/boxlen
  beta = 1.0

!write(*,*),x0,'x0 condinit'
!write(*,*),r0,'r0 condinit'


!mettre les min et max!!!!!!!!!!!!
  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/centre
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /centre
#if NDIM == 3
     zz=x(i,3)-z0
     rc3=sqrt(xx**2+yy**2+zz*2)
     IF (rc3 .lt. 2*r0) THEN
        q(i,id) = rho_0*(r0/rc)**2
        q(i,iu) = v_inf*(xx/rc3)*(1-R_star/rc3)**beta
        q(i,iv) = v_inf*(yy/rc3)*(1-R_star/rc3)**beta
        q(i,iw) = v_inf*(zz/rc3)*(1-R_star/rc3)**beta
        q(i,ip) = P_0*(r0/rc)**(10/3)

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
        IF (rc .lt. 2*r0) THEN
           q(i,id) = integration2D_rho(xx,yy,rho_0,r0,delta_x)
      
           ! mettre la precision en argument!!!!!!!
           !q(i,id) = min(rho_0*(r0/rc)**2,4.d0*rho_0)
           q(i,iu) = integration2D_vx(xx,yy,v_inf,R_star,delta_x)
           q(i,iv) = integration2D_vy(xx,yy,v_inf,R_star,delta_x)
          ! q(i,iu) = v_inf*(xx/rc)*max((1-R_star/rc)**beta,0.d0)
          ! q(i,iv) = v_inf*(yy/rc)*max((1-R_star/rc)**beta,0.d0)
          !q(i,ip) = min(P_0*(r0/rc)**(10/3),2.d0**(2.d0*gamma)*P_0)
           q(i,ip) = integration2D_P(xx,yy,P_0,r0,delta_x)

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


end subroutine condinit2

!regrouper les 2 fonctions ensemble?????????????

function integration2D_rho(xx,yy,rho_0,r0,delta_x)
implicit none


real(kind=8)::integration2D_rho
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j
real(kind=8)::density2D
real(kind=8)::integral
real(kind=8)::xx,yy   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::rho_0
real(kind=8)::r0 
real(kind=8)::delta_x!=1.725/2**7

!write(*,*),delta_x,'delta x integration2D'
volume= delta_x*delta_x! (size of box/number of cells)**2, to improve according to number of dimensions

integral=0.0d0

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         integral =integral +delta_x/n*delta_x/n*density2D((xx+(i-5-0.5d0)*delta_x/n),(yy+(j-5-0.5d0)*delta_x/n),rho_0,r0) 
    end do
end do  

integration2D_rho=integral/volume

return
end function


function density2D(x,y,rho_0,r0) ! function which returns the density for a given position
implicit none

real(kind=8) :: density2D
real(kind=8)::x,y
real(kind=8)::rho_0
real(kind=8)::r0 
real(kind=8)::r ! distance from center of the mask


r=sqrt(x**2+y**2)

!takes into account the issue for r=0 and sets a threshold to density
if (r .lt. 0.5d0*r0)then
   density2D=min(rho_0*(r0/r)**2,4.d0*rho_0)
else 
   density2D=rho_0*(r0/r)**2
endif

return
end function



function integration2D_P(xx,yy,P_0,r0,delta_x)
implicit none


real(kind=8)::integration2D_P
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j
real(kind=8)::pressure2D
real(kind=8)::integral
real(kind=8)::xx,yy   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::P_0
real(kind=8)::r0
real(kind=8)::delta_x!=1.725/2**7

volume= delta_x*delta_x! (size of box/number of cells)**2, to improve according to number of dimensions

integral=0.0d0

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         integral =integral +delta_x/n*delta_x/n*pressure2D((xx+(i-5-0.5d0)*delta_x/n),(yy+(j-5-0.5d0)*delta_x/n),P_0,r0) 
    end do
end do  


integration2D_P=integral/volume

return
end function


function pressure2D(x,y,P_0,r0) ! function which returns the density for a given position
implicit none

real(kind=8) :: pressure2D
real(kind=8)::x,y
real(kind=8)::P_0 
real(kind=8)::r0
real(kind=8)::r ! distance from center of the mask
real(kind=8)::gamma =5.d0/3.d0 ! a ameliorer


r=sqrt(x**2+y**2)

if (r .lt. 0.5d0*r0)then
   pressure2D=min(P_0*(r0/r)**(10/3),2.0d0**(2.0d0*gamma)*P_0)
else
   pressure2D=P_0*(r0/r)**(10/3)
endif

return
end function


function integration2D_vx(xx,yy,v_inf,R_star,delta_x)
implicit none


real(kind=8)::integration2D_vx
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j
real(kind=8)::velocity2D_x
real(kind=8)::integral
real(kind=8)::xx,yy   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::v_inf
real(kind=8)::R_star
real(kind=8)::delta_x!=1.725/2**7

volume= delta_x*delta_x! (size of box/number of cells)**2, to improve according to number of dimensions

integral=0.0d0

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
      integral =integral +delta_x/n*delta_x/n*velocity2D_x((xx+(i-5-0.5d0)*delta_x/n),(yy+(j-5-0.5d0)*delta_x/n),v_inf,R_star)
      
    end do
end do  

integration2D_vx=integral/volume

return
end function integration2D_vx



function integration2D_vy(xx,yy,v_inf,R_star,delta_x)
implicit none


real(kind=8)::integration2D_vy
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j
real(kind=8)::velocity2D_y
real(kind=8)::integral
real(kind=8)::xx,yy   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::delta_x!=1.725d0/128.d0!,delta_y=1.725d0/128d0 ! size of the cell
real(kind=8)::v_inf
real(kind=8)::R_star

volume= delta_x*delta_x! (size of box/number of cells)**2, to improve according to number of dimensions

integral=0.0d0

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         integral =integral +delta_x/n*delta_x/n*velocity2D_y((xx+(i-5-0.5d0)*delta_x/n),(yy+(j-5-0.5d0)*delta_x/n),v_inf,R_star)
         
    end do
end do  

integration2D_vy=integral/volume

return
end function integration2D_vy




!faire un tableau pour avoir les vitesses sur les 3 axes (voire toutes les autres variables d'ailleurs

function velocity2D_x(x,y,v_inf,R_star) ! function which returns the velocity for a given position
implicit none
real(kind=8):: velocity2D_x
real(kind=8)::x,y
real(kind=8)::R_star 
real(kind=8)::v_inf  
real(kind=8)::r ! distance from center of the mask
real(kind=8)::beta =1.d0 ! a ameliorer


r=sqrt(x**2+y**2)
velocity2D_x=v_inf*(x/r)*max((1-R_star/r)**beta,0.d0)

return
end function velocity2D_x

function velocity2D_y(x,y,v_inf,R_star) ! function which returns the velocity for a given position
implicit none
real(kind=8):: velocity2D_y
real(kind=8)::x,y
real(kind=8)::R_star 
real(kind=8)::v_inf   
real(kind=8)::r ! distance from center of the mask
real(kind=8)::beta =1.d0 ! a ameliorer


r=sqrt(x**2+y**2)
velocity2D_y=v_inf*(y/r)*max((1-R_star/r)**beta,0.d0)

return
end function velocity2D_y


function fract(x,y,r0,delta_x)

!function which computes the fraction of a cell being occupied by the mask
implicit none
real(kind=8):: fract
real(kind=8):: delta_x!=1.725d0/128.d0
real(kind=8):: x,y,r
real(kind=8):: x_s,y_s! postion on the center of the slice
real(kind=8):: r0
integer :: i,j
integer :: n=10 ! number of slices in each cell in each direction
real :: compteur


compteur=0.d0

do i=1,n
   do j=1,n
      x_s=x+(i-5.0d0-0.5d0)*delta_x/n
      y_s=y+(j-5.0d0-0.5d0)*delta_x/n
      r=sqrt(x_s**2+y_s**2)
      if (r .lt. 2.0d0*r0)then
         compteur=compteur+1.0d0
      endif   
    end do
end do  

fract=compteur/(n**2)
return
end function fract
