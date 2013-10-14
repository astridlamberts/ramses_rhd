!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn,coef_x01,coef_y01,coef_z01,coef_r01,coef_x02,coef_y02,coef_z02,coef_r02)
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
  integer:: i,j,id,iu,iv,iw,ip
  real(dp):: x01,y01,z01,x02,y02,z02,rc,rs,xx,yy,zz,r01,r02,rcut
  real(dp):: R_star1,R_star2, beta,rho_amb,P_amb,P0,rho0
  integer::  ivar
  real(dp):: eta,nu, alpha, v_w,v_inf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: coef_x01  ! position of the centre of the first star (scaled to the size of the box)
  real(dp):: coef_y01
  real(dp):: coef_z01
  real(dp):: coef_r01  ! radius of the mask of the first star


  real(dp):: coef_x02  ! position of the centre of the second star (scaled to the size of the box)
  real(dp):: coef_y02
  real(dp):: coef_z02
  real(dp):: coef_r02




  !coef_x0=0.5d0
  !coef_y0=0.5d0
  !coef_z0=0.5d0
  !coef_r0=0.195d0


  IF (nn .eq. 0) RETURN


  id=1; iu=2; iv=3; iw=4; ip=ndim+2

!run parameters
   eta = 0.065
   nu = 0.85
   alpha = 0.69

 ! x0   = coef_x0*boxlen
 ! y0   = coef_y0*boxlen
 ! z0   = coef_z0*boxlen
 ! r0   = coef_r0*alpha

  x01   = 0.35d0*boxlen
  y01   = 0.5d0*boxlen
  z01   = 0.5d0*boxlen
  r01   = 0.1d0*alpha !0.195
 

  x02   = 0.65d0*boxlen
  y02   = 0.5d0*boxlen
  z02   = 0.5d0*boxlen
  r02   = 0.1d0*alpha !0.195

  R_star1 = 0.0539 !(if alpha=0.69)
  R_star2 = 0.0539 !(if alpha=0.69)

  beta = 1.0
  rho_amb = 0.1*eta*alpha**(-2)*nu**(-1)
  P_amb = 45.d0*nu*eta*alpha**(-2)

  P0 = 361*nu*eta*alpha**(-2)
  rho0 = 540*eta*alpha**(-2)*nu**(-1)



  v_w=nu
  v_inf = 1.36*v_w




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
        q(i,id) = rho0*(r01/rc)**2
        q(i,iu) = v_inf*(xx/rc3)*(1-R_star1/rc3)**beta
        q(i,iv) = v_inf*(yy/rc3)*(1-R_star1/rc3)**beta
        q(i,iw) = v_inf*(zz/rc3)*(1-R_star1/rc3)**beta
        q(i,ip) = P0*(r01/rc)**(10/3)
      ENDIF
#endif

#if NDIM==2   
        !write(*,*)'on est a 2D'
        IF (rc .lt. 2*r01) THEN
           q(i,id) = min(rho0*(r01/rc)**2,4.d0*rho0)
           q(i,iu) = v_inf*(xx/rc)*max((1-R_star1/rc)**beta,0.d0)
           q(i,iv) = v_inf*(yy/rc)*max((1-R_star1/rc)**beta,0.d0)
           q(i,ip) = min(P0*(r01/rc)**(10/3),2.d0**(2.d0*gamma)*P0)
        ENDIF
#endif
   ENDDO     


! Second star :


  DO i=1,nn ! all the cells
     xx=x(i,1)-x02 ! position cell/center
     yy=x(i,2)-y02
     rc=sqrt(xx**2+yy**2)! 2d distance /center


#if NDIM == 3
     zz=x(i,3)-z02
     rc3=sqrt(xx**2+yy**2+zz*2)
     IF (rc3 .lt. 2*r02) THEN
        q(i,id) = rho0*(r01/rc)**2
        q(i,iu) = v_inf*(xx/rc3)*(1-R_star2/rc3)**beta
        q(i,iv) = v_inf*(yy/rc3)*(1-R_star2/rc3)**beta
        q(i,iw) = v_inf*(zz/rc3)*(1-R_star2/rc3)**beta
        q(i,ip) = P0*(r02/rc)**(10/3)
      ENDIF
#endif

#if NDIM==2   
        !write(*,*)'on est a 2D'
        IF (rc .lt. 2*r02) THEN
           q(i,id) = min(rho0*(r01/rc)**2,4.d0*rho0)
           q(i,iu) = v_inf*(xx/rc)*max((1-R_star2/rc)**beta,0.d0)
           q(i,iv) = v_inf*(yy/rc)*max((1-R_star2/rc)**beta,0.d0)
           q(i,ip) = min(P0*(r02/rc)**(10/3),2.d0**(2.d0*gamma)*P0)
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
  real(dp)::x01,y01,z01,rc,xcen,ycen,r01,x02,y02,z02,r02
  real(dp)::coef_x01=0.35d0
  real(dp)::coef_y01=0.5d0
  real(dp)::coef_z01=0.55d0
  real(dp)::coef_r01=0.1d0

  real(dp)::coef_x02=0.65d0
  real(dp)::coef_y02=0.5d0
  real(dp)::coef_z02=0.5d0
  real(dp)::coef_r02=0.1d0


  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar

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
     call condinit2(xx,u_mask,dx_loc,ngrid,coef_x01,coef_y01,coef_z01,coef_r01)
     call condinit2(xx,u_mask,dx_loc,ngrid,coef_x02,coef_y02,coef_z02,coef_r02)
     
    ! Scatter variables
     alpha=0.69d0
     x01   = coef_x01*boxlen! position of the centre of the mask
     y01   = coef_y01*boxlen!
     z01   = coef_z01*boxlen
     r01   = coef_r01*alpha!

     x02   = coef_x02*boxlen! position of the centre of the mask
     y02   = coef_y02*boxlen!
     z02   = coef_z02*boxlen
     r02   = coef_r02*alpha


! mettre 3D!
     do ivar=1,nvar
         do i=1,ngrid
            xcen=xx(i,1)-x01 ! position cell/centre          !!!
            ycen=xx(i,2)-y01
            rc=sqrt(xcen**2+ycen**2)! 2d distance /centre          !!!
            if (rc .lt. 2*r01) then
               uold(ind_cell(i),ivar)=u_mask(i,ivar)
            endif
             
            xcen=xx(i,1)-x02 ! position cell/centre          !!!
            ycen=xx(i,2)-y02
            rc=sqrt(xcen**2+ycen**2)! 2d distance /centre          !!!
            if (rc .lt. 2*r02) then
               uold(ind_cell(i),ivar)=u_mask(i,ivar)
            endif
         end do
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask


subroutine condinit2(x,u,dx,nn,coef_x0,coef_y0,coef_z0,coef_r0)
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
  real(dp):: R_star, beta,rho_amb,P_amb,P0,rho0
  integer::  ivar
  real(dp):: eta,nu, alpha, v_w,v_inf
  real(dp):: coef_x0,coef_y0,coef_z0,coef_r0
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  IF (nn .eq. 0) RETURN


! faire le cas 1D!

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

!run parameters
  eta = 0.065
   nu = 0.85
   alpha = 0.69


  x0   = coef_x0*boxlen
  y0   = coef_y0*boxlen
  r0   = coef_r0*alpha !0.195
  z0   = coef_z0*boxlen
  R_star = 0.0539 !(if alpha=0.69)
  beta = 1.0
  rho_amb = 0.1*eta*alpha**(-2)*nu**(-1)
   P_amb = 45.d0*nu*eta*alpha**(-2)
  P0 = 361*nu*eta*alpha**(-2)
  rho0 = 540*eta*alpha**(-2)*nu**(-1)
  v_w=nu
  v_inf = 1.36*v_w


  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 ! position cell/centre
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)! 2d distance /centre
#if NDIM == 3
     zz=x(i,3)-z0
     rc3=sqrt(xx**2+yy**2+zz*2)
     IF (rc3 .lt. 2*r0) THEN
        q(i,id) = rho0*(r0/rc)**2
        q(i,iu) = v_inf*(xx/rc3)*(1-R_star/rc3)**beta
        q(i,iv) = v_inf*(yy/rc3)*(1-R_star/rc3)**beta
        q(i,iw) = v_inf*(zz/rc3)*(1-R_star/rc3)**beta
        q(i,ip) = P0*(r0/rc)**(10/3)

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
        !write(*,*)'on est a 2D'
        IF (rc .lt. 2*r0) THEN
           q(i,id) = integration2D(xx,yy)! mettre la precision du decoupage en argument!!!!
           !q(i,id) = min(rho0*(r0/rc)**2,4.d0*rho0)
           q(i,iu) = v_inf*(xx/rc)*max((1-R_star/rc)**beta,0.d0)
           q(i,iv) = v_inf*(yy/rc)*max((1-R_star/rc)**beta,0.d0)
           q(i,ip) = min(P0*(r0/rc)**(10/3),2.d0**(2.d0*gamma)*P0)

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



subroutine integration2D(xx,yy)


! just for 1 cell for the moment

integer ::n=10 ! number of slices in a cell, in a direction

integer :: i,j,k
!real(kind=8)::x,y,z
real(kind=8)::density

real(kind=8)::integral=0.0d0
real(kind=8)::dx=0.1d0,dy=0.1d0! siez of a slice

real(kind=8)::xx,yy   ! position of the center of the cell 

do i=1,n
   do j=1,n
! computing the integral in the middle of the slice
         integral =integral +dx*dy*density2D(xx+(i-5-0.5d0)*dx,yy+(j-5-0.5d0)*dy) 
    end do
end do  

write(*,*),integral

end subroutine integration2D



function density2D(x,y) ! function which returns the density for a given position
real(kind=8)::x,y
real(kind=8)::rho_0 != 540.0d0 
real(kind=8)::r0    = 0.1d0
real(kind=8)::r ! distance from center of the mask

!rho_0=540.0d0
rho_0=8.0d0
!density=a*x*y*z

r=sqrt(x**2+y**2)
density=rho_0*(r0/r)**2
!density=rho_0*r

return
end function density2D
