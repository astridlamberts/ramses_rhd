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
  real(dp):: x01,y01,z01,x02,y02,z02,rc,rs,xx,yy,zz,r0_1,r0_2
  real(dp):: rho_amb,P_amb
  real(dp):: R_star1,v_inf1,Mach_1
  real(dp):: R_star2,v_inf2,Mach_2
  real(dp)::rho_1,P_1,v_1,rho_2,P_2,v_2  ! variables for r =a
  integer :: ivar
  real(dp)::beta
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: integration_rho,integration_p,integration_vx,integration_vy,integration_vz ! function which computes the mean value of a variable in the mask
  real(dp):: fract
  real(dp):: M_dot1,M_dot2 ! mass loss rate
  real(dp):: delta_x
  real(dp)::a ! physical separation of the system, expressed in a AU
  real(dp)::a1=1.0d0 ! separation of the system in simulation units
  real(dp)::pi=3.141459265
  real(dp)::epsilon_star1,epsilon_star2
  real(dp)::R_p,R_rho
  real(dp)::period
  real(dp)::eta
  real(dp)::integration_rho_pulsar,integration_P_pulsar


  IF (nn .eq. 0) RETURN

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  delta_x=boxlen/2.0d0**(levelmin)
  
  beta = 1.0d0  ! a ameliorer

!run parameters
   
call read_namelist(a,epsilon_star1,epsilon_star2,r0_1,r0_2,M_dot1,M_dot2,v_inf1,v_inf2&
     ,Mach_1,Mach_2,R_rho,R_P,period)

!positions of the stars
  x01   = 0.5d0*boxlen!-(a/2.0d0)
  y01   = 0.5d0*boxlen
  z01   = 0.5d0*boxlen
  
  x02   = 0.5d0*boxlen+a1 !+(a/2.0d0)
  y02   = 0.5d0*boxlen
  z02   = 0.5d0*boxlen


! compute  variables at 'a'
  v_1 = v_inf1*(1.0d0-epsilon_star1)**beta 
  rho_1 = M_dot1/(4.0d0*pi*a**2*v_1)

  v_2 = v_inf2!*(1.0d0-epsilon_star2)**beta 
  rho_2 = M_dot2/(4.0d0*pi*a**2*v_2)
!
 ! write(*,*),r0_2,'r02 ini'

#if NDIM==2
 rho_1 = M_dot1/(2.0d0*pi*a*v_1)
 rho_2 = M_dot2/(2.0d0*pi*a*v_2)
#endif

  P_1 = v_1**2* rho_1/(Mach_1**2*gamma)
  P_2 = v_2**2* rho_2/(Mach_2**2*gamma)

! ambient medium

  P_amb= P_1/R_P
  rho_amb=rho_1/R_rho

!computing eta_w=v1*Mdot1/v2*Mdot2. The position of the stagnation point should only depend on this parametre
!  eta=(M_dot1*v_inf1)/(M_dot2*v_inf2)
  


 ! There is an initial domain filled with a densité rho_amb and sound speed cs_amb. Afterwards, we add a mask (for r <r0) with a given density, pressure and speed. 

q(:,id) = rho_amb
q(:,iu) = 0.0d0
q(:,iv) = 0.0d0
#if NDIM == 3
  q(:,iw) = 0.0d0
#endif
q(:,ip) = P_amb



! First star!!!!!!!!!!!
zz = 0.0d0  !dummy initialisation

  DO i=1,nn ! all the cells
     xx=x(i,1)-x01 ! position cell/center
     yy=x(i,2)-y01
     rc=sqrt(xx**2+yy**2)! 2d distance /center
!has to be checked
#if NDIM == 3 
     zz=x(i,3)-z01
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0_1) THEN
        q(i,iu) =  fract(xx,yy,zz,r0_1,delta_x)*integration_vx(xx,yy,zz,v_inf1,epsilon_star1&
            ,beta ,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,iu)
        q(i,iv) =  fract(xx,yy,zz,r0_1,delta_x)*integration_vy(xx,yy,zz,v_inf1,epsilon_star1&
             ,beta,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,iv)
        q(i,iw) =  fract(xx,yy,zz,r0_1,delta_x)*integration_vz(xx,yy,zz,v_inf1,epsilon_star1&
             ,beta,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,iw)
        q(i,id) =  fract(xx,yy,zz,r0_1,delta_x)*integration_rho(xx,yy,zz,M_dot1,r0_1,delta_x&
             ,a,v_inf1,epsilon_star1,beta)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,id)
        q(i,ip) =  fract(xx,yy,zz,r0_1,delta_x)*integration_P(xx,yy,zz,P_1,rho_1,M_dot1,a,r0_1&
            ,v_inf1,epsilon_star1,beta,gamma,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,ip)


      ENDIF
#endif

#if NDIM==2   
!      write(*,*),rc,r0_1,'rc,r_01'
        IF (rc .lt. 1.2d0*r0_1) THEN
           !write(*,*),rc,r0_1,'rc,r0_1' 
       q(i,iu) =  fract(xx,yy,zz,r0_1,delta_x)*integration_vx(xx,yy,zz,v_inf1,epsilon_star1&
             ,beta,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,iu)
        q(i,iv) =  fract(xx,yy,zz,r0_1,delta_x)*integration_vy(xx,yy,zz,v_inf1,epsilon_star1&
             ,beta,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,iv)
        q(i,id) =  fract(xx,yy,zz,r0_1,delta_x)*integration_rho(xx,yy,zz,M_dot1,r0_1,delta_x,a&
             ,v_inf1,epsilon_star1,beta)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,id)
        q(i,ip) =  fract(xx,yy,zz,r0_1,delta_x)*integration_P(xx,yy,zz,P_1,rho_1,M_dot1,a,r0_1&
           ,v_inf1 ,epsilon_star1,beta,gamma,delta_x)+(1.d0-fract(xx,yy,zz,r0_1,delta_x))*q(i,ip)

	ENDIF
#endif
   ENDDO     


!pulsar
zz = 0.0d0  

  DO i=1,nn 
     xx=x(i,1)-x02 
     yy=x(i,2)-y02
     rc=sqrt(xx**2+yy**2)
!has to be checked
#if NDIM == 3
     zz=x(i,3)-z02
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0_2) THEN
        q(i,iu) = fract(xx,yy,zz,r0_2,delta_x)*v_inf2*(xx/rc)+(1.0d0-fract(xx,yy,zz,r0_2,delta_x)&
             *q(i,iu))
        q(i,iv) = fract(xx,yy,zz,r0_2,delta_x)*v_inf2*(yy/rc)+(1.0d0-fract(xx,yy,zz,r0_2,delta_x)&
             *q(i,iv))
        q(i,iw) = fract(xx,yy,zz,r0_2,delta_x)*v_inf2*(zz/rc)+(1.0d0-fract(xx,yy,zz,r0_2,delta_x)&
             *q(i,iw))
        q(i,id) =  fract(xx,yy,zz,r0_2,delta_x)*integration_rho_pulsar(xx,yy,zz,M_dot2,r0_2&
             ,delta_x,a,v_inf2)+(1.d0-fract(xx,yy,zz,r0_2,delta_x))*q(i,id)
        q(i,ip) =  fract(xx,yy,zz,r0_2,delta_x)*integration_P_pulsar(xx,yy,zz,P_2,rho_2,M_dot2,a&
             ,r0_2,v_inf2,gamma,delta_x)+(1.d0-fract(xx,yy,zz,r0_2,delta_x))*q(i,ip)
 !       write(*,*),q(i,ip),q(i,id),'d,p'
     
      ENDIF
#endif

#if NDIM==2   

     ! write(*,*),rc,r0_2,'rc,r0_2'
        IF (rc .lt. 1.2d0*r0_2) THEN
   !        write(*,*),'pulsar'
        q(i,iu) = fract(xx,yy,zz,r0_2,delta_x)*v_inf2*(xx/rc)+(1.0d0-fract(xx,yy,zz,r0_2,delta_x)&
             *q(i,iu))
        q(i,iv) = fract(xx,yy,zz,r0_2,delta_x)*v_inf2*(yy/rc)+(1.0d0-fract(xx,yy,zz,r0_2,delta_x)&
             *q(i,iv))
        q(i,id) =  fract(xx,yy,zz,r0_2,delta_x)*integration_rho_pulsar(xx,yy,zz,M_dot2,r0_2&
             ,delta_x,a,v_inf2)+(1.d0-fract(xx,yy,zz,r0_2,delta_x))*q(i,id)
        q(i,ip) =  fract(xx,yy,zz,r0_2,delta_x)*integration_P_pulsar(xx,yy,zz,P_2,rho_2,M_dot2,a&
             ,r0_2,v_inf2,gamma,delta_x)+(1.d0-fract(xx,yy,zz,r0_2,delta_x))*q(i,ip)
!    write(*,*),q(i,ip),q(i,id),'d,p'
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
  use hydro_parameters
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
  real(dp)::x01,y01,z01,rc1,xcen1,ycen1,zcen1,r0_1,r0_2,x02,y02,z02,xcen2,ycen2,zcen2,rc2
  real(dp)::fract
  real(dp)::v_inf1,v_inf2
  logical::error,ok_file1,ok_file2,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar
  real(dp)::a ! physical separation of the system, in AU
  real(dp)::a1=1.0d0  ! separation of the system in the simulation, it's the leght scale
  real(dp)::M_dot1,M_dot2,Mach_1,Mach_2,v_1,rho_1,P_1,v_2,rho_2,P_2
  real(dp):: epsilon_star1,epsilon_star2
  real(dp)::beta,R_rho,R_P
  real(dp)::boxlen_0 !size of small box
  real(dp)::delta_x!
  real(dp)::force,velocity_fluid
  real(dp)::period 

   beta = 1.0d0 ! a ameliorer!  
   zcen2 = 0.0d0
   zcen1 = 0.0d0 


call read_namelist(a,epsilon_star1,epsilon_star2,r0_1,r0_2,M_dot1,M_dot2,v_inf1,v_inf2&
     ,Mach_1,Mach_2,R_rho,R_P,period)

!write(*,*),a,epsilon_star1,epsilon_star2,r0_1,r0_2,M_dot1,M_dot2,v_inf1,v_inf2,Mach_1,Mach_2,R_rho,R_P,period
!values for r=a
!star
  v_1 = v_inf1*(1.0d0-epsilon_star1)**beta 
  rho_1 = M_dot1/(4.0d0*3.1415*a**2*v_1)
!pulsar
  v_2 = v_inf2!*(1.0d0-epsilon_star2)**beta 
  rho_2 = M_dot2/(4.0d0*3.1415*a**2*v_2)

 
#if NDIM==2
 rho_1 = M_dot1/(2.0d0*3.1415*a*v_1)
 rho_2 = M_dot2/(2.0d0*3.1415*a*v_2)
#endif

  P_1 = v_1**2* rho_1/(Mach_1**2*gamma)
  P_1 = v_2**2* rho_2/(Mach_2**2*gamma)

 !positions of the stars

  x01   = 0.5d0*boxlen!-(a/2.0d0)
  y01   = 0.5d0*boxlen
  z01   = 0.5d0*boxlen
 ! write(*,*),t,period,cos(2.0D0*3.1415/period*t),'t,T,cos'
  
  x02   = 0.5d0*boxlen+a1*cos(2.0D0*3.1415/period *t) !+(a/2.0d0)
  y02   = 0.5d0*boxlen+a1*sin(2.0d0*3.1415/period *t)
  z02   = 0.5d0*boxlen


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
   ! Call reset of othe mask routine, changes uold in the mask
     call condinit2(xx,u_mask,dx_loc,ngrid,x01,y01,z01,r0_1,a,M_dot1,v_inf1,Mach_1,epsilon_star1&
          ,beta,delta_x)
    ! write(*,*),delta_x,'reset_mask'
     call condinit2_pulsar(xx,u_mask,dx_loc,ngrid,x02,y02,z02,r0_2,a,M_dot2,v_inf2,Mach_2,delta_x)
     
         do i=1,ngrid
            xcen1=xx(i,1)-x01 
            ycen1=xx(i,2)-y01
            xcen2=xx(i,1)-x02 
            ycen2=xx(i,2)-y02
            rc1=sqrt(xcen1**2+ycen1**2)
            rc2=sqrt(xcen2**2+ycen2**2)
! a regrouper avec le reste du 3D?
#if NDIM ==3 
            zcen1=xx(i,3)-z01
            rc1=sqrt(xcen1**2+ycen1**2+zcen1**2) 
            zcen2=xx(i,3)-z02
            rc2=sqrt(xcen2**2+ycen2**2+zcen2**2) 

#endif
         !   write(*,*),xcen1,xx(i,1),x01,'xcen,xx,x01'

! in the star
#if NDIM==2

            if (rc1 .lt. 1.2d0*r0_1) then
 
               write(*,*),rc1,r0_1,'rc1,r01,etoile'
              uold(ind_cell(i),1)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,1)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*uold(ind_cell(i),1)
              !write(*,*),rc1,r0_1,xcen1,ycen1,zcen1,'rc1,r0_1,xcen,ycen,zcen'


!#if NDIM==2

              velocity_fluid=sqrt(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2)&
                   /uold(ind_cell(i),1)
              

              ! the force works so there is an extra energy density related to the speed of the flow
              uold(ind_cell(i),4)= fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,4)&
                   +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),4)&
                   +force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1)*velocity_fluid*dtnew(ilevel))

              ! the wind is accelerated by the force
              uold(ind_cell(i),2)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,2)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),2)&
                    +xcen1/rc1*dtnew(ilevel)*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1))

             uold(ind_cell(i),3)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,3)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),3)&
                    +ycen1/rc1*dtnew(ilevel)*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1))          

 !in the pulsar
           else if (rc2 .lt. 1.2d0*r0_2)then 
              write(*,*),rc1,rc2,r0_1,r0_2,'rc1,rc2,r0,pulsar'
              do ivar=1,nvar
                 uold(ind_cell(i),ivar)=fract(xcen2,ycen2,zcen2,r0_2,delta_x)*u_mask(i,ivar)&
                    +(1.d0-fract(xcen2,ycen2,zcen2,r0_2,delta_x))*uold(ind_cell(i),ivar)
              enddo   

! the ambient medium
           elseif(( rc1 .gt. 1.2d0*r0_1) .and. (rc2 .gt. 1.2d0*r0_2))then!
write(*,*),rc1,rc2,'autre'

                  velocity_fluid=sqrt(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2)&
                   /uold(ind_cell(i),1)
              uold(ind_cell(i),4)=uold(ind_cell(i),4)&
              +force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1)*velocity_fluid*dtnew(ilevel)

              uold(ind_cell(i),2)=uold(ind_cell(i),2)&
              +xcen1/rc1*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1 )*dtnew(ilevel)

             uold(ind_cell(i),3)=uold(ind_cell(i),3)&
              +ycen1/rc1*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1 )*dtnew(ilevel)
  

             endif !!!
#endif
! has to be checked
#if NDIM==3 

             if (rc1 .lt. 1.2d0*r0_1) then
 
              uold(ind_cell(i),1)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,1)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*uold(ind_cell(i),1)

              velocity_fluid=sqrt(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2)&
                   /uold(ind_cell(i),1)
              

              ! the force works so there is an extra energy density related to the speed of the flow
              uold(ind_cell(i),5)= fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,5)&
                   +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),5)&
                   +force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1)*velocity_fluid*dtnew(ilevel))

              ! the wind is accelerated by the force
              uold(ind_cell(i),2)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,2)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),2)&
                    +xcen1/rc1*dtnew(ilevel)*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1))

             uold(ind_cell(i),3)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,3)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),3)&
                    +ycen1/rc1*dtnew(ilevel)*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1))          

             uold(ind_cell(i),4)=fract(xcen1,ycen1,zcen1,r0_1,delta_x)*u_mask(i,4)&
                    +(1.d0-fract(xcen1,ycen1,zcen1,r0_1,delta_x))*(uold(ind_cell(i),4)&
                    +zcen1/rc1*dtnew(ilevel)*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1))          


! in the pulsar
           else if (rc2 .lt. 1.2d0*r0_2)then 
              do ivar=1,nvar
                 uold(ind_cell(i),ivar)=fract(xcen2,ycen2,zcen2,r0_2,delta_x)*u_mask(i,ivar)&
                    +(1.d0-fract(xcen2,ycen2,zcen2,r0_2,delta_x))*uold(ind_cell(i),ivar)
              enddo   

! the ambient medium
           elseif(( rc1 .gt. 1.2d0*r0_1) .and. (rc2 .gt. 1.2d0*r0_2))then!


             velocity_fluid=sqrt(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2)&
                   /uold(ind_cell(i),1)
              uold(ind_cell(i),5)=uold(ind_cell(i),5)&
              +force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1)*velocity_fluid*dtnew(ilevel)

              uold(ind_cell(i),2)=uold(ind_cell(i),2)&
              +xcen1/rc1*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1 )*dtnew(ilevel)

             uold(ind_cell(i),3)=uold(ind_cell(i),3)&
              +ycen1/rc1*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1 )*dtnew(ilevel)
  
             uold(ind_cell(i),4)=uold(ind_cell(i),4)&
              +zcen1/rc1*force(v_inf1,a,epsilon_star1,rho_1,v_1,rc1 )*dtnew(ilevel)
             
             endif  !!

#endif


         !  endif 
     end do
   end do
        ! End loop over cells
 end do
     ! End loop over grids

 end subroutine reset_mask


subroutine condinit2(x,u,dx,nn,x0,y0,z0,r0,a,M_dot,v_inf,Mach,epsilon_star,beta,delta_x)
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
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0
  real(dp):: epsilon_star, beta
  integer::  ivar
  real(dp):: v_inf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: integration_rho,integration_P
  real(dp):: integration_vx,integration_vy,integration_vz
  real(dp):: delta_x
  real(dp):: boxlen_0
  real(dp):: M_dot,Mach
  real(dp)::v_a,rho_a,P_a,a
  IF (nn .eq. 0) RETURN

  v_a = v_inf*(1.0d0-epsilon_star)**beta 
  rho_a = M_dot/(4*3.1415d0*a**2*v_a)
  
#if NDIM==2
  rho_a = M_dot/(2*3.1415d0*a*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach**2*gamma)
! write(*,*),P_a,v_a,rho_a,v_inf,epsilon_star,M_dot,'P,v,d,vinf,eps,Mdot etoile'
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

        q(i,iu) =  integration_vx(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
        q(i,iv) =  integration_vy(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
        q(i,iw) =  integration_vz(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
        q(i,id) =  integration_rho(xx,yy,zz,M_dot,r0,delta_x,a,v_inf,epsilon_star,beta)
  q(i,ip) = integration_P(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,epsilon_star,beta,gamma,delta_x)
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
    
        q(i,iu) = integration_vx(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
        q(i,iv) = integration_vy(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
        q(i,id) = integration_rho(xx,yy,zz,M_dot,r0,delta_x,a,v_inf,epsilon_star,beta)
        q(i,ip) = integration_P(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,epsilon_star,beta,gamma,delta_x)

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



subroutine condinit2_pulsar(x,u,dx,nn,x0,y0,z0,r0,a,M_dot,v_inf,Mach,delta_x)
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
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0
  real(dp):: epsilon_star, beta
  integer::  ivar
  real(dp):: v_inf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: integration_rho_pulsar,integration_P_pulsar
  real(dp):: integration_vx,integration_vy,integration_vz
  real(dp):: delta_x
  real(dp):: boxlen_0
  real(dp):: M_dot,Mach
  real(dp)::v_a,rho_a,P_a,a
  IF (nn .eq. 0) RETURN

! write(*,*),delta_x,'pulsar'
  v_a = v_inf
  rho_a = M_dot/(4*3.1415d0*a**2*v_a)
  
#if NDIM==2
  rho_a = M_dot/(2*3.1415d0*a*v_a)
#endif

  P_a = v_a**2* rho_a/(Mach**2*gamma)
!  write(*,*),P_a,v_a,rho_a,v_inf,epsilon_star,M_dot,'P,v,d,vinf,eps,Mdot'

  id=1; iu=2; iv=3; iw=4; ip=ndim+2

  zz=0.0d0 

!  write(*,*),r0,'rO pulsar'
  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,iu) = v_inf*(xx/rc)
        q(i,iv) = v_inf*(yy/rc)
        q(i,iw) = v_inf*(zz/rc)
        q(i,id) =  integration_rho_pulsar(xx,yy,zz,M_dot,r0,delta_x,a,v_inf)
       q(i,ip) = integration_P_pulsar(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,gamma,delta_x)
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


!write(*,*),r0,'test'

#if NDIM==2   
        IF (rc .lt. 1.2d0*r0) THEN
!           write(*,*),rc,r0,xx,yy,'r,x,y,r0 condinit'  
        q(i,id) =  integration_rho_pulsar(xx,yy,zz,M_dot,r0,delta_x,a,v_inf)
        q(i,ip) = integration_P_pulsar(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,gamma,delta_x)
!        write(*,*),xx,yy,rc,q(i,iu),q(i,iv),q(i,id),q(i,ip),'x,y,r,vx,vy,d,o'
      !  write(*,*),q(i,id)

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


end subroutine condinit2_pulsar



function integration_rho(xx,yy,zz,M_dot,r0,delta_x,a,v_inf,epsilon_star,beta)
  use amr_parameters
implicit none


real(kind=8)::integration_rho
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::density
real(kind=8)::integral
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::M_dot,a,v_inf,epsilon_star,beta
real(kind=8)::r0 
real(kind=8)::delta_x


volume= delta_x**ndim

integral=0.0d0

#if NDIM==2
   do i=1,n
       do j=1,n
! computing the integral in the middle of the slice, be careful with size of the box!
         integral = integral+(delta_x/n)**ndim*density((xx+(i-5.0d0-0.5d0)*delta_x/n)&
              ,(yy+(j-5.0d0-0.5d0)*delta_x/n),zz,M_dot,r0,a,v_inf,epsilon_star,beta)
       end do
   end do  
#endif

#if NDIM==3
   do i=1,n
       do j=1,n
          do k=1,n
             integral =integral +(delta_x/n)**ndim*density((xx+(i-5.0d0-0.5d0)*delta_x/n)&
                  ,(yy+(j-5.0d0-0.5d0)*delta_x/n),(zz+(k-5.0d0-0.5d0)*delta_x/n)&
                  ,M_dot,r0,a,v_inf,epsilon_star,beta)
          end do
       end do
   end do
#endif

integration_rho=integral/volume

return
end function

function density(x,y,z,M_dot,r0,a,v_inf,epsilon_star,beta) 
  use amr_parameters
implicit none

real(kind=8) :: density,velocity
real(kind=8)::x,y,z
real(kind=8)::r0 
real(kind=8)::r ! distance from center of the mask
real(kind=8)::M_dot,a,v_inf,epsilon_star,beta

r=sqrt(x**2+y**2+z**2)
 
#if NDIM==2
    density=M_dot/(2.0d0*3.1415d0*r*a*velocity(r,v_inf,epsilon_star,beta))
#endif
 
#if NDIM==3
    density=M_dot/(4.0d0*3.1415d0*(r*a)**2*velocity(r,v_inf,epsilon_star,beta))
#endif

return
end function


function integration_P(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,epsilon_star,beta,gamma,delta_x)
  use amr_parameters
implicit none


real(kind=8)::integration_P
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::density
real(kind=8)::integral
real(kind=8)::xx,yy,zz   
real(kind=8)::volume
real(kind=8)::P_a,rho_a! values at 'a'
real(kind=8)::gamma,M_dot,v_inf,epsilon_star,beta,a
real(kind=8)::delta_x,r0

volume= delta_x**ndim

integral=0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         integral =integral +delta_x/n*delta_x/n*P_a*(density((xx+(i-5.0d0-0.5d0)*delta_x/n)&
          ,(yy+(j-5.0d0-0.5d0)*delta_x/n),zz,M_dot,r0,a,v_inf,epsilon_star,beta)/rho_a)**gamma 
      end do
   end do
#endif

#if NDIM==3
    do k=1,n
       do i=1,n
          do j=1,n
             integral =integral +(delta_x/n)**ndim*P_a*(density((xx+(i-5.0d0-0.5d0)*delta_x/n)&
          ,(yy+(j-5.0d0-0.5d0)*delta_x/n),(zz+(k-5.0d0-0.5d0)*delta_x/n),M_dot,r0,a,v_inf&
          ,epsilon_star,beta)/rho_a)**gamma 
          end do
       end do
    end do
#endif

integration_P=integral/volume

return
end function

function integration_vx(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
  use amr_parameters
implicit none


real(kind=8)::integration_vx
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::velocity
real(kind=8)::integral
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::v_inf
real(kind=8)::epsilon_star
real(kind=8)::delta_x
real(kind=8)::beta
real(kind=8)::r



volume= delta_x**ndim

integral=0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
         r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2)
         integral =integral +delta_x/n*delta_x/n*(xx/r)*velocity(r,v_inf,epsilon_star,beta)
    end do
   end do
#endif

#if NDIM==3 
   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2&
                  +(zz+(k-5.0d0-0.5d0)*delta_x/n)**2)
            integral =integral +(delta_x/n)**ndim*(xx/r)*velocity(r,v_inf,epsilon_star,beta)
         end do
      end do
   enddo
#endif

integration_vx=integral/volume

return
end function 



function integration_vy(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
  use amr_parameters
implicit none


real(kind=8)::integration_vy
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::velocity
real(kind=8)::integral
real(kind=8)::xx,yy,zz  
real(kind=8)::volume
real(kind=8)::delta_x ! size of the cell
real(kind=8)::v_inf
real(kind=8)::epsilon_star
real(kind=8)::beta
real(kind=8)::r

volume= delta_x**ndim

integral=0.0d0

#if NDIM==2
   do i=1,n
      do j=1,n
             r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2)
         integral =integral +(delta_x/n)**ndim*(yy/r)*velocity(r,v_inf,epsilon_star,beta)
      end do
    end do
#endif

#if NDIM==3
   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2&
                  +(zz+(k-5.0d0-0.5d0)*delta_x/n)**2)
         integral =integral +(delta_x/n)**ndim*(yy/r)*velocity(r,v_inf,epsilon_star,beta)
         end do
      end do
   end do 
#endif

integration_vy=integral/volume

return
end function integration_vy


function integration_rho_pulsar(xx,yy,zz,M_dot,r0,delta_x,a,v_inf)
  use amr_parameters
implicit none


real(kind=8)::integration_rho_pulsar
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::density
real(kind=8)::integral
real(kind=8)::xx,yy,zz   ! position of the center of the cell 
real(kind=8)::volume
real(kind=8)::M_dot,a,v_inf,epsilon_star,beta
real(kind=8)::r0 
real(kind=8)::delta_x
real(kind=8)::r

volume= delta_x**ndim

integral=0.0d0
!write(*,*),xx,yy,sqrt(xx**2+yy**2),'x,y,r,density'

#if NDIM==2
   do i=1,n
       do j=1,n
          r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2)
          integral= integral+(delta_x/n)**ndim*(1.0d0/r)
!          write(*,*),r,xx,yy,'r,x,y'
       end do
   end do  
integral=integral*M_dot/(2*3.1415*a*v_inf)
#endif
!write(*,*),xx,yy,r,integral,'x,y,r,integral'

#if NDIM==3
   do i=1,n
       do j=1,n
          do k=1,n
             r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2&
                  +(zz+(k-5.0d0-0.5d0)*delta_x/n)**2)
             integral =integral +(delta_x/n)**ndim*1.0d0/(r**2)
          end do
       end do
   end do
integral=integral*M_dot/(4*3.1415*v_inf*a**2)
#endif


!write(*,*),xx,yy,zz,M_dot,r0,delta_x,a,v_inf,sqrt(xx**2+yy**2),M_dot/(volume*2*3.14*a*v_inf*sqrt(xx**2+yy**2)),'x,y,z,Mdot,d0,deltax,a,vinf,r,d'

integration_rho_pulsar=integral/volume
!write(*,*),integration_rho_pulsar,integral,volume
!write(*,*),M_dot,v_inf
return
end function













function integration_P_pulsar(xx,yy,zz,P_a,rho_a,M_dot,a,r0,v_inf,gamma,delta_x)
  use amr_parameters
implicit none
real(kind=8)::integration_P_pulsar
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::density
real(kind=8)::integral
real(kind=8)::xx,yy,zz   
real(kind=8)::volume
real(kind=8)::P_a,rho_a! values at 'a'
real(kind=8)::gamma,M_dot,v_inf,a
real(kind=8)::delta_x,r0
real(kind=8)::r

volume= delta_x**ndim

integral=0.0d0
!write(*,*),P_a,'Pa pression'

#if NDIM==2
   do i=1,n
      do j=1,n
         r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2)
!         integral =integral +(delta_x/n)**ndim*(1.0d0/r)**gamma
         integral =integral+(a/r)**gamma
     end do
   end do
!integral=integral*P_a*(M_dot/(2.0d0*3.1415d0*a*v_inf*rho_a))**(-gamma)   
integral=integral*P_a
#endif

#if NDIM==3
    do k=1,n
       do i=1,n
          do j=1,n
             r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2&
                  + (zz+(k-5.0d0-0.5d0)*delta_x/n)**2)
!             integral =integral +(delta_x/n)**ndim**(1.0d0/r**2)**(-gamma)
             integral =integral+(a/r)**gamma
         end do
       end do
    end do
integral=integral*P_a
!integral=integral*P_a*(M_dot/(4*3.1415*a**2*v_inf*rho_a))**gamma 
#endif


!write(*,*),integral/volume,r,P_a,gamma,'P,r,Pa,gam'

integration_P_pulsar=integral/volume


return
end function




function integration_vz(xx,yy,zz,v_inf,epsilon_star,beta,delta_x)
  use amr_parameters
implicit none


real(kind=8)::integration_vz
integer ::n=10 ! number of slices in a cell, in a direction
integer :: i,j,k
real(kind=8)::velocity
real(kind=8)::integral
real(kind=8)::xx,yy,zz  
real(kind=8)::volume
real(kind=8)::delta_x ! size of the cell
real(kind=8)::v_inf
real(kind=8)::epsilon_star
real(kind=8)::beta
real(kind=8)::r

volume= delta_x**ndim

integral=0.0d0


   do i=1,n
      do j=1,n
         do k=1,n
            r=sqrt((xx+(i-5.0d0-0.5d0)*delta_x/n)**2+(yy+(j-5.0d0-0.5d0)*delta_x/n)**2 &
                  +(zz+(k-5.0d0-0.5d0)*delta_x/n)**2)
            integral =integral +(delta_x/n)**ndim*(zz/r)*velocity(r,v_inf,epsilon_star,beta)
         end do
      end do
   end do 


integration_vz=integral/volume

return
end function integration_vz

function velocity(r,v_inf,epsilon_star,beta) 
implicit none
real(kind=8):: velocity
real(kind=8)::epsilon_star 
real(kind=8)::v_inf   
real(kind=8)::r 
real(kind=8)::beta  

velocity=max((v_inf*(1.0d0-epsilon_star/r)**beta),v_inf/10)

return
end function velocity


function force(v_inf,a,epsilon_star,rho_a,v_a,rc)

implicit none
real(kind=8)::force
real(kind=8)::v_inf
real(kind=8)::a
real(kind=8)::epsilon_star
real(kind=8)::rho_a
real(kind=8)::v_a
real(kind=8)::rc

force=(a*v_inf*v_a*rho_a*epsilon_star)/(rc**3)

return
end function force




function fract(x,y,z,r0,delta_x)
 use amr_parameters

!function which computes the fraction of a cell being occupied by the mask
implicit none
real(kind=8):: fract
real(kind=8):: delta_x
real(kind=8):: x,y,z,r
real(kind=8):: x_s,y_s,z_s! postion on the center of the slice
real(kind=8):: r0
integer :: i,j,k
integer :: n=10 
real :: compteur


compteur=0.d0
#if NDIM==2
   do i=1,n
      do j=1,n
         x_s=x+(i-5.0d0-0.5d0)*delta_x/n
         y_s=y+(j-5.0d0-0.5d0)*delta_x/n
         r=sqrt(x_s**2+y_s**2)
         if (r .lt. r0)then
            compteur=compteur+1
         endif
      end do
   end do
#endif

#if NDIM==3
   do k=1,n
      do i=1,n
         do j=1,n
            x_s=x+(i-5.0d0-0.5d0)*delta_x/n
            y_s=y+(j-5.0d0-0.5d0)*delta_x/n
            z_s=z+(k-5.0d0-0.5d0)*delta_x/n
            r=sqrt(x_s**2+y_s**2+z_s**2)
            if (r .lt. r0)then
               compteur=compteur+1
            endif
         end do
      end do
   end do   
#endif

fract=compteur/(n**ndim)

return
end function fract


!subroutine which reads the additionnal namelist

!subroutine which reads the additionnal namelist
subroutine read_namelist(a,epsilon_star1,epsilon_star2,r0_1,r0_2,M_dot1,M_dot2,v_inf1,v_inf2&
     ,Mach_1,Mach_2,R_rho,R_P,period)

implicit none
!mpi?????
 logical::nml_ok
 real(kind=8)::M_dot1,M_dot2,v_inf1,v_inf2,r0_1,r0_2,epsilon_star1,epsilon_star2,Mach_1,Mach_2
 real(kind=8)::R_rho,R_P,a,period

 a=1.0d0

 M_dot1 = 1.0d0
 Mach_1 = 30.0d0
 v_inf1 = 422.0d0
 R0_1   = 0.092d0
 epsilon_star1 = 0.046d0

 M_dot2 = 1.0d0
 Mach_2 = 30.0d0
 v_inf2 = 422.0d0
 R0_2   = 0.092d0
 epsilon_star2 = 0.046d0 

 R_rho = 1000.0d0
 R_P   = 1.0d0
 period = 0.0001

! namelist definitions

 namelist/first_star/M_dot1,Mach_1,v_inf1,r0_1,epsilon_star1

 namelist/second_star/M_dot2,Mach_2,v_inf2,r0_2,epsilon_star2

 namelist/ambient_medium/a,R_rho,R_P,period

 INQUIRE(file='param.nml',exist=nml_ok)
  if(.not. nml_ok)then
    ! if(myid==1)then
    !    write(*,*)'File '//TRIM(infile)//' does not exist'
    ! endif
     write(*,*),'il manque le fichier param.nml'
     call clean_stop
  end if

  open(unit=10,file='param.nml')
  rewind(10)
  read(10,NML=first_star)
  rewind(10)
  read(10,NML=second_star)
  rewind(10)
  read(10,NML=ambient_medium)
  close(10)

!conversion according to the value of a
  v_inf1=v_inf1/a
  v_inf2=v_inf2/a
  
  r0_1  = r0_1/a
  r0_2  = r0_2/a

  epsilon_star1 = epsilon_star1/a
  epsilon_star2 = epsilon_star2/a
 

!  write(*,*),M_dot1,Mach_1,v_inf1,r0_1,epsilon_star1,'Mdot1,Mach1,vinf1,r01,eps1'
!  write(*,*),M_dot2,Mach_2,v_inf2,r0_2,epsilon_star2,a,'Mdot2,Mach2,vinf2,r02,eps2'
end subroutine read_namelist
