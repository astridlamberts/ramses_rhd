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
  integer:: i,j
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: rho_amb,P_amb,v_a,P_a,rho_a1,rho_a2,h,lor,lor1,lor2,P_a1,P_a2
  integer :: ivar,ncell
  real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
  real(dp)::a1=1.0d0 ! separation of the system in simulation units
  integer::index_pas_scal1=1
  integer::index_pas_scal2=2
  ncell=ncoarse+twotondim*ngridmax


  IF (nn .eq. 0) RETURN
  

  ! initialize the medium and masks
 lor1=(1.-vinf1**2)**(-1./2.)
 lor2=(1.-vinf2**2)**(-1./2.)

  rho_a1 = Mdot1*a1**2/(4.0d0*3.1415*vinf1*lor1)
  rho_a2 = Mdot2*a1**2/(4.0d0*3.1415*vinf2*lor2)
  
#if NDIM==2
  rho_a1 = Mdot1*a1/(2.0d0*3.1415*vinf1*lor1)
  rho_a2 = Mdot2*a1/(2.0d0*3.1415*vinf2*lor2)
#endif
  
! We assume the mach number given in the log file is the classical mach number
  P_a1 = vinf1**2* rho_a1/(gamma*Mach1**2-vinf1**2*gamma/(gamma-1))
  P_a2 = vinf2**2* rho_a2/(gamma*Mach2**2-vinf2**2*gamma/(gamma-1))

!  P_a1 = vinf1**2*lor1* rho_a1/(gamma*Mach1**2-vinf1**2*lor1**2*gamma/(gamma-1))
!  P_a2 = vinf2**2*lor2* rho_a2/(gamma*Mach2**2-vinf2**2*lor2**2*gamma/(gamma-1))
  P_amb=P_a1/R_P
  rho_amb=rho_a1/R_rho
  
  q(:,1) = rho_amb
  q(:,2) = 0.0d0
  q(:,3) = 0.0d0
  q(:,4) = 0.0d0
  q(:,5) = P_amb

!there are no passive scalars in the initial medium
  if (nvar .gt. 5) then
     q(:,6) =0.0d0 
     q(:,7) =0.0d0 
  endif


!initial position of the bodies
  
  call init_rot

  
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
 
  do i=1,nn
     
     ! Convert primitive to conservative variables
     ! specific enthalpy
     h=1.0d0+gamma/(gamma-1.0d0)*q(i,5)/q(i,1)
     ! Lorentz factor
     lor=(1.-(q(i,2)**2+q(i,3)**2+q(i,4)**2))**(-1.d0/2.d0)
     
     ! proper density  -> density in lab frame
     u(i,1)= q(i,1)*lor
     ! proper 3-velocity -> momentum in lab frame
     u(i,2)=lor**2*q(i,1)*h*q(i,2)
     u(i,3)=lor**2*q(i,1)*h*q(i,3)
     u(i,4)=lor**2*q(i,1)*h*q(i,4)
     ! isotropic gas pressure -> total fluid energy
     u(i,5)=lor**2*q(i,1)*h-q(i,5)
     
     !passive scalars
     do ivar=6,nvar
        u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)*lor
     end do
!     write(*,*),q(i,1),q(i,2),q(i,3),q(i,4),q(i,5),'q',i
!     write(*,*),u(i,1),u(i,2),u(i,3),u(i,4),u(i,5),'u'
  end do



end subroutine condinit



subroutine init_star(q,x,nn,vinf,rstar,Mdot,Mach,x0,y0,z0,r0,beta,index_pas_scal)

  use amr_parameters
  use wind_parameters
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
  real(dp)::pi
  real(dp)::x0,y0,z0
  integer:: i,j,index_pas_scal
  real(dp)::a1=1.0d0 ! value of the separation in a AU
 

  pi=acos(-1.0d0)
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


#if NDIM == 3 
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
#endif
  
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,2) =  frac(xx,yy,zz,r0,dx)* int_vx(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,2)
        q(i,3) =  frac(xx,yy,zz,r0,dx)* int_vy(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,3)
        q(i,4) =  frac(xx,yy,zz,r0,dx)* int_vz(xx,yy,zz,vinf,rstar,beta,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,4)
        q(i,1) =  frac(xx,yy,zz,r0,dx)*int_rho(xx,yy,zz,Mdot,r0,dx,vinf,rstar,beta)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,2)
        q(i,5) =  frac(xx,yy,zz,r0,dx)*int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,5)
        if (nvar  .gt. 5) then
           q(i,5+index_pas_scal)=frac(xx,yy,zz,r0,dx)
        endif
     ENDIF
  ENDDO     

end subroutine init_star


subroutine init_puls(q,x,nn,vinf,r0,Mdot,Mach,x0,y0,z0,index_pas_scal) 
  use wind_parameters
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx 
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp):: rc,rs,xx,yy,zz
  real(dp):: vinf,r0,Mdot,Mach,beta,lor
  real(dp)::rho_a,P_A,rho_m,p_m ! variables for r =a, r=r0
  real(dp),dimension(1:nvector,1:nvar)::q   ! Primitive variables
  real(dp):: int_rhoP,int_pP
  real(dp):: frac
  real(dp)::pi=3.14159265
  real(dp)::x0,y0,z0
  integer:: i,j
  real(dp)::a1=1.0d0
  integer::index_pas_scal

  dx=boxlen/2.0d0**(levelmin)
 
  lor=(1.-vinf**2)**(-1./2.)

 !compute variables at 'a'
  rho_a = Mdot*a1**2/(4.0d0*pi*vinf*lor)

#if NDIM==2
 rho_a = Mdot*a1/(2.0d0*pi*vinf*lor)
#endif

  P_a = vinf**2* rho_a/(Mach**2*gamma-vinf**2*gamma/(gamma-1))
! P_a = vinf**2* lor*rho_a/(Mach**2*gamma-vinf**2*lor**2*gamma/(gamma-1))
! check if 3D velocity is ok!!!!


  zz = 0.0d0  

  DO i=1,nn 
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
#endif
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,2) = frac(xx,yy,zz,r0,dx)*vinf*(xx/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,2)
        q(i,3) = frac(xx,yy,zz,r0,dx)*vinf*(yy/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,3)
        q(i,4) = frac(xx,yy,zz,r0,dx)*vinf*(zz/rc)+(1.0d0-frac(xx,yy,zz,r0,dx))*q(i,4)
        q(i,1) = frac(xx,yy,zz,r0,dx)*int_rhoP(xx,yy,zz,Mdot,r0,dx,vinf)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,1)
        q(i,5) = frac(xx,yy,zz,r0,dx)*int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,dx)+(1.d0-frac(xx,yy,zz,r0,dx))*q(i,5)
        if (nvar  .gt. 5) then      
           q(i,5+index_pas_scal) = frac(xx,yy,zz,r0,dx)
        endif
     ENDIF
!write(*,*),rc,q(i,1),q(i,2),q(i,5),i
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
  real(dp)::v1,v2,rho1,rho2,lor1,lor2
  real(dp)::delta_x
  real(dp)::pex,pey,pez
!  real(dp)::ax02,ay02,r
  real(dp)::dt ! local timestep
  real(dp)::omega
  integer::index_pas_scal1=1
  integer::index_pas_scal2=2

  ncell=ncoarse+twotondim*ngridmax
  
  if (rotation .eqv. .true.)then
     if ((vx01 .eq. 0) .and. (vy01 .eq. 0)) then
        !if this is a restart with no initial rotation
        call init_rot
     endif
        
     if (ilevel .eq. nlevelmax) then
           ! The leapfrog is called only for finelevels otherwise there are too many updates  
        call leapfrog(x01,y01,z01,vx01,vy01,vz01,x02,y02,z02,vx02,vy02,vz02,Gstar,M1,M2,dtnew(ilevel))
     endif
     else
        x01=xcen+M2/(M1+M2)*a1*(1-exc)*(cos(node)*cos(peri+theta)-sin(node)*sin(peri+theta)*cos(inc))
        y01=ycen+M2/(M1+M2)*a1*(1-exc)*(sin(node)*cos(peri+theta)+cos(node)*sin(peri+theta)*cos(inc))
        z01=zcen+M2/(M1+M2)*a1*(1-exc)*sin(peri+theta)*sin(inc)
        
        x02=x01-a1*(1-exc)*(cos(node)*cos(peri+theta)-sin(node)*sin(peri+theta)*cos(inc))
        y02=y01-a1*(1-exc)*(sin(node)*cos(peri+theta)+cos(node)*sin(peri+theta)*cos(inc))
        z02=z01-a1*(1-exc)*sin(peri+theta)*sin(inc)
        
     endif
    
   
!inutile ici?
!values for r=a (mettre dans une fonction)
  !v1 = vinf1*(1.0d0-rstar1)**beta1 
     lor1=(1.-vinf1**2)**(-1./2.)
     rho1 = Mdot1*a1**2/(4.0d0*3.1415*vinf1*lor1)
  
  !v2  = vinf2*(1.0d0-rstar2)**beta2
     lor2=(1.-vinf2**2)**(-1./2.) 
     rho2 = Mdot2*a1**2/(4.0d0*3.1415*vinf2*lor2)
  
#if NDIM==2
     rho1 = Mdot1*a1/(2.0d0*3.1415*vinf1*lor1)
     rho2 = Mdot2*a1/(2.0d0*3.1415*vinf2*lor2)
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
               call reset_star(xc1,yc1,zc1,r01,dx,u_mask,i,ind_cell(i),ilevel,vinf1,rstar1,rho1,v1,rc1)!
            else if (type1=='puls')then
               call reset_puls(xc1,yc1,zc1,r01,dx,u_mask,i,ind_cell(i))
            endif

         else if (rc2 .lt. 1.2d0*r02)then 
            !2nd body: a star or a pulsar
            if (type2=='star')then
               call reset_star(xc2,yc2,zc2,r02,dx,u_mask,i,ind_cell(i),ilevel,vinf2,rstar2,rho2,v2,rc2)
            else if (type2=='puls')then    
               call reset_puls(xc2,yc2,zc2,r02,dx,u_mask,i,ind_cell(i))
            endif
         endif
      end do
   end do
   ! End loop over cells
end do
! End loop over grids

 end subroutine reset_mask




subroutine reset_star(xc,yc,zc,r0,delta_x,u_mask,i,ind_cell,ilevel,vinf,rstar,rho,v,rc)
                                       
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
      uold(ind_cell,2)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,2)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,2)
      uold(ind_cell,3)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,3)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,3)
      uold(ind_cell,4)=frac(xc,yc,zc,r0,delta_x)*u_mask(i,4)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,4)
      uold(ind_cell,5)= frac(xc,yc,zc,r0,delta_x)*u_mask(i,5)+(1.d0-frac(xc,yc,zc,r0,delta_x))*uold(ind_cell,5)

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


subroutine init2_star(x,u,dx,nn,x0,y0,z0,r0,Mdot,vinf,Mach,rstar,beta,delta_x,index_pas_scal)
    
!This routine computes the hydro variables in the mask region, in the case of a star
!it returns the vector u and requires all the parameters of the star as input

  use amr_parameters
  use hydro_parameters
  use hydro_commons  !!!sur?
  use wind_parameters
  implicit none

  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer:: i,j
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
  real(dp)::h,lor
  IF (nn .eq. 0) RETURN

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
#endif
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,1) =  int_rho(xx,yy,zz,Mdot,r0,delta_x,vinf,rstar,beta)
        q(i,2) =  int_vx(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,3) =  int_vy(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,4) =  int_vz(xx,yy,zz,vinf,rstar,beta,delta_x)
        q(i,5) =  int_P(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,rstar,beta,gamma,delta_x)
        if (nvar .gt. 5 ) then
           q(i,5+1)=0.0d0
           q(i,5+2)=0.0d0
           q(i,5+index_pas_scal)= 1.0d0
        endif
     
        ! Convert primitive to conservative variables
        ! specific enthalpy
        h=1.0d0+gamma/(gamma-1.0d0)*q(i,5)/q(i,1)
        ! Lorentz factor
        lor=(1.-(q(i,2)**2+q(i,3)**2+q(i,4)**2))**(-1.d0/2.d0)
        
        ! proper density  -> density in lab frame
        u(i,1)= q(i,1)*lor
        ! proper 3-velocity -> momentum in lab frame
        u(i,2)=lor**2*q(i,1)*h*q(i,2)
        u(i,3)=lor**2*q(i,1)*h*q(i,3)
        u(i,4)=lor**2*q(i,1)*h*q(i,4)
        ! isotropic gas pressure -> total fluid energy
        u(i,5)=lor**2*q(i,1)*h-q(i,5)
        
        !passive scalars
        do ivar=6,nvar
           u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)*lor
        end do

     ENDIF

   ENDDO 


end subroutine init2_star



subroutine init2_puls(x,u,dx,nn,x0,y0,z0,r0,Mdot,vinf,Mach,delta_x,index_pas_scal)
 
  use amr_parameters
  use hydro_parameters
  use hydro_commons 
  use wind_parameters 
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  
  integer:: i,j
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0
  real(dp):: rstar, beta
  integer::  ivar
  real(dp):: vinf
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp):: int_rhoP,int_Pp
  real(dp):: delta_x
  real(dp):: Mdot,Mach
  real(dp)::rho_a,p_a,p_m,rho_m,vinfm
  real(dp)::a1=1.0d0
  integer ::index_pas_scal
  real(dp)::h,lor
 
  IF (nn .eq. 0) RETURN
  lor=(1.-vinf**2)**(-1./2.)
  rho_a = Mdot*a1**2/(4.0d0*3.141592*vinf*lor)


#if NDIM==2
  rho_a = Mdot*a1/(2.0d0*3.141592*vinf*lor)
#endif

  P_a = vinf**2* rho_a/(Mach**2*gamma)
!  P_a = vinf**2*lor* rho_a/(gamma*Mach**2-vinf**2*lor**2*gamma/(gamma-1))
 ! correction to take into account the acceleration of the wind at short distance in 3D case

#if NDIM==3
!  write(*,*),vinf,'ini'
  rho_m=  Mdot/(4.0d0*3.141592*vinf*r0**2)
  P_m=P_a*(rho_a/rho_m)**(-gamma)
  vinf=sqrt(vinf**2-P_m/rho_m*gamma/(gamma-1.d0))
 ! write(*,*),vinf,P_m,rho_m
#endif

  zz=0.0d0 
  
  DO i=1,nn ! all the cells
     xx=x(i,1)-x0 
     yy=x(i,2)-y0
     rc=sqrt(xx**2+yy**2)
#if NDIM == 3
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2+zz**2)
#endif
     IF (rc .lt. 1.2d0*r0) THEN
        q(i,1) = int_rhoP(xx,yy,zz,Mdot,r0,delta_x,vinf)
        q(i,2) = vinf*(xx/rc)
        q(i,3) = vinf*(yy/rc)
        q(i,4) = vinf*(zz/rc)
        q(i,5) = int_pP(xx,yy,zz,P_a,rho_a,Mdot,r0,vinf,gamma,delta_x)
        if (nvar .gt. 5) then
           q(i,5+1)=0.0d0
           q(i,5+2)=0.0d0
           q(i,5+index_pas_scal)=1.0d0
        endif
        

        ! Convert primitive to conservative variables
        ! specific enthalpy
        h=1.0d0+gamma/(gamma-1.0d0)*q(i,5)/q(i,1)
        ! Lorentz factor
        lor=(1.-(q(i,2)**2+q(i,3)**2+q(i,4)**2))**(-1.d0/2.d0)
        
        ! proper density  -> density in lab frame
        u(i,1)= q(i,1)*lor
        ! proper 3-velocity -> momentum in lab frame
        u(i,2)=lor**2*q(i,1)*h*q(i,2)
        u(i,3)=lor**2*q(i,1)*h*q(i,3)
        u(i,4)=lor**2*q(i,1)*h*q(i,4)
        ! isotropic gas pressure -> total fluid energy
        u(i,5)=lor**2*q(i,1)*h-q(i,5)
        
        !passive scalars
        do ivar=6,nvar
           u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)*lor
        end do
     ENDIF
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
real(kind=8)::r,lor

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
lor=(1.-vinf**2)**(-1./2.)
int=int*Mdot/(2*3.1415*vinf*lor)


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
lor=(1.-vinf**2)**(-1./2.)
int=int*Mdot/(4*3.1415*vinf*lor)
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

implicit none

 logical::nml_ok

 namelist/first_body/type1,Mdot1,Mach1,vinf1,r01,rstar1,beta1,M1
 namelist/second_body/type2,Mdot2,Mach2,vinf2,r02,rstar2,beta2,M2
 namelist/ambient_medium/a,rotation,direction,xcen,ycen,zcen, R_P,R_rho,exc,Gstar,peri,inc,node,&
     &theta,x01,y01,z01

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
 !!!!!!!!a faire pour le barycentre!!!!!!!
  xcen = xcen*boxlen
  ycen = ycen*boxlen
  zcen = zcen*boxlen


  x01 = x01*boxlen
  y01 = y01*boxlen
  z01 = z01*boxlen

 Gstar=Gstar/a**3

end subroutine  read_wind_params


subroutine init_rot
use wind_parameters
use amr_commons
implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif


real (kind=8):: pi
real (kind=8):: a1=1.0d0
real (kind=8)::ax1,ay1,az1,ax2,ay2,az2
real (kind=8)::r,r1,r2,dthetadt
real (kind=8)::period
     
pi=acos(-1.0d0)

!inclination
inc=inc/180.d0*pi
!argument of pericenter
peri=peri/180.d0*pi
!longitude of ascending node
node=node/180.d0*pi
!initial position with respect to periastron
theta=theta/180.d0*pi

x01=xcen+M2/(M1+M2)*a1*(1-exc)*(cos(node)*cos(peri+theta)-sin(node)*sin(peri+theta)*cos(inc))
y01=ycen+M2/(M1+M2)*a1*(1-exc)*(sin(node)*cos(peri+theta)+cos(node)*sin(peri+theta)*cos(inc))
z01=zcen+M2/(M1+M2)*a1*(1-exc)*sin(peri+theta)*sin(inc)

x02=x01-a1*(1-exc)*(cos(node)*cos(peri+theta)-sin(node)*sin(peri+theta)*cos(inc))
y02=y01-a1*(1-exc)*(sin(node)*cos(peri+theta)+cos(node)*sin(peri+theta)*cos(inc))
z02=z01-a1*(1-exc)*sin(peri+theta)*sin(inc)

!if (rotation .eqv. .true.) then
   r=sqrt((x01-x02)**2+(y01-y02)**2+(z01-z02)**2)
   
   ax1=-Gstar*M2/r**2*(x01-x02)/r
   ay1=-Gstar*M2/r**2*(y01-y02)/r
   az1=-Gstar*M2/r**2*(z01-z02)/r
   
   ax2=-Gstar*M1/r**2*(x02-x01)/r
   ay2=-Gstar*M1/r**2*(y02-y01)/r
   az2=-Gstar*M1/r**2*(z02-z01)/r 
   
   dthetadt=sqrt(a1*(1-exc**2)*Gstar*(M1+M2)/(r**4))
   
   
   r1=sqrt((x01-xcen)**2.0d0+(y01-ycen)**2.0d0+(z01-zcen)**2.0d0)
   r2=sqrt((x02-xcen)**2.0d0+(y02-ycen)**2.0d0+(z02-zcen)**2.0d0)
   
   if (direction .eq.'d') then 
      ! direct rotation
      
      vx01=dthetadt*(-cos(node)*sin(peri+theta)-sin(node)*cos(peri+theta)*cos(inc))*r1
      vy01=dthetadt*(-sin(node)*sin(peri+theta)+cos(node)*cos(peri+theta)*cos(inc))*r1
      vz01=dthetadt*cos(peri+theta)*sin(inc)*r1
      
      vx02=-dthetadt*(-cos(node)*sin(peri+theta)-sin(node)*cos(peri+theta)*cos(inc))*r2
      vy02=-dthetadt*(-sin(node)*sin(peri+theta)+cos(node)*cos(peri+theta)*cos(inc))*r2
      vz02=-dthetadt*cos(peri+theta)*sin(inc)*r2
   else if (direction .eq. 'i') then
      vx01=-dthetadt*(-cos(node)*sin(peri+theta)-sin(node)*cos(peri+theta)*cos(inc))*r1
      vy01=-dthetadt*(-sin(node)*sin(peri+theta)+cos(node)*cos(peri+theta)*cos(inc))*r1
      vz01=-dthetadt*cos(peri+theta)*sin(inc)*r1
      
      vx02=dthetadt*(-cos(node)*sin(peri+theta)-sin(node)*cos(peri+theta)*cos(inc))*r2
      vy02=dthetadt*(-sin(node)*sin(peri+theta)+cos(node)*cos(peri+theta)*cos(inc))*r2
      vz02=dthetadt*cos(peri+theta)*sin(inc)*r2
   else  
      write(*,*), "please specify the rotation direction"
   endif
!endif
end subroutine init_rot




subroutine leapfrog(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,G,M1,M2,dt)
implicit none



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


end subroutine leapfrog
