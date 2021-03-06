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
  

  ! initialize the medium and mask

 lor2=(1.-vinf2**2)**(-1./2.)


  rho_a2 = Mdot2*a1**2/(4.0d0*3.1415*vinf2*lor2)
  
#if NDIM==2
  rho_a2 = Mdot2*a1/(2.0d0*3.1415*vinf2*lor2)
#endif
  
! We assume the mach number given in the log file is the classical mach number

  P_a2 = vinf2**2* rho_a2/(gamma*Mach2**2-vinf2**2*gamma/(gamma-1))


!  P_a2 = vinf2**2*lor2* rho_a2/(gamma*Mach2**2-vinf2**2*lor2**2*gamma/(gamma-1))
  P_amb=P_a2/R_P
  rho_amb=rho_a2/R_rho
  
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


 
     call init_puls(q,x,nn,vinf2,r02,Mdot2,Mach2,x02,y02,z02,index_pas_scal2)
 
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
  
! no rotation possible    
   
!inutile ici?
!values for r=a (mettre dans une fonction)
  
  !v2  = vinf2*(1.0d0-rstar2)**beta2
     lor2=(1.-vinf2**2)**(-1./2.) 
     rho2 = Mdot2*a1**2/(4.0d0*3.1415*vinf2*lor2)
  
#if NDIM==2
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
 
        call init2_puls(xx,u_mask,dx_loc,ngrid,x02,y02,z02,r02,Mdot2,vinf2,Mach2,dx,index_pas_scal2)

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

         if (rc2 .lt. 1.2d0*r02)then 
            !2nd body: a star or a pulsar
            call reset_puls(xc2,yc2,zc2,r02,dx,u_mask,i,ind_cell(i))
         endif
      end do
   end do
   ! End loop over cells
end do
! End loop over grids

 end subroutine reset_mask





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

!#if NDIM==3
!  write(*,*),vinf,'ini'
!  rho_m=  Mdot/(4.0d0*3.141592*vinf*r0**2)
!  P_m=P_a*(rho_a/rho_m)**(-gamma)
!  vinf=sqrt(vinf**2-P_m/rho_m*gamma/(gamma-1.d0))
 ! write(*,*),vinf,P_m,rho_m
!#endif

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
     &theta,x01,y01,z01,x02,y02,z02

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


