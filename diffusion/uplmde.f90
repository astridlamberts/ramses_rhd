!################################################################
!################################################################
!################################################################ 
!################################################################
real(dp) function condana(d,T)
  use amr_commons
  use hydro_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the coefficient of thermal conductivity
  ! kappa as a function of mass density and electron temperature.
  ! The energy flux is given by F = -kappa grad T.
  ! Units are supposed to be in cgs here (as in units.f90)
  !--------------------------------------------------------------
  real(dp)::d,T
 
  condana=chi_params(1)*d**chi_params(2)*T**chi_params(3)

end function condana
!################################################################
!################################################################
!################################################################ 
!################################################################
real(dp) function relaxana(d,T)
  use amr_commons
  use hydro_commons
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the electron-ion equipartition rate 
  ! omega_ei as a function of mass density and electron 
  ! temperature. omega_ei = 1/tau_ei where tau_ei is the 
  ! electron-ion equipartition time in seconds.
  ! Units are supposed to be in cgs here (as in units.f90)
  !--------------------------------------------------------------
  real(dp)::d,T
 
  relaxana=omega_params(1)*d**omega_params(2)*T**omega_params(3)

end function relaxana
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_cg
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !=========================================================
  ! Iterative solver with Conjugate Gradient method 
  ! to solve A x = b
  !   r     : stored in unew(i,1)
  !   p     : stored in unew(i,2)
  ! Ap & z  : stored in unew(i,3)
  !   x     : stored in uold(i,ndim+3)
  !   b     : stored in unew(i,ndim+3)
  !=========================================================
  integer::i,idim,info,ind,iter,iskip,itermax,ilevel,ivar
  real(dp)::error,error_ini,epsilon,sum,ntot
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  real(dp)::r2_old,alpha_cg,beta_cg,etot_norm
  real(dp)::r2,pAp,rhs_norm,r2_all,pAp_all,rhs_norm_all,r3,r3_all,nleaf_all,nleaf_tot
  character(LEN=80)::filename
 
  if(verbose)write(*,111)

  ! Set constants
  epsilon=epsilon_diff

  !====================================
  ! Convert from rho.T to T
  !====================================
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              uold(active(ilevel)%igrid(i)+iskip,ndim+3)=uold(active(ilevel)%igrid(i)+iskip,ndim+3)/uold(active(ilevel)%igrid(i)+iskip,1)
              hilbert_key(active(ilevel)%igrid(i)+iskip)=uold(active(ilevel)%igrid(i)+iskip,ndim+3)
           end if
        end do
     end do
  end do
  !===================================================
  ! Compute thermal coefficient :
  ! Conductivity         : stored in divu(indcell(i))
  ! Collisions frequency : stored in enew(indcell(i)) 
  !===================================================
  call cmp_thermal_coefficient
  ! Update boundaries
  do ilevel=levelmin,nlevelmax
     call make_virtual_fine_dp(uold(1,ndim+3),ilevel)
     call make_virtual_fine_dp(hilbert_key(1),ilevel)
     call make_virtual_fine_dp(enew(1),ilevel)
     call make_virtual_fine_dp(divu(1),ilevel)
  end do
  if(simple_boundary)call make_boundary_diffusion
  
  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0.0
  etot_norm=0.0
  nleaf_tot=0
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              rhs_norm =rhs_norm +(uold(active(ilevel)%igrid(i)+iskip,ndim+3)*0.5**(ilevel)/(gamma-1.0))**2
              etot_norm=etot_norm+(uold(active(ilevel)%igrid(i)+iskip,ndim+3)*0.5**(ilevel)/(gamma-1.0))
              nleaf_tot=nleaf_tot+1
           end if
        end do
     end do
  end do

  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nleaf_tot,nleaf_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
  nleaf_tot=nleaf_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(nleaf_tot))

  !==============================================
  ! Compute r = b - Ax and store it into unew(i,1)
  ! Also set p = r and store it into unew(i,2)
  !==============================================
  call cmp_residual_diffusion_cg
  
  !====================================
  ! Main iteration loop
  !====================================
  
  iter=0; itermax=1000
  error=1.0D0; error_ini=1.0D0
  
  ! do while(error>epsilon*rhs_norm.and.iter<itermax)
  do while(error>epsilon*error_ini.and.iter<itermax)
     
     iter=iter+1
    
     !============================================
     ! Compute z = Mr and stored it into unew(i,3)
     ! with Preconditionner M=1/diag(A)
     !============================================
     call cmp_precond_diffusion
     do ilevel=levelmin,nlevelmax
        call make_virtual_fine_dp(unew(1,3),ilevel)
     end do

     !====================================
     ! Compute scalar r.z
     !====================================
     r2=0.0d0  
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                 r2=r2+unew(active(ilevel)%igrid(i)+iskip,1)*unew(active(ilevel)%igrid(i)+iskip,3)
              end if
           end do
        end do
     end do
     ! Compute global norm
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(r2,r2_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     r2=r2_all
#endif
     
     !====================================
     ! Compute residual norm
     !====================================
     r3=0.0d0
     etot_norm=0.0
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                 r3=r3+(unew(active(ilevel)%igrid(i)+iskip,1))**2
                 etot_norm=etot_norm+(unew(active(ilevel)%igrid(i)+iskip,1))
              end if
           end do
        end do
     end do 
     ! Compute global norm
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(r3,r3_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     r3=r3_all
#endif
     
     !====================================
     ! Compute beta factor
     !====================================
     if(iter==1)then
        beta_cg=0.
     else
        beta_cg=r2/r2_old
     end if
     r2_old=r2
     
     !====================================
     ! Recurrence on p = z + beta*p
     !====================================
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                 unew(active(ilevel)%igrid(i)+iskip,2)=unew(active(ilevel)%igrid(i)+iskip,3)+beta_cg*unew(active(ilevel)%igrid(i)+iskip,2)
              end if
           end do
        end do
     end do     
     ! Update boundaries
     do ilevel=levelmin,nlevelmax
        call make_virtual_fine_dp(unew(1,2),ilevel)
     end do
     if(simple_boundary)call make_boundary_diffusion
   
     !==============================================
     ! Compute q = Ap and store it into unew(i,3)
     !==============================================
     call cmp_Ap_diffusion_cg

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                 pAp=pAp+ unew(active(ilevel)%igrid(i)+iskip,2)*unew(active(ilevel)%igrid(i)+iskip,3)
              end if
           end do
        end do
     end do 
     ! Compute global sum
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(pAp,pAp_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     pAp=pAp_all
#endif

     !====================================
     ! Compute alpha factor
     !====================================
     alpha_cg = r2/pAp

     !====================================
     ! Recurrence on x = x + alpha*p
     !====================================
     ivar=ndim+3
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
              if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                 uold(active(ilevel)%igrid(i)+iskip,ivar)=uold(active(ilevel)%igrid(i)+iskip,ivar)+alpha_cg*unew(active(ilevel)%igrid(i)+iskip,2)
              end if
           end do
        end do
     end do

     !====================================
     ! Recurrence on r
     !====================================
     do ilevel=levelmin,nlevelmax
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(ilevel)%ngrid
             if(son(active(ilevel)%igrid(i)+iskip) == 0)then
                unew(active(ilevel)%igrid(i)+iskip,1)=unew(active(ilevel)%igrid(i)+iskip,1)-alpha_cg*unew(active(ilevel)%igrid(i)+iskip,3)
             end if
          end do
        end do
     end do

     error=SQRT(r3)
!    if(nstep_coarse > 1500)write(10,*)iter,error/error_ini
     if(iter==1)error_ini=error
     if(verbose) write(*,112)iter,error/rhs_norm,error/error_ini,etot_norm
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do
  ! End main iteration loop
   
! close(10)
  if(myid==1)write(*,115)ilevel-1,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Diffusion fail to converge...'
  end if

  !=============================
  ! Update energy value
  !=============================
  if(.not.static) then
     call cmp_total_energy
  end if

  !====================================
  ! Convert from T to rho.T
  !====================================
  ivar=ndim+3
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              uold(active(ilevel)%igrid(i)+iskip,ivar)=uold(active(ilevel)%igrid(i)+iskip,ivar)*uold(active(ilevel)%igrid(i)+iskip,1)
           end if
        end do
     end do
  end do

  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     ! Update boundaries
     call make_virtual_fine_dp(uold(1,ndim+3),ilevel)
     call make_virtual_fine_dp(uold(1,ndim+2),ilevel)
  end do
!if(simple_boundary)call make_boundary_hydro(ilevel)


111 format('   Entering diffusion_cg')
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x),e23.15)
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine diffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_diffusion_cg
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Diffusion The residual is stored in unew(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::i3,j3,k3,i0,j0,k0
  integer::id1,id2,ig1,ig2,ih1,ih2,nx_loc
  real(dp)::dx2,scale,oneoverfour,dx2_loc,dx,dx_loc,surf_loc,vol_loc
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:ndim),save::nug,nud
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right
  real(dp),dimension(1:nvector),save::residu
  real(dp)::Cg,Cd,Cv,enewold,nu_d,nu_c,nu_g,nu_harmo,usquare,rho,wdt,eps
  real(dp),dimension(1:3)::skip_loc
 
  logical,dimension(1:nvector,1:ndim),save::exist_nborg,exist_nbord

  ! Set constants
  Cv=1.0D0/(gamma-1.0D0)
  
  !=================================================
  ! Set unew(i,1) to 0 for leaf cells
  !=================================================
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              unew(active(ilevel)%igrid(i)+iskip,1)=0.0
           end if
        end do
     end do
  end do

  do ilevel=nlevelmax,levelmin,-1
    
     dx=0.5D0**ilevel
       
     ! Rescaling factors
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     surf_loc=dx_loc**(ndim-1)
     
     iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
     iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)


     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        
     ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        
        ! Gather neighboring grids
        do i=1,ngrid
           igridn(i,0)=ind_grid(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
              ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
              igridn(i,2*idim-1)=son(ind_left (i,idim))
              igridn(i,2*idim  )=son(ind_right(i,idim))
              ! Check if grids sits near boundaries
              exist_nborg(i,idim)=son(ind_left (i,idim))>0
              exist_nbord(i,idim)=son(ind_right(i,idim))>0
           end do
        end do
       
        ! Loop over cells
        do ind=1,twotondim

           ! Compute central cell index
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather neighboring temperature
           do idim=1,ndim
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then                    
                    if(igridn(i,ig1)>0)then
                       if(son(igridn(i,ig1)+ih1)>0)then
                          phig(i,idim)=0.0
                       else
                          phig(i,idim)=uold(igridn(i,ig1)+ih1,ndim+3)
                       endif
                    else 
                       phig(i,idim)=uold(ind_left(i,idim),ndim+3)
                    end if
                 
                 end if
              end do
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          phid(i,idim)=0.0
                       else
                          phid(i,idim)=uold(igridn(i,ig2)+ih2,ndim+3)
                       endif
                    else 
                       phid(i,idim)=uold(ind_right(i,idim),ndim+3)
                    end if

                 end if
              end do
           end do
              
           ! Compute residual using 6 neighbors temperature
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho = uold(ind_cell(i),1)
                 wdt = enew(ind_cell(i))*dtnew(ilevel)
                 residu(i)=Cv*rho*vol_loc*z_plasma*(1.0D0+(wdt/(1.0+wdt)))*uold(ind_cell(i),ndim+3)/a_plasma
              end if
           end do
            
           do idim=1,ndim
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig1)>0)then 
                       if(son(igridn(i,ig1)+ih1)>0)then
                          nug(i,idim)=0.0
                       else
                          nug(i,idim)=divu(igridn(i,ig1)+ih1)
                       end if
                    else 
                       nug(i,idim)=divu(ind_left(i,idim))
                    end if
                 end if
              end do
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          nud(i,idim) =0.0
                       else
                          nud(i,idim)=divu(igridn(i,ig2)+ih2)
                       end if
                    else 
                       nud(i,idim)=divu(ind_right(i,idim))
                    end if
                 end if
              end do
           end do  
       
              
           do idim=1,ndim      
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax          
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
                 
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    nu_c=divu(ind_cell(i))
                    nu_d=nud(i,idim)
                    nu_g=nug(i,idim)
                    
                    if(igridn(i,ig1)> 0)then
                       if(son(igridn(i,ig1)+ih1)>0)then
                          Cg=0.0
                       else
                          Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
                       end if
                    else 
                       Cg=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim ))/1.5
                    end if
                    
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          Cd=0.0
                       else
                          Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
                       end if
                    else
                       Cd=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))/1.5
                    end if
                    
                    Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                    Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                    
                    residu(i)=residu(i)+(Cg+Cd)*uold(ind_cell(i),ndim+3)-Cg*phig(i,idim)-Cd*phid(i,idim)
                 end if
              end do
           end do
              
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 usquare=0.0
                 do idim=1,ndim
                    usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
                 end do
                 
                 rho = uold(ind_cell(i),1)
                 wdt = enew(ind_cell(i))*dtnew(ilevel)
                 eps = (uold(ind_cell(i),ndim+2)-rho*usquare/2.0D0)/rho
                 
                 residu(i)=Cv*rho*vol_loc*z_plasma*(1.0D0+(1.0-z_plasma)*wdt)/(a_plasma*(1.0D0+wdt))*uold(ind_cell(i),ndim+3) &
                      &    + wdt/(1.0+wdt)*eps*z_plasma*rho*vol_loc - residu(i)
              end if
           end do

           ! Store results in unew(i,1)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),1)= unew(ind_cell(i),1) + residu(i)
              end if
           end do
           
           ! Store results in unew(i,2)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),2)= unew(ind_cell(i),1)
              end if
           end do
           
        end do
        ! End loop over cells

        ! Update residual at level ilevel-1
        i3min=0; i3max=1
        j3min=0; j3max=0
        if(ndim>1)j3max=1
        k3min=0; k3max=0
        if(ndim>2)k3max=1

        ! Loop over dimensions
        do idim=1,ndim
           i0=0; j0=0; k0=0
           if(idim==1)i0=1
           if(idim==2)j0=1
           if(idim==3)k0=1

           ! Loop over left boundary cells
           do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nborg(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_g=divu(ind_left(i,idim))
                       Cd=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim))
                       Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_left(i,idim),1)=unew(ind_left(i,idim),1) + Cd*(uold(ind_cell(i),ndim+3)-uold(ind_left(i,idim),ndim+3))/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over left boundary cells
              
           ! Loop over right boundary cells
           do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nbord(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_d=divu(ind_right(i,idim))
                       Cg=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))
                       Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_right(i,idim),1)=unew(ind_right(i,idim),1) - Cg*(uold(ind_right(i,idim),ndim+3)-uold(ind_cell(i),ndim+3))/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over right boundary cells
              
        end do
        ! End loop over dimensions

     end do
     ! End loop over grids

  end do
  ! End loop over levels

end subroutine cmp_residual_diffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_diffusion_cg
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into unew(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::i3,j3,k3,i0,j0,k0
  integer::id1,id2,ig1,ig2,ih1,ih2, nx_loc
  real(dp)::oneoverfour,dx2,nu_harmo,dx2_loc,scale,dx,dx_loc,surf_loc,vol_loc
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:ndim),save::nug,nud
  real(dp),dimension(1:nvector),save::residu
  real(dp)::Cv,Cg,Cd,nu_c,nu_d,nu_g,rho,wdt
  real(dp),dimension(1:3)::skip_loc

  logical,dimension(1:nvector,1:ndim),save::exist_nborg,exist_nbord
 
  ! Set constants
  Cv=1.0D0/(gamma-1.0D0)
  
  !=================================================
  ! Set unew(i,3) to 0 for leaf cells
  !=================================================
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              unew(active(ilevel)%igrid(i)+iskip,3)=0.0
           end if
        end do
     end do
  end do

  do ilevel=nlevelmax,levelmin,-1
    
     dx=(0.5D0**ilevel)

     ! Rescaling factors
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     surf_loc=dx_loc**(ndim-1)
     vol_loc=dx_loc**ndim
     
     iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
     iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        
        ! Gather neighboring grids
        do i=1,ngrid
           igridn(i,0)=ind_grid(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
              ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
              igridn(i,2*idim-1)=son(ind_left (i,idim))
              igridn(i,2*idim  )=son(ind_right(i,idim))
              ! Check if grids sits near boundaries
              exist_nborg(i,idim)=son(ind_left (i,idim))>0
              exist_nbord(i,idim)=son(ind_right(i,idim))>0
           end do
        end do
       
        ! Loop over cells
        do ind=1,twotondim
           
           ! Compute central cell index
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           
           
           ! Gather neighboring p value
           do idim=1,ndim
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig1)>0)then
                       if(son(igridn(i,ig1)+ih1)>0)then
                          phig(i,idim)=0.0
                       else
                          phig(i,idim)=unew(igridn(i,ig1)+ih1,2)
                       endif
                    else 
                       phig(i,idim)=unew(ind_left(i,idim),2)
                    end if
                 end if
              end do
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          phid(i,idim)=0.0
                       else
                          phid(i,idim)=unew(igridn(i,ig2)+ih2,2)
                       endif
                    else 
                       phid(i,idim)=unew(ind_right(i,idim),2)
                    end if
                 end if
              end do
           end do
           
           ! Compute Ap using neighbors p values
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho = uold(ind_cell(i),1)
                 wdt = enew(ind_cell(i))*dtnew(ilevel)
                 residu(i)=Cv*rho*vol_loc*z_plasma*(1.0D0+(wdt/(1.0D0+wdt)))*unew(ind_cell(i),2)/a_plasma
              end if
           end do

           do idim=1,ndim
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig1)>0)then 
                       if(son(igridn(i,ig1)+ih1)>0)then
                           nug(i,idim)=0.0
                       else
                          nug(i,idim)=divu(igridn(i,ig1)+ih1)
                       end if
                    else 
                       nug(i,idim)=divu(ind_left(i,idim))
                    end if
                 end if
              end do
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          nud(i,idim) =0.0
                       else
                          nud(i,idim)=divu(igridn(i,ig2)+ih2)
                       end if
                    else 
                       nud(i,idim)=divu(ind_right(i,idim))
                    end if                 
                 end if
              end do
           end do
              
           do idim=1,ndim      
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax          
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    nu_c=divu(ind_cell(i))
                    nu_d=nud(i,idim)
                    nu_g=nug(i,idim)
                    if(igridn(i,ig1)>0)then
                       if(son(igridn(i,ig1)+ih1)>0)then
                          Cg=0.0
                       else
                          Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
                       end if
                    else 
                       Cg=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim ))/1.5
                    end if
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          Cd=0.0
                       else
                          Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
                       end if
                    else
                       Cd=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))/1.5
                    end if                   
                    Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                    Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                    residu(i)=residu(i)+(Cg+Cd)*unew(ind_cell(i),2)-Cg*phig(i,idim)-Cd*phid(i,idim)
                 end if
              end do
           end do
           
           ! Store results in unew(i,3)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),3)= unew(ind_cell(i),3)+residu(i)
              end if
           end do
           
        end do
        ! End loop over cells
        
        ! Update p at level ilevel-1
        i3min=0; i3max=1
        j3min=0; j3max=0
        if(ndim>1)j3max=1
        k3min=0; k3max=0
        if(ndim>2)k3max=1

        ! Loop over dimensions
        do idim=1,ndim
           i0=0; j0=0; k0=0
           if(idim==1)i0=1
           if(idim==2)j0=1
           if(idim==3)k0=1

           ! Loop over left boundary cells
           do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nborg(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_g=divu(ind_left(i,idim))
                       Cd=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim))
                       Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_left(i,idim),3)=unew(ind_left(i,idim),3) - Cd*(unew(ind_cell(i),2)-unew(ind_left(i,idim),2))/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over left boundary cells
              
           ! Loop over right boundary cells
           do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nbord(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_d=divu(ind_right(i,idim))
                       Cg=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))
                       Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_right(i,idim),3)=unew(ind_right(i,idim),3) + Cg*(unew(ind_right(i,idim),2)-unew(ind_cell(i),2))/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over right boundary cells
              
        end do
        ! End loop over dimensions

     end do
     ! End loop over grids

  end do
  ! End loop over levels
  
end subroutine cmp_Ap_diffusion_cg
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_precond_diffusion
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes z = r/diag(a) for the Conjugate Gradient
  ! Diffusion Solver and store the result into unew(i,3).
  ! r is stored into unew(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::i3,j3,k3,i0,j0,k0
  integer::id1,id2,ig1,ig2,ih1,ih2,nx_loc
  real(dp)::oneoverfour,dx2,nu_harmo,dx2_loc,scale,dx,dx_loc,surf_loc,vol_loc
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:ndim),save::nug,nud
  real(dp),dimension(1:nvector),save::residu
  real(dp)::Cv,Cg,Cd,nu_c,nu_d,nu_g,rho,wdt
  real(dp),dimension(1:3)::skip_loc

  logical,dimension(1:nvector,1:ndim),save::exist_nborg,exist_nbord

  ! Set constants
  Cv=1.0D0/(gamma-1.0D0)

  !=================================================
  ! Set unew(i,3) to 0 for leaf cells
  !=================================================
  do ilevel=levelmin,nlevelmax
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              unew(active(ilevel)%igrid(i)+iskip,3)=0.0
           end if
        end do
     end do
  end do
     
  do ilevel=nlevelmax,levelmin,-1

     dx=0.5D0**ilevel
     
     ! Rescaling factors
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     surf_loc=dx_loc**(ndim-1)
     vol_loc=dx_loc**ndim

     iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
     iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
     iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
     iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
     
     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector

        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        
        ! Gather neighboring grids
        do i=1,ngrid
           igridn(i,0)=ind_grid(i)
        end do
        do idim=1,ndim
           do i=1,ngrid
              ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
              ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
              igridn(i,2*idim-1)=son(ind_left (i,idim))
              igridn(i,2*idim  )=son(ind_right(i,idim))
              
              ! Check if grids sits near boundaries
              exist_nborg(i,idim)=son(ind_left (i,idim))>0
              exist_nbord(i,idim)=son(ind_right(i,idim))>0
           end do
        end do
        
        
        ! Loop over cells
        do ind=1,twotondim
           
           ! Compute central cell index
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           
           
           ! Gather neighboring conduction coefficients
           do idim=1,ndim
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig1)>0)then 
                       nug(i,idim)=divu(igridn(i,ig1)+ih1)
                    else 
                       nug(i,idim)=divu(ind_left(i,idim))
                    end if
                 end if
              end do
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    if(igridn(i,ig2)>0)then
                       nud(i,idim)=divu(igridn(i,ig2)+ih2)
                    else 
                       nud(i,idim)=divu(ind_right(i,idim))
                    end if
                 end if
              end do
           end do
              
           ! Compute diagonal using neighbors conduction coefficients
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho = uold(ind_cell(i),1)
                 wdt = enew(ind_cell(i))*dtnew(ilevel)
                 residu(i)=Cv*rho*vol_loc*z_plasma*(1.0D0+(wdt/(1.0D0+wdt)))/a_plasma
              end if
           end do
               
           do idim=1,ndim      
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax          
              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax
                    
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    nu_c=divu(ind_cell(i))
                    nu_d=nud(i,idim)
                    nu_g=nug(i,idim)                    
                    if(igridn(i,ig1)>0)then
                       if(son(igridn(i,ig1)+ih1)>0)then
                          Cg=0.0
                       else
                          Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
                       end if
                    else 
                       Cg=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim ))/1.5
                    end if
                    if(igridn(i,ig2)>0)then
                       if(son(igridn(i,ig2)+ih2)>0)then
                          Cd=0.0
                       else
                          Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
                       end if
                    else
                       Cd=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))/1.5
                    end if
                    Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                    Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                    residu(i)=residu(i)+(Cg+Cd)
                 end if
              end do
              
           end do
           ! End loop over dimensions

           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),3)=unew(ind_cell(i),3)+residu(i)
              end if
           end do
              
           ! Store results in unew(i,3)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),3)=unew(ind_cell(i),1)/unew(ind_cell(i),3)
              end if
           end do
        end do
        ! End loop over cells

        ! Update diag at level ilevel-1
        i3min=0; i3max=1
        j3min=0; j3max=0
        if(ndim>1)j3max=1
        k3min=0; k3max=0
        if(ndim>2)k3max=1

        ! Loop over dimensions
        do idim=1,ndim
           i0=0; j0=0; k0=0
           if(idim==1)i0=1
           if(idim==2)j0=1
           if(idim==3)k0=1

           ! Loop over left boundary cells
           do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nborg(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_g=divu(ind_left(i,idim))
                       Cd=nu_harmo(nu_g,nu_c,ind_cell(i),ind_left(i,idim))
                       Cd=Cd*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_left(i,idim),3)=unew(ind_left(i,idim),3) + Cd/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over left boundary cells
              
           ! Loop over right boundary cells
           do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max

              ! Compute central cell index
              ind=1+i3+2*j3+4*k3
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do

              do i=1,ngrid 
                 if(son(ind_cell(i)) == 0 )then
                    if(.not. exist_nbord(i,idim))then
                       nu_c=divu(ind_cell(i))
                       nu_d=divu(ind_right(i,idim))
                       Cg=nu_harmo(nu_d,nu_c,ind_cell(i),ind_right(i,idim))
                       Cg=Cg*dtnew(ilevel)*surf_loc/dx_loc
                       unew(ind_right(i,idim),3)=unew(ind_right(i,idim),3) + Cg/1.5
                    end if
                 endif
              end do
           
           end do
           end do
           end do
           ! End loop over right boundary cells
              
        end do
        ! End loop over dimensions

     end do
     ! End loop over grids

  end do
  ! End loop over levels
  
end subroutine cmp_precond_diffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_diffusion
  use amr_commons
  use hydro_commons
  implicit none
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ilevel
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::dx,dx_loc,scale, ek_bound,condana,relaxana
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector)::cond,relax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_kappa,dd,t2
  if(.not. simple_boundary)return

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_d*scale_l**2/scale_t*scale_v**2/scale_T2

  do ilevel=levelmin,nlevelmax

     ! Mesh size at level ilevel
     dx=0.5D0**ilevel
     
     ! Rescaling factors
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     
     ! Set position of cell centers relative to grid center
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
     
     ! Loop over boundaries
     do ibound=1,nboundary
        
        ! Compute direction of reference neighbors
        boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
        if(boundary_dir==1)inbor=2
        if(boundary_dir==2)inbor=1
        if(boundary_dir==3)inbor=4
        if(boundary_dir==4)inbor=3
        if(boundary_dir==5)inbor=6
        if(boundary_dir==6)inbor=5
        
        ! Compute index of reference cells
        ! Zero flux
        if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
        if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
        if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
        if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
        if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
        if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
        ! Zero flux
        if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
        if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
        if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
        if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
        if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
        if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
        ! Imposed boundary
        if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
        if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
        if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
        if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
        if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
        if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
        
        ! Loop over grids by vector sweeps
        ncache=boundary(ibound,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
           end do
           
           ! Gather neighboring reference grid
           do i=1,ngrid
              ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
           end do
           
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              ! Gather neighboring reference cell
                 iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
                 do i=1,ngrid
                    ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
                 end do
                 
                 ! Zero flux boundary conditions
                 if((boundary_type(ibound)/10).ne.2)then
                    
                    ! Gather reference hydro variables
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          uu(i,ndim+3)=uold(ind_cell_ref(i),ndim+3)
                       end if
                    end do
                    ! Scatter to boundary region
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          uold(ind_cell(i),ndim+3)=uu(i,ndim+3)
                          hilbert_key(ind_cell(i))=uu(i,ndim+3)
                       end if
                    end do
                    
                    ! Gather reference hydro variables
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          uu(i,1)=enew(ind_cell_ref(i))
                       end if
                    end do
                    ! Scatter to boundary region
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          enew(ind_cell(i))=uu(i,1)   
                       end if
                    end do
                    
                    ! Gather reference hydro variables
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          uu(i,2)=divu(ind_cell_ref(i))
                       end if
                    end do
                    ! Scatter to boundary region
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          divu(ind_cell(i))=uu(i,2)   
                       end if
                    end do
                    
                    ! Gather reference hydro variables
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          uu(i,3)=unew(ind_cell_ref(i),2)
                       end if
                    end do
                    ! Scatter to boundary region
                    do i=1,ngrid
                       if(son(ind_cell(i)) == 0)then
                          unew(ind_cell(i),2)=uu(i,3)
                       end if
                    end do
                    
                    ! Imposed boundary conditions
                 else
                    
                    ! Compute cell center in code units
                    do idim=1,ndim
                       do i=1,ngrid
                          if(son(ind_cell(i)) == 0)then
                             xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
                          end if
                       end do
                    end do
                    
                    ! Rescale position from code units to user units
                    do idim=1,ndim
                       do i=1,ngrid
                          if(son(ind_cell(i)) == 0)then
                             xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
                          end if
                       end do
                    end do
                    
                    call boundana(xx,uu,dx_loc,ibound,ngrid)
                                   
                    ! Scatter variables
                    do i=1,ngrid 
                       if(son(ind_cell(i)) == 0)then
                          dd=max(uu(i,1),smallr)
                          t2=uu(i,ndim+3)/uu(i,1) ! Warning: divide by density was NOT perfomed before
                          uold(ind_cell(i),ndim+3)=uu(i,ndim+3)/uu(i,1)
                          hilbert_key(ind_cell(i))=uu(i,ndim+3)/uu(i,1)
                          unew(ind_cell(i),2)=0.0D0
                          t2=t2*scale_t2
                          dd=dd*scale_d
                          ! Compute thermal conductivity
                          divu(ind_cell(i))=condana(dd,t2)/scale_kappa
                          ! Compute equipartition frequency
                          enew(ind_cell(i))=relaxana(dd,t2)*scale_t
                       end if
                    end do
                 
                 end if
            
           end do
           ! End loop over cells
           
        end do
        ! End loop over grids
        
     end do
     ! End loop over boundaries
   
  end do
  ! End loop over levels

111 format('   Entering make_boundary_diffusion for level ',I2)

end subroutine make_boundary_diffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_thermal_coefficient
  use amr_commons
  use hydro_commons 
  implicit none

  integer::ilevel,igrid,ncache,ngrid,iskip
  integer::ind,idim,ivar,i
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  
  real(dp)::condana,relaxana
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_kappa,dd,t2

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_d*scale_l**2/scale_t*scale_v**2/scale_T2

  do ilevel=levelmin,nlevelmax
     ncache=active(ilevel)%ngrid
     if(ncache>0)then
        ! Loop over grids by vector sweeps
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
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0) then
                    dd=max(uold(ind_cell(i),1),smallr)
                    t2=uold(ind_cell(i),ndim+3) ! Warning: divide by density was perfomed before
                    t2=t2*scale_t2
                    dd=dd*scale_d
                    ! Compute thermal conductivity
                    divu(ind_cell(i))= condana(dd,t2)/scale_kappa
                    ! Compute equipartition frequency
                    enew(ind_cell(i))= relaxana(dd,t2)*scale_t
                 end if
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids        
     end if
  end do
  ! End loop over levels

end subroutine cmp_thermal_coefficient
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_total_energy
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ind,iskip,i,idim,nx_loc
  real(dp)::enewold,usquare,Cv,wdt,eps,ekin,teini,tiini,tfin,tefin,rho,tifin
  real(dp)::dx,dx_loc,scale,vol_loc!,boxlen

  Cv=1.0/(gamma-1.0)
  nx_loc=(icoarse_max-icoarse_min+1)   
  scale=boxlen/dble(nx_loc)
           
  do ilevel=levelmin,nlevelmax
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           if(son(active(ilevel)%igrid(i)+iskip) == 0)then
              usquare=0.0
              do idim=1,ndim
                 usquare=usquare+(uold(active(ilevel)%igrid(i)+iskip,idim+1)/uold(active(ilevel)%igrid(i)+iskip,1))**2
              end do
              rho   = uold(active(ilevel)%igrid(i)+iskip,1)
              ekin  = rho*usquare/2.0
              wdt   = enew(active(ilevel)%igrid(i)+iskip)*dtnew(ilevel)
              eps   = (uold(active(ilevel)%igrid(i)+iskip,ndim+2)-ekin)/(rho)
              teini = hilbert_key(active(ilevel)%igrid(i)+iskip)
              tiini = eps*a_plasma/Cv-teini*z_plasma
              tefin = uold(active(ilevel)%igrid(i)+iskip,ndim+3)
              tifin = 1.0/(1.0+wdt)*tiini + wdt/(1.0+wdt)*tefin
              tfin  = (tifin+z_plasma*tefin)/a_plasma
              eps   = Cv*tfin
              uold(active(ilevel)%igrid(i)+iskip,ndim+2)=rho*eps+ekin
          end if
        end do
     end do   
  end do
     
end subroutine cmp_total_energy
!################################################################
!################################################################
!################################################################ 
!################################################################
function nu_harmo(xxx,ddd,i,j)
  use hydro_commons
  use amr_commons
  implicit none
  integer ::i,j
  real(dp),INTENT(IN)::ddd,xxx
  real(dp)::nu_harmo,Tbl,Tbr


  Tbr=max(1.0d-15,hilbert_key(i))**chi_params(3)
  Tbl=max(1.0d-15,hilbert_key(j))**chi_params(3)
  nu_harmo=Tbr/ddd+Tbl/xxx
  nu_harmo=(Tbr+Tbl)/nu_harmo

  return 
end function nu_harmo
!################################################################
!################################################################
!################################################################ 
!################################################################
