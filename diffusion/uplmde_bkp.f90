!-------------------------------------------------------------
!These routines compute the thermal conductivity (condana)
!and the inverse of time relaxation (relaxana)
!-------------------------------------------------------------
subroutine condana(xx,cond,ng)
  use amr_commons
  use hydro_commons
   implicit none

  integer::ng,i
  real(dp)::opac
  real(dp),dimension(1:nvector)::cond
  real(dp),dimension(1:nvector,1:nvar)::xx
 

  if(chi_type=='analytic') then

     do i=1,ng

        cond(i)=chi_params(1) &
             &   *max(xx(i,1),1.0d-15)**chi_params(2) &
             &   *max(xx(i,ndim+3),1.0d-15)**chi_params(3)
     enddo
!!$     do i=1,ng
!!$
!!$        cond(i)=chi_params(1) &
!!$             &   *xx(i,1)**chi_params(2) &
!!$             &   *xx(i,ndim+3)**chi_params(3)
!!$     enddo
  endif

  if(chi_type=='user') then
     !-----------------------------------------
     ! Add here user's relaxation's coefficient 
     do i=1,ng
        cond(i)=1.0d0
     enddo
  endif

!!$     !------------
!!$     !Pour tube 1d
!!$     !------------
!!$     
!!$     !    cond(i)=4d-4
!!$     
!!$     !pour un bon precurseur thermique
!!$     
!!$     !    cond(i)=1d-2
!!$     
!!$     !---------------------------------
!!$     !Conductivite pour test de top hat
!!$     !---------------------------------
!!$     !cond(i)=10.0/xx(i,1)
!!$
!!$     !    if(xx(i,1)==10.0) then
!!$     !       opac=2000.0
!!$     !    else 
!!$     !       opac=0.2d+0
!!$     !    end if
!!$     
!!$     !     cond(i)=(4.0*3d+2*7.56d-15 &
!!$     !          &  *xx(i,ndim+3)**3) &
!!$     !          &  /(3.0*opac*49)
!!$     
!!$     !cond(i)=(16.0*1.03d+3*xx(i,ndim+3)**3)/(3.0*opac)
!!$     
!!$     !--------------------------------
!!$     ! Conductivite pour test de sedov

!!$! cond(i)=2.5D0*max(xx(i,ndim+3),1.0d-15)**0.5

end subroutine condana
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine relaxana(xx,relax,ng)
  use amr_commons
  use hydro_commons
  implicit none

  integer::ng,i
  real(dp),dimension(1:nvector)::relax
  real(dp),dimension(1:nvector,1:nvar)::xx

  if(omega_type=='analytic') then

     do i=1,ng
        relax(i)=omega_params(1) &
             &   *max(xx(i,1),1.0d-15)**omega_params(2) &
             &   *max(xx(i,ndim+3),1.0d-15)**omega_params(3)
     enddo
!!$     do i=1,ng
!!$        relax(i)=omega_params(1) &
!!$             &   *xx(i,1)**omega_params(2) &
!!$             &   *xx(i,ndim+3)**omega_params(3)
!!$     enddo
  endif

  if(omega_type=='user') then
     !-----------------------------------------
     ! Add here user's relaxation's coefficient 
     !-----------------------------------------
     do i=1,ng
        relax(i)=1.0
     enddo
     !.......................

  endif

!-----------------------------
!Peu de couplage, omega faible
!-----------------------------
 !    relax(i)=1.0d10

end subroutine relaxana

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_cg(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !=========================================================
  ! Iterative solver with Conjugate Gradient method 
  ! to solve A x = b
  !   r     : stored in unew(i,1)
  !   p     : stored in unew(i,2)
  ! Ap & z  : stored in unew(i,3)
  !   x     : stored in uold(i,ndim+3)
  !   b     : stored in unew(i,ndim+3)
  !=========================================================
  integer::i,idim,info,ind,iter,iskip,itermax
  real(dp)::error,error_ini,epsilon,sum
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  real(dp)::r2_old,alpha_cg,beta_cg
  real(dp)::r2,pAp,rhs_norm,r2_all,pAp_all,rhs_norm_all,r3,r3_all
  character(LEN=80)::filename
 
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  scale=boxlen/dble(icoarse_max-icoarse_min+1)
  oneoversix=1.0D0/dble(twondim)
  epsilon=epsilon_diff

  !====================================
  ! Convert from rho.T to T
  !====================================
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        uold(active(ilevel)%igrid(i)+iskip,ndim+3)= &
             & uold(active(ilevel)%igrid(i)+iskip,ndim+3)/uold(active(ilevel)%igrid(i)+iskip,1)
        hilbert_key(active(ilevel)%igrid(i)+iskip)=uold(active(ilevel)%igrid(i)+iskip,ndim+3)
     end do
  end do

  !===================================================
  ! Compute thermal coefficient :
  ! Conductivity         : stored in divu(indcell(i))
  ! Collisions frequency : stored in enew(indcell(i)) 
  !===================================================
  call cmp_thermal_coefficient(ilevel)
  ! Update boundaries
  call make_virtual_fine_dp(uold(1,ndim+3),ilevel)
  call make_virtual_fine_dp(hilbert_key(1),ilevel)
  call make_virtual_fine_dp(enew(1),ilevel)
  call make_virtual_fine_dp(divu(1),ilevel)
  if(simple_boundary)call make_boundary_diffusion(ilevel)
  
  !===============================
  ! Compute right-hand side norm
  !===============================
  rhs_norm=0.0
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rhs_norm=rhs_norm+uold(active(ilevel)%igrid(i)+iskip,ndim+3)**2
     end do
  end do
  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*numbtot(1,ilevel)))
   
  !==============================================
  ! Compute r = b - Ax and store it into unew(i,1)
  ! Also set p = r and store it into unew(i,2)
  !==============================================
  call cmp_residual_diffusion_cg(ilevel)

  !====================================
  ! Main iteration loop
  !====================================
  
  iter=0; itermax=100000
  error=1.0D0; error_ini=1.0D0
  
! do while(error>epsilon*rhs_norm.and.iter<itermax)
  do while(error>epsilon*error_ini.and.iter<itermax)

     iter=iter+1

     !====================================
     ! Compute z stored in unew(i,3)
     ! with Preconditionner M=1/diag(A)
     !====================================
     call cmp_precond_diffusion(ilevel)
     call make_virtual_fine_dp(unew(1,3),ilevel)
 
     !====================================
     ! Compute scalar r.z
     !====================================
     r2=0.0d0
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           r2=r2+unew(active(ilevel)%igrid(i)+iskip,1)*unew(active(ilevel)%igrid(i)+iskip,3)
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
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           r3=r3+unew(active(ilevel)%igrid(i+iskip),1)**2
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
     ! Recurrence on p
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,2)=unew(active(ilevel)%igrid(i)+iskip,3) &
                &              +beta_cg*unew(active(ilevel)%igrid(i)+iskip,2)
        end do
     end do
     ! Update boundaries
     call make_virtual_fine_dp(unew(1,2),ilevel)
     if(simple_boundary)call make_boundary_diffusion(ilevel)

     !==============================================
     ! Compute q = Ap and store it into unew(i,3)
     !==============================================
     call cmp_Ap_diffusion_cg(ilevel)

     !====================================
     ! Compute p.Ap scalar product
     !====================================
     pAp=0.0d0
 
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           pAp=pAp+unew(active(ilevel)%igrid(i)+iskip,2)*unew(active(ilevel)%igrid(i)+iskip,3)
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
     ! Recurrence on x
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,ndim+3)= &
                &           uold(active(ilevel)%igrid(i)+iskip,ndim+3) &
                & +alpha_cg*unew(active(ilevel)%igrid(i)+iskip,2)
        end do
     end do

     !====================================
     ! Recurrence on r
     !====================================
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,1)=unew(active(ilevel)%igrid(i)+iskip,1) &
                &             -alpha_cg*unew(active(ilevel)%igrid(i)+iskip,3)
        end do
     end do
     ! Compute error
     error=DSQRT(r3/dble(twotondim*numbtot(1,ilevel)))
     if(iter==1)error_ini=error
     if(verbose)write(*,112)iter,error/rhs_norm,error/error_ini

  end do

  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Diffusion failed to converge...'
  end if

  !=============================
  ! Update energy value
  !=============================
  if(.not.static) then
     call cmp_total_energy(ilevel)
  end if

  !====================================
  ! Convert from T to rho.T
  !====================================
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        uold(active(ilevel)%igrid(i)+iskip,ndim+3)= &
             &        uold(active(ilevel)%igrid(i)+iskip,ndim+3) &
             &       *uold(active(ilevel)%igrid(i)+iskip,1)
     end do
  end do

  ! Update boundaries
  call make_virtual_fine_dp(uold(1,ndim+3),ilevel)
  call make_virtual_fine_dp(uold(1,ndim+2),ilevel)
!if(simple_boundary)call make_boundary_hydro(ilevel)

111 format('   Entering diffusion_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine diffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_diffusion_cg(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in unew(i,1).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2,nx_loc
  real(dp)::dx2,scale,oneoversix,dx_loc,dx2_loc
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


  ! Set constants
  dx2=(0.5D0**ilevel)**2
  scale=boxlen/dble(icoarse_max-icoarse_min+1)
  Cv=1.0D0/(gamma-1.0D0)

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx2_loc=dx2*scale*scale

  !  if(cosmo)scale=scale/boxlen
  !  fourpi=4.D0*ACOS(-1.0D0)*scale
  !  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  !  fact=oneoversix*fourpi*dx2

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
        end do
     end do
     
     ! Interpolate potential from upper level
 !    do idim=1,ndim
 !       call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel)
 !       call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel)
 !    end do


     ! Loop over cells
     do ind=1,twotondim
        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=uold(igridn(i,ig1)+ih1,ndim+3)
              else
                 phig(i,idim)=phi_left(i,id1,idim)
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=uold(igridn(i,ig2)+ih2,ndim+3)
              else
                 phid(i,idim)=phi_right(i,id2,idim)
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Compute residual using 6 neighbors potential
        do i=1,ngrid
           rho = uold(ind_cell(i),1)
           wdt = enew(ind_cell(i))*dtnew(ilevel)
           residu(i)=Cv*rho*(1.0D0+(wdt/(1.0+wdt)))*uold(ind_cell(i),ndim+3)
        end do

        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              nug(i,idim)=divu(igridn(i,ig1)+ih1)
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              nud(i,idim)=divu(igridn(i,ig2)+ih2)
           end do
        end do
        
        do idim=1,ndim      
           do i=1,ngrid

              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax

              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax

              nu_c=divu(ind_cell(i))
              nu_d=nud(i,idim)
              nu_g=nug(i,idim)

              Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
              Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
           
              Cg=Cg*dtnew(ilevel)/dx2_loc
              Cd=Cd*dtnew(ilevel)/dx2_loc

              residu(i)=residu(i)+(Cg+Cd) &
                   &    *uold(ind_cell(i),ndim+3) &
                   &    -Cg*phig(i,idim) &
                   &    -Cd*phid(i,idim)
           end do
        end do

        do i=1,ngrid
           usquare=0.0
           do idim=1,ndim
              usquare=usquare+(uold(ind_cell(i),idim+1)/uold(ind_cell(i),1))**2
           end do

           rho = uold(ind_cell(i),1)
           wdt = enew(ind_cell(i))*dtnew(ilevel)
           eps = (uold(ind_cell(i),ndim+2)-rho*usquare/2.0D0)/rho

           residu(i)=Cv*rho*1.0D0/(1.0D0+wdt)*uold(ind_cell(i),ndim+3) &
                   & + wdt/(1.0+wdt)*eps*rho &
                   & -residu(i)


        end do

        ! Store results in unew(i,1)
        do i=1,ngrid
           unew(ind_cell(i),1)=residu(i)
        end do

        ! Store results in unew(i,2)
        do i=1,ngrid
           unew(ind_cell(i),2)=residu(i)
        end do
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_residual_diffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Ap_diffusion_cg(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes Ap for the Conjugate Gradient
  ! Poisson Solver and store the result into unew(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2, nx_loc
  real(dp)::oneoversix,dx2,nu_harmo,dx_loc,dx2_loc,scale
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:ndim),save::nug,nud
  real(dp),dimension(1:nvector),save::residu
  real(dp)::Cv,Cg,Cd,nu_c,nu_d,nu_g,rho,wdt
  real(dp),dimension(1:3)::skip_loc

  ! Set constants
  oneoversix=1.0D0/dble(twondim)
  dx2=(0.5D0**ilevel)**2
  Cv=1.0D0/(gamma-1.0D0)


  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx2_loc=dx2*scale*scale


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
           igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
           igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim

        ! Gather neighboring potential
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig1)>0)then
                 phig(i,idim)=unew(igridn(i,ig1)+ih1,2)
              else
                 phig(i,idim)=0.
              end if
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
                 phid(i,idim)=unew(igridn(i,ig2)+ih2,2)
              else
                 phid(i,idim)=0.
              end if
           end do
        end do

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

         do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              nug(i,idim)=divu(igridn(i,ig1)+ih1)
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              nud(i,idim)=divu(igridn(i,ig2)+ih2)
           end do
        end do

        ! Compute Ap using neighbors potential
        do i=1,ngrid
           rho = uold(ind_cell(i),1)
           wdt = enew(ind_cell(i))*dtnew(ilevel)
           residu(i)=Cv*rho*(1.0D0+(wdt/(1.0D0+wdt))) &
                &    *unew(ind_cell(i),2)
        end do

        do idim=1,ndim
           do i=1,ngrid

              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax

              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax

              nu_c=divu(ind_cell(i))
              nu_d=nud(i,idim)
              nu_g=nug(i,idim)

              Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
              Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
           
              Cg=Cg*dtnew(ilevel)/dx2_loc
              Cd=Cd*dtnew(ilevel)/dx2_loc
             

              residu(i)=residu(i)+(Cg+Cd)*unew(ind_cell(i),2) &
                   &   -Cg*phig(i,idim) &
                   &   -Cd*phid(i,idim)
           end do
        end do

        ! Store results in unew(i,3)
        do i=1,ngrid
           unew(ind_cell(i),3)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_Ap_diffusion_cg
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_diffusion(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::dx,dx_loc,scale, ek_bound
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector)::cond,relax

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

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
                 uu(i,ndim+3)=uold(ind_cell_ref(i),ndim+3)
              end do
              ! Scatter to boundary region
              do i=1,ngrid
                 uold(ind_cell(i),ndim+3)=uu(i,ndim+3)
                 hilbert_key(ind_cell(i))=uu(i,ndim+3)
              end do       
              
              ! Gather reference hydro variables
              do i=1,ngrid
                 uu(i,1)=enew(ind_cell_ref(i))
              end do
              ! Scatter to boundary region
              do i=1,ngrid
                 enew(ind_cell(i))=uu(i,1)               
              end do
             
              ! Gather reference hydro variables
              do i=1,ngrid
                 uu(i,2)=divu(ind_cell_ref(i))
              end do
              ! Scatter to boundary region
              do i=1,ngrid
                 divu(ind_cell(i))=uu(i,2)               
              end do

              ! Gather reference hydro variables
              do i=1,ngrid
                 uu(i,3)=unew(ind_cell_ref(i),2)
              end do
              ! Scatter to boundary region
              do i=1,ngrid
                 unew(ind_cell(i),2)=uu(i,3)
              end do

              call condana(uu,cond,ngrid)
              call relaxana(uu,relax,ngrid)

              do i=1,ngrid
                 divu(ind_cell(i))=cond(i)
                 enew(ind_cell(i))=relax(i)
              end do
              
           ! Imposed boundary conditions
           else
              
              ! Compute cell center in code units
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
              
              call boundana(xx,uu,dx_loc,ibound,ngrid)
              call condana(uu,cond,ngrid)
              call relaxana(uu,relax,ngrid)

              ! Scatter variables
              do i=1,ngrid
                 uold(ind_cell(i),ndim+3)=uu(i,ndim+3)/uu(i,1)
                 hilbert_key(ind_cell(i))=uu(i,ndim+3)/uu(i,1)
                 unew(ind_cell(i),2)=0.0D0
                 divu(ind_cell(i))=cond(i)
                 enew(ind_cell(i))=relax(i)
              end do
                 
           end if
              
        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_diffusion for level ',I2)

end subroutine make_boundary_diffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_thermal_coefficient(ilevel)
  use amr_commons
  use hydro_commons 
  implicit none

  integer::ilevel,igrid,ncache,ngrid,iskip
  integer::ind,idim,ivar,i
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector)::cond,relax


  ncache=active(ilevel)%ngrid
  if(ncache==0)return

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
        ! Take conservative hydro variables
        do ivar=1,nvar
           do i=1,ngrid
              uu(i,ivar)=uold(ind_cell(i),ivar)
           end do
        end do

        ! Compute thermal conductivity
        call condana(uu,cond,ngrid)
   
        ! Compute "collisions" frequency
        call relaxana(uu,relax,ngrid)
      
        do i=1,ngrid
           divu(ind_cell(i))=cond(i)
           enew(ind_cell(i))=relax(i) 
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine cmp_thermal_coefficient
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_precond_diffusion(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes z for the Conjugate Gradient
  ! Diffusion Solver and store the result into unew(i,3).
  !------------------------------------------------------------------
  integer::i,idim,igrid,ngrid,ncache,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2,nx_loc
  real(dp)::oneoversix,dx2,nu_harmo,dx2_loc,scale
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector,1:ndim),save::phig,phid
  real(dp),dimension(1:nvector,1:ndim),save::nug,nud
  real(dp),dimension(1:nvector),save::residu
  real(dp)::Cv,Cg,Cd,nu_c,nu_d,nu_g,rho,wdt
  real(dp),dimension(1:3)::skip_loc

  ! Set constants
  oneoversix=1.0D0/dble(twondim)
  dx2=(0.5D0**ilevel)**2
  Cv=1.0D0/(gamma-1.0D0)


  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx2_loc=dx2*scale*scale

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
           igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
           igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

         do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              nug(i,idim)=divu(igridn(i,ig1)+ih1)
           end do
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              nud(i,idim)=divu(igridn(i,ig2)+ih2)
           end do
        end do

        ! Compute Ap using neighbors potential
        do i=1,ngrid
           rho = uold(ind_cell(i),1)
           wdt = enew(ind_cell(i))*dtnew(ilevel)
           residu(i)=Cv*rho*(1.0D0+(wdt/(1.0D0+wdt))) 
        end do

        do idim=1,ndim
           do i=1,ngrid
              
              id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
              ih1=ncoarse+(id1-1)*ngridmax

              id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
              ih2=ncoarse+(id2-1)*ngridmax

              nu_c=divu(ind_cell(i))
              nu_d=nud(i,idim)
              nu_g=nug(i,idim)

              Cd=nu_harmo(nu_d,nu_c,ind_cell(i),igridn(i,ig2)+ih2)
              Cg=nu_harmo(nu_g,nu_c,ind_cell(i),igridn(i,ig1)+ih1)
                         
              Cg=Cg*dtnew(ilevel)/dx2_loc
              Cd=Cd*dtnew(ilevel)/dx2_loc

              residu(i)=residu(i)+(Cg+Cd)
           end do
        end do

        do i=1,ngrid
           residu(i)=unew(ind_cell(i),1)/residu(i)
        end do

        ! Store results in unew(i,3)
        do i=1,ngrid
           unew(ind_cell(i),3)=residu(i)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine cmp_precond_diffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_total_energy(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ind,iskip,i,idim
  real(dp)::enewold,usquare,CV,wdt,eps,ekin,teini,tiini,tfin,tefin,rho
  
  Cv=1.0/(gamma-1.0)

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        active(ilevel)%igrid(i)=active(ilevel)%igrid(i)+iskip
     end do


     do i=1,active(ilevel)%ngrid
        
        usquare=0.0
        do idim=1,ndim
           usquare=usquare+(uold(active(ilevel)%igrid(i),idim+1)/uold(active(ilevel)%igrid(i),1))**2
        end do

        rho   = uold(active(ilevel)%igrid(i),1)
        ekin  = rho*usquare/2.0
        wdt   = enew(active(ilevel)%igrid(i))*dtnew(ilevel)
        eps   = (uold(active(ilevel)%igrid(i),ndim+2)-ekin)/rho
        teini = hilbert_key(active(ilevel)%igrid(i))
        tiini = eps/Cv-teini
        tefin = uold(active(ilevel)%igrid(i),ndim+3)
        tfin  = 1.0/(1.0+wdt)*tiini + (1.0+2.0*wdt)/(1.0+wdt)*tefin
        eps   = Cv*tfin

        uold(active(ilevel)%igrid(i),ndim+2)=rho*eps+ekin

     end do
     do i=1,active(ilevel)%ngrid
        active(ilevel)%igrid(i)=active(ilevel)%igrid(i)-iskip
        
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

!--------------------------------------------------

