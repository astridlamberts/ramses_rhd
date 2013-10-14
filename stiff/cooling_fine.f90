subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(active(ilevel)%ngrid==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if(ilevel==levelmin)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf
  real(dp)::scale_nH,scale_T2,scale_dt
  real(kind=8)::dtcool
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,dT2dt,ekk

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute specific internal energy
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        T2(i)=T2(i)-ekk(i)
     end do

     ! Compute T2 in Kelvin
     scale_T2 = (gamma-1.0)*(100.*boxlen*1.E5/dble(nx)/aexp)**2*mH/kB
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in cm-3
     scale_nH = omega_m*rhoc*(h0/100.)**2/aexp**3*X/mH
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute cooling time step in second
     scale_dt = aexp**2/(h0*1.E5/3.08d24)
     dtcool=dtnew(ilevel)*scale_dt

     ! Compute net cooling at constant nH
     call solve_cooling(nH,T2,dtcool,dT2dt,nleaf)

!!$     ! Damp excessive cooling
!!$     do i=1,nleaf
!!$        dT2dt(i) = MAX(dT2dt(i),-T2(i)/2.0)
!!$     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=nH(i)/scale_nH
     end do

     ! Compute net energy sink in code units
     do i=1,nleaf
        dT2dt(i) = dT2dt(i)*nH(i)/scale_T2
     end do

     ! Update total fluid energy
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        T2(i)=T2(i)+dT2dt(i)
     end do
     do i=1,nleaf
        uold(ind_leaf(i),ndim+2)=T2(i)
     end do

  end do
  ! End loop over cells

end subroutine coolfine1



