subroutine output_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  real(dp)::lscale,vscale
  real(sp)::gamma_sp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_hydro'

  ilun=ncpu+myid+10

  lscale=dble(icoarse_max-icoarse_min+1)/boxlen
  if(cosmo)then
     vscale=dble(icoarse_max-icoarse_min+1)/(boxlen*100.)
  else
     vscale=dble(icoarse_max-icoarse_min+1)/boxlen
  end if
  gamma_sp=gamma
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)nvar
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)gamma_sp
  ! Output primitive variables
  do ilevel=1,nlevelmax
     write(ilun)ilevel
     write(ilun)numbl(myid,ilevel)
     if(numbl(myid,ilevel)>0)then
        ncache=numbl(myid,ilevel)
        allocate(ind_grid(1:ncache),xsp(1:ncache),xdp(1:ncache))
        ! Loop over level grids
        igrid=headl(myid,ilevel)
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              ind_grid(i)=ind_grid(i)+iskip
           end do
           ! Output density (in units of mean background density)
           do i=1,ncache
              xsp(i)=uold(ind_grid(i),1)
           end do
           write(ilun)xsp
           ! Output velocity (in comoving km s-1)
           do ivar=2,ndim+1
              do i=1,ncache
                 xsp(i)=uold(ind_grid(i),ivar)/uold(ind_grid(i),1)/vscale
              end do
              write(ilun)xsp
           end do
           ! Compute kinetic energy
           xdp=0.0
           do ivar=2,ndim+1
              do i=1,ncache
                 xdp(i)=xdp(i)+0.5*uold(ind_grid(i),ivar)**2 &
                      &           /uold(ind_grid(i),1)
              end do
           end do
           ! Output pressure (in consistant units)
           do i=1,ncache
              xdp(i)=uold(ind_grid(i),ndim+2)-xdp(i)
              xdp(i)=xdp(i)-uold(ind_grid(i),ndim+4)
              xdp(i)=xdp(i)/uold(ind_grid(i),ndim+3)
           end do
           do i=1,ncache
              xsp(i)=xdp(i)/vscale**2 ! Pressure
           end do
           write(ilun)xsp
           ! Output stiffen gas parameters
           do i=1,ncache
              xdp(i)=1.0d0/uold(ind_grid(i),ndim+3)+1.0d0
           end do
           do i=1,ncache
              xsp(i)=xdp(i) ! Gamma
           end do
           write(ilun)xsp
           do i=1,ncache
              gamma=1.0d0/uold(ind_grid(i),ndim+3)+1.0d0
              xdp(i)=uold(ind_grid(i),ndim+4)/gamma/uold(ind_grid(i),ndim+3)
           end do
           do i=1,ncache
              xsp(i)=xdp(i)/vscale**2 ! Pinf
           end do
           write(ilun)xsp
           ! Output remaining primitive variables
           do ivar=ndim+5,nvar
              do i=1,ncache
                 xsp(i)=uold(ind_grid(i),ivar)/uold(ind_grid(i),1)
              end do
              write(ilun)xsp
           end do
           do i=1,ncache
              ind_grid(i)=ind_grid(i)-iskip
           end do
        end do
        deallocate(ind_grid,xsp,xdp)
     end if
  end do
  close(ilun)

end subroutine output_hydro

subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  real(dp)::lscale,vscale
  real(sp)::gamma_sp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)nvar
  do ilevel=1,nlevelmax
     write(ilun)ilevel
     write(ilun)numbl(myid,ilevel)
     if(numbl(myid,ilevel)>0)then
        ncache=numbl(myid,ilevel)
        allocate(ind_grid(1:ncache),xdp(1:ncache))
        ! Loop over level grids
        igrid=headl(myid,ilevel)
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ncache
              ind_grid(i)=ind_grid(i)+iskip
           end do
           do ivar=1,nvar
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i),ivar)
              end do
              write(ilun)xdp
           end do
           do i=1,ncache
              ind_grid(i)=ind_grid(i)-iskip
           end do
        end do
        deallocate(ind_grid, xdp)
     end if
  end do
  close(ilun)
     
end subroutine backup_hydro





