subroutine write_screen
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
  include 'mpif.h'
  !
  integer::igrid,jgrid,ind,icpu,info
  integer::i,icell,ncell,ilevel,ncache
  integer::icellmin
  real(dp)::dx,scale,smallp,ddd,ppp,rrr,eee,uuu,aaa,bbb

  integer     ,dimension(:),allocatable::ind_grid,ind_cell,ind_sort,ll,ll_all
  real(kind=8),dimension(:),allocatable::rr,et,ei,dd,uu,mm,gg,dtot,aa,bb
  real(kind=8),dimension(:),allocatable::rr_all,et_all,ei_all,aa_all,bb_all
  real(kind=8),dimension(:),allocatable::dd_all,uu_all,mm_all,gg_all,dtot_all

  integer,dimension(1:ncpu)::iskip,ncell_loc,ncell_all

  if(ndim>1)return

  call MPI_BARRIER(MPI_COMM_WORLD,info)
  
  ncell=0
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Count leaf cells
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))== 0)then
                 ncell=ncell+1
              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

  ncell_loc=0
  ncell_all=0
  ncell_loc(myid)=ncell
  call MPI_ALLREDUCE(ncell_loc,ncell_all,ncpu,MPI_INTEGER,MPI_SUM,&
       & MPI_COMM_WORLD,info)

  ncell=0
  iskip=0
  do icpu=1,ncpu
     iskip(icpu)=ncell
     ncell=ncell+ncell_all(icpu)
  end do

  if(myid==1)write(*,114)ncell

  if(ncell>0)then

  allocate(rr(1:ncell),mm(1:ncell),dd(1:ncell),dtot(1:ncell))
  allocate(et(1:ncell),ei(1:ncell),aa(1:ncell),bb(1:ncell))
  allocate(uu(1:ncell),ll(1:ncell),gg(1:ncell))
  allocate(rr_all(1:ncell),mm_all(1:ncell),dd_all(1:ncell),dtot_all(1:ncell))
  allocate(et_all(1:ncell),ei_all(1:ncell),aa_all(1:ncell),bb_all(1:ncell))
  allocate(uu_all(1:ncell),ll_all(1:ncell),gg_all(1:ncell))
  rr=0.0D0; mm=0.0D0; dd=0.0D0; dtot=0.0D0; et=0.0D0
  ei=0.0D0; uu=0.0D0; gg=0.0D0; ll=0; aa=0.0; bb=0.0
  rr_all=0.0D0; mm_all=0.0D0; dd_all=0.0D0; dtot_all=0.0D0; et_all=0.0D0
  ei_all=0.0D0; uu_all=0.0D0; gg_all=0.0D0; ll_all=0;aa_all=0.0;bb_all=0.0

  icell=iskip(myid)
  do ilevel=1,nlevelmax
     icellmin=icell
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        icell=icellmin
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 icell=icell+1
                 rr(icell)=xg(ind_grid(i),1)+(dble(ind)-1.5D0)*dx
                 ll(icell)=ilevel
              end if
           end do
        end do
        if(hydro)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1
                    dd(icell)=uold(ind_cell(i),1)
                    mm(icell)=dd(icell)
                    uu(icell)=uold(ind_cell(i),2)/dd(icell)
                    et(icell)=uold(ind_cell(i),3)
                    ei(icell)=et(icell)/dd(icell)-0.5d0*uu(icell)**2
                    ei(icell)=ei(icell)*dd(icell)
                    aa(icell)=uold(ind_cell(i),4)
                    bb(icell)=uold(ind_cell(i),5)
                 end if
              end do
           end do
        end if
        if(poisson)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1
                    dtot(icell)=rho(ind_cell(i))
                    gg(icell)=f(ind_cell(i),1)
                 end if
              end do
           end do
        end if
        deallocate(ind_grid, ind_cell)
     end if
  end do

  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mm,mm_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dd,dd_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtot,dtot_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(et,et_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ei,ei_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(uu,uu_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(gg,gg_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ll,ll_all,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(aa,aa_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(bb,bb_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rr=rr_all; mm=mm_all; dd=dd_all; dtot=dtot_all; et=et_all
  ei=ei_all; uu=uu_all; gg=gg_all; ll=ll_all; aa=aa_all; bb=bb_all


  if(myid==1)then
     ! Sort radius
     allocate(ind_sort(1:ncell))
     call quick_sort(rr,ind_sort,ncell)
     ! Write results to screen
     smallp=smallc**2/gamma
     scale=dble(icoarse_max-icoarse_min+1)/boxlen
     do i=1,ncell
        rrr=(rr(i)-dble(icoarse_min))/scale
        ddd=MAX(dd(ind_sort(i)),smallr)
        eee=ei(ind_sort(i))/scale**2
        uuu=uu(ind_sort(i))/scale
        aaa=aa(ind_sort(i))
        bbb=bb(ind_sort(i))/scale**2
        ppp=(eee-bbb)/aaa
        if(ABS(uuu)<1.d-15*scale)uuu=0.0D0
        if(ABS(bbb)<1.d-15*ppp)bbb=0.0D0
        write(*,113) &
             & ll(ind_sort(i)),  &
             & rrr, &
             & ddd, &
             & uuu, &
             & ppp, &
             & 1.0/aaa+1.0, &
             & bbb/(1.0/aaa+1.0)/aaa
     end do
     deallocate(ind_sort)
  end if

  ! Deallocate local arrays
  deallocate(mm,rr,dd,dtot,et,ei,uu,ll,gg,aa,bb,aa_all,bb_all)
  deallocate(mm_all,rr_all,dd_all,dtot_all,et_all,ei_all,uu_all,ll_all,gg_all)

  end if
 
  call MPI_BARRIER(MPI_COMM_WORLD,info)

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,9(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')

end subroutine write_screen
