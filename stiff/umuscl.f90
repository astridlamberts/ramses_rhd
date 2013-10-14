! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela's PLMDE scheme
!              with various slope limiters.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gin         => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
! ----------------------------------------------------------------
subroutine unsplit(uin,gin,flux,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gin 

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2       ),save::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq  

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp  
  
  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx  

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primative variables, compute sound speeds  
  if(debug)write(*,*)'Before ctoprim'
  call ctoprim(uin,qin,cin,gin,dt,ngrid)

  ! Compute TVD slopes
  if(debug)write(*,*)'Before uslope'
  call uslope(qin,dq,ngrid)

  ! Compute 3D traced-states in all three directions
  if(debug)write(*,*)'Before trace'
  if(scheme=='muscl')then
     if(ndim==1)call trace1d (qin,dq,cin,qm,qp,dx      ,dt,ngrid)
     if(ndim==2)call trace2d (qin,dq,cin,qm,qp,dx,dy   ,dt,ngrid)
     if(ndim==3)call trace3d (qin,dq,cin,qm,qp,dx,dy,dz,dt,ngrid)
  endif
  if(scheme=='plmde')then
     if(ndim==1)call tracex  (qin,dq,cin,qm,qp,dx      ,dt,ngrid)
     if(ndim==2)call tracexy (qin,dq,cin,qm,qp,dx,dy   ,dt,ngrid)
     if(ndim==3)call tracexyz(qin,dq,cin,qm,qp,dx,dy,dz,dt,ngrid)
  endif

  ! Solve for 1D flux in X direction
  if(debug)write(*,*)'Before cmpflxm along x'
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,fx,tx,ngrid)
  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  if(debug)write(*,*)'Before cmpflxm along y'
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,fx,tx,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
#endif

  if(debug)write(*,*)'unsplit complete'

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(q,dq,c,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! Local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip
  real(dp)::dtdx,csq
  real(dp)::r, u, p, a
  real(dp)::drx, dux, dpx, dax
  real(dp)::sr0, su0, sp0, sa0
  
  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   =  q(l,i,j,k,ir)
              u   =  q(l,i,j,k,iu)
              p   =  q(l,i,j,k,ip)
              csq =  c(l,i,j,k)**2

              ! TVD slopes in X direction
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dpx = dq(l,i,j,k,ip,1)
              
              ! Source terms (including transverse derivatives)
              sr0 = -u*drx - (dux)*r
              su0 = -u*dux - (dpx)/r
              sp0 = -u*dpx - (dux)*r*csq

              ! Right state
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))

              ! Left state
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

           end do
        end do
     end do
  end do

  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 sa0 = -u*dax             ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
              end do
           end do
        end do
     end do
  end do

end subroutine trace1d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(q,dq,c,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip
  real(dp)::dtdx, dtdy, csq
  real(dp)::r, u, v, p, a
  real(dp)::drx, dux, dvx, dpx, dax
  real(dp)::dry, duy, dvy, dpy, day
  real(dp)::sr0, su0, sv0, sp0, sa0
  
  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   = q(l,i,j,k,ir)
              u   = q(l,i,j,k,iu)
              v   = q(l,i,j,k,iv)
              p   = q(l,i,j,k,ip)
              csq = c(l,i,j,k)**2

              ! TVD slopes in all directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dpx = dq(l,i,j,k,ip,1)
              
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dpy = dq(l,i,j,k,ip,2)
              
              ! source terms (with transverse derivatives)
              sr0 = -u*drx-v*dry - (dux+dvy)*r
              su0 = -u*dux-v*duy - (dpx    )/r
              sv0 = -u*dvx-v*dvy - (dpy    )/r
              sp0 = -u*dpx-v*dpy - (dux+dvy)*r*csq

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(l,i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(l,i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(l,i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(l,i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(l,i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(l,i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))

           end do
        end do
     end do
  end do

  ! passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sa0 = -u*dax-v*day       ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
              end do
           end do
        end do
     end do
  end do

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(q,dq,c,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  integer ::i, j, k, l, n
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx, dtdy, dtdz, csq
  real(dp)::r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::sr0, su0, sv0, sw0, sp0, sa0
  
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values
              r   = q(l,i,j,k,ir)
              u   = q(l,i,j,k,iu)
              v   = q(l,i,j,k,iv)
              w   = q(l,i,j,k,iw)
              p   = q(l,i,j,k,ip)
              csq = c(l,i,j,k)**2

              ! TVD slopes in all 3 directions
              drx = dq(l,i,j,k,ir,1)
              dpx = dq(l,i,j,k,ip,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dwx = dq(l,i,j,k,iw,1)
              
              dry = dq(l,i,j,k,ir,2)
              dpy = dq(l,i,j,k,ip,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,j,k,iw,2)
              
              drz = dq(l,i,j,k,ir,3)
              dpz = dq(l,i,j,k,ip,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-v*dry-w*drz - (dux+dvy+dwz)*r
              sp0 = -u*dpx-v*dpy-w*dpz - (dux+dvy+dwz)*r*csq
              su0 = -u*dux-v*duy-w*duz - (dpx        )/r
              sv0 = -u*dvx-v*dvy-w*dvz - (dpy        )/r
              sw0 = -u*dwx-v*dwy-w*dwz - (dpz        )/r

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - half*drx + sr0*dtdx*half
              qp(l,i,j,k,ip,1) = p - half*dpx + sp0*dtdx*half
              qp(l,i,j,k,iu,1) = u - half*dux + su0*dtdx*half
              qp(l,i,j,k,iv,1) = v - half*dvx + sv0*dtdx*half
              qp(l,i,j,k,iw,1) = w - half*dwx + sw0*dtdx*half
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))

              ! Left state at left interface
              qm(l,i,j,k,ir,1) = r + half*drx + sr0*dtdx*half
              qm(l,i,j,k,ip,1) = p + half*dpx + sp0*dtdx*half
              qm(l,i,j,k,iu,1) = u + half*dux + su0*dtdx*half
              qm(l,i,j,k,iv,1) = v + half*dvx + sv0*dtdx*half
              qm(l,i,j,k,iw,1) = w + half*dwx + sw0*dtdx*half
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

              ! Top state at bottom interface
              qp(l,i,j,k,ir,2) = r - half*dry + sr0*dtdy*half
              qp(l,i,j,k,ip,2) = p - half*dpy + sp0*dtdy*half
              qp(l,i,j,k,iu,2) = u - half*duy + su0*dtdy*half
              qp(l,i,j,k,iv,2) = v - half*dvy + sv0*dtdy*half
              qp(l,i,j,k,iw,2) = w - half*dwy + sw0*dtdy*half
              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))

              ! Bottom state at top interface
              qm(l,i,j,k,ir,2) = r + half*dry + sr0*dtdy*half
              qm(l,i,j,k,ip,2) = p + half*dpy + sp0*dtdy*half
              qm(l,i,j,k,iu,2) = u + half*duy + su0*dtdy*half
              qm(l,i,j,k,iv,2) = v + half*dvy + sv0*dtdy*half
              qm(l,i,j,k,iw,2) = w + half*dwy + sw0*dtdy*half
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))

              ! Back state at front interface
              qp(l,i,j,k,ir,3) = r - half*drz + sr0*dtdz*half
              qp(l,i,j,k,ip,3) = p - half*dpz + sp0*dtdz*half
              qp(l,i,j,k,iu,3) = u - half*duz + su0*dtdz*half
              qp(l,i,j,k,iv,3) = v - half*dvz + sv0*dtdz*half
              qp(l,i,j,k,iw,3) = w - half*dwz + sw0*dtdz*half
              qp(l,i,j,k,ir,3) = max(smallr, qp(l,i,j,k,ir,3))

              ! Front state at back interface
              qm(l,i,j,k,ir,3) = r + half*drz + sr0*dtdz*half
              qm(l,i,j,k,ip,3) = p + half*dpz + sp0*dtdz*half
              qm(l,i,j,k,iu,3) = u + half*duz + su0*dtdz*half
              qm(l,i,j,k,iv,3) = v + half*dvz + sv0*dtdz*half
              qm(l,i,j,k,iw,3) = w + half*dwz + sw0*dtdz*half
              qm(l,i,j,k,ir,3) = max(smallr, qm(l,i,j,k,ir,3))

           end do
        end do
     end do
  end do

  ! Passive scalars
  do n = ndim+3, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 a   = q(l,i,j,k,n)       ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sa0 = -u*dax-v*day-w*daz     ! Source terms
                 qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half  ! Right state
                 qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half  ! Left state
                 qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half  ! Bottom state
                 qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half  ! Upper state
                 qp(l,i,j,k,n,3) = a - half*daz + sa0*dtdz*half  ! Front state
                 qm(l,i,j,k,n,3) = a + half*daz + sa0*dtdz*half  ! Back state
              end do
           end do
        end do
     end do
  end do

end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,tmp,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp)::entho
  real(dp),dimension(1:nvector),save::rhoetot
  real(dp),dimension(1:nvector,1:nvar),save::qleft,qright,qgdnv

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           
           ! Mass density
           do l = 1, ngrid
              qleft (l,1) = qm(l,i,j,k,1,xdim)
              qright(l,1) = qp(l,i,j,k,1,xdim)
           end do
           
           ! Normal velocity
           do l = 1, ngrid
              qleft (l,2) = qm(l,i,j,k,ln,xdim)
              qright(l,2) = qp(l,i,j,k,ln,xdim)
           end do
           
           ! Pressure
           do l = 1, ngrid
              qleft (l,3) = qm(l,i,j,k,ndim+2,xdim)
              qright(l,3) = qp(l,i,j,k,ndim+2,xdim)
           end do
           
           ! Alpha
           do l = 1, ngrid
              qleft (l,4) = qm(l,i,j,k,ndim+3,xdim)
              qright(l,4) = qp(l,i,j,k,ndim+3,xdim)
           end do
           
           ! Beta
           do l = 1, ngrid
              qleft (l,5) = qm(l,i,j,k,ndim+4,xdim)
              qright(l,5) = qp(l,i,j,k,ndim+4,xdim)
           end do
           
           ! Tangential velocity 1
#if NDIM>1
           do l = 1, ngrid
              qleft (l,6) = qm(l,i,j,k,lt1,xdim)
              qright(l,6) = qp(l,i,j,k,lt1,xdim)
           end do
#endif
           ! Tangential velocity 2
#if NDIM>2
           do l = 1, ngrid
              qleft (l,7) = qm(l,i,j,k,lt2,xdim)
              qright(l,7) = qp(l,i,j,k,lt2,xdim)
           end do
#endif           
           ! Other advected quantities
           do n = ndim+5, nvar
              do l = 1, ngrid
                 qleft (l,n) = qm(l,i,j,k,n,xdim)
                 qright(l,n) = qp(l,i,j,k,n,xdim)
              end do
           end do
          
           ! Solve Riemann problem
           if(riemann.eq.'acoustic')then
              call riemann_acoustic(qleft,qright,qgdnv,ngrid)
           else
              call riemann_approx(qleft,qright,qgdnv,ngrid)
           end if
           
           ! Compute fluxes
           
           ! Mass density
           do l = 1, ngrid 
              flx(l,i,j,k,1) = qgdnv(l,1)*qgdnv(l,2)
           end do
           
           ! Normal momentum
           do l = 1, ngrid
              flx(l,i,j,k,ln) = flx(l,i,j,k,1)*qgdnv(l,2)+qgdnv(l,3)
           end do

           ! Transverse momentum 1
#if NDIM>1
           do l = 1, ngrid
              flx(l,i,j,k,lt1) = flx(l,i,j,k,1)*qgdnv(l,6)
           end do
#endif

           ! Transverse momentum 2
#if NDIM>2
           do l = 1, ngrid
              flx(l,i,j,k,lt2) = flx(l,i,j,k,1)*qgdnv(l,7)
           end do
#endif           
           ! Alpha
           do l = 1, ngrid
              flx(l,i,j,k,ndim+3) = qgdnv(l,4)
           end do
           ! Beta
           do l = 1, ngrid
              flx(l,i,j,k,ndim+4) = qgdnv(l,5)
           end do

           ! Total energy
           do l = 1, ngrid
              rhoetot(l)=qgdnv(l,4)*qgdnv(l,3)+qgdnv(l,5) + half*qgdnv(l,1)*qgdnv(l,2)**2
           end do
           do idim = 1,ndim-1
              do l = 1, ngrid
                 rhoetot(l)=rhoetot(l)+half*qgdnv(l,1)*qgdnv(l,idim+5)**2
              end do
           end do
           do l = 1, ngrid
              flx(l,i,j,k,ndim+2) = qgdnv(l,2)*(rhoetot(l)+qgdnv(l,3))
           end do

           ! Other advected quantities
           do n = ndim+5, nvar
              do l = 1, ngrid
                 flx(l,i,j,k,n) = flx(l,i,j,k,1)*qgdnv(l,n)
              end do
           end do

           ! Temporary Godunov states
           do l = 1, ngrid
              tmp(l,i,j,k,1) = qgdnv(l,2)   ! Normal velocity
           end do
           do l = 1,ngrid
              if(qgdnv(l,2)>zero)then       ! Internal energy flux
                 tmp(l,i,j,k,2) = (qleft (l,4)*qleft (l,3)+qleft (l,5))*qgdnv(l,2)
              else
                 tmp(l,i,j,k,2) = (qright(l,4)*qright(l,3)+qright(l,5))*qgdnv(l,2)
              end if
           end do

        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,c,gin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  integer ::i, j, k, l, n, idim
  real(dp)::eint, smalle
  real(dp),dimension(1:nvector),save::eken

  smalle = smallc**2/gamma/(gamma-one)

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do

           ! Debug
           if(debug)then
              do l = 1, ngrid
                 if(uin(l,i,j,k,1).le.smallr)then
                    write(*,*)'negative density'
                    write(*,*)uin(l,i,j,k,1)
                    stop
                 end if
              end do
           end if

           ! Compute velocities
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = uin(l,i,j,k,idim+1)/uin(l,i,j,k,1)
                end do
           end do

           ! Compute specific kinetic energy
           eken(1:ngrid) = zero
           do idim = 1, ndim
              do l = 1, ngrid
                 eken(l) = eken(l) + half*q(l,i,j,k,idim+1)**2
              end do
           end do

           ! Compute pressure
           do l = 1, ngrid
              eint = uin(l,i,j,k,ndim+2)-uin(l,i,j,k,1)*eken(l)
              q(l,i,j,k,ndim+2)=(eint-uin(l,i,j,k,ndim+4))/uin(l,i,j,k,ndim+3)
           end do

           ! Compute sound speed
           do l = 1, ngrid
              gamma = 1.0/uin(l,i,j,k,ndim+3)+1.0
              pinf = uin(l,i,j,k,ndim+4)/uin(l,i,j,k,ndim+3)/gamma
              c(l,i,j,k)=max(gamma*(q(l,i,j,k,ndim+2)+pinf)/q(l,i,j,k,1),smallc**2)
              c(l,i,j,k)=sqrt(c(l,i,j,k))
           end do

           ! Stiffen gas parameters
           do l = 1, ngrid
              q(l,i,j,k,ndim+3)=uin(l,i,j,k,ndim+3)
              q(l,i,j,k,ndim+4)=uin(l,i,j,k,ndim+4)
           end do

           ! Gravity predictor step
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) + gin(l,i,j,k,idim)*dt*half
              end do
           end do

        end do
     end do
  end do

  ! Passive scalar
  do n = ndim+5, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/uin(l,i,j,k,1)
              end do
           end do
        end do
     end do
  end do
 
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(q,dq,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi

!!$              ! ultrabee
!!$              if (ndim>0) then        
!!$                 do l = 1, ngrid
!!$                    dcen = q(l,i,j,k,2)*slope_type
!!$                    dlft = two/(zero+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
!!$                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
!!$                    dsgn = sign(one, dlft)
!!$                    slop = min(abs(dlft),abs(drgt))
!!$                    dlim = slop
!!$                    if((dlft*drgt)<=zero)dlim=zero
!!$                    dq(l,i,j,k,n,1) = dsgn*dlim
!!$                 end do
!!$              end if
!!$
!!$              ! superbee
!!$              if (ndim>0) then        
!!$                 do l = 1, ngrid
!!$                    dcen = q(l,i,j,k,2)*slope_type
!!$                    dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
!!$                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
!!$                    dsgn = sign(one, dlft)
!!$                    slop = min(abs(dlft),abs(drgt))
!!$                    dlim = slop
!!$                    if((dlft*drgt)<=zero)dlim=zero
!!$                    dq(l,i,j,k,n,1) = dsgn*dlim
!!$                 end do
!!$              end if

              if (ndim>0) then        
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              end if

#if NDIM>1              
              if (ndim>1) then
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end if
#endif
#if NDIM>2
              if (ndim>2) then
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              endif
#endif
           end do
        end do
     end do     
  end do
  
end subroutine uslope
