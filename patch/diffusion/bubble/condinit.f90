SUBROUTINE condinit(x,u,dx,nn)
  USE amr_parameters
  USE hydro_parameters
  USE poisson_parameters
  IMPLICIT NONE
  INTEGER ::nn                            ! Number of cells
  REAL(dp)::dx                            ! Cell size
  REAL(dp),DIMENSION(1:nvector,1:nvar)::u ! Conservative variables
  REAL(dp),DIMENSION(1:nvector,1:ndim)::x ! Cell center position.
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
  !==============================================================
  INTEGER::ivar,i,j,id,iu,iv,iw,ip
  REAL(dp)::x0,lambday,ky,lambdaz,kz
  REAL(dp)::rho0,p0,v0,v1,v2,r,y0,ie,y1,T1,T2,mu,yb,y1b,H
  REAL(dp),DIMENSION(1:nvector,1:nvar),SAVE::q   ! Primitive variables

  ! Call built-in initial condition generator
!  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........
  id=1; iu=2; iv=3; iw=4; ip=ndim+2 ; ie=ndim+3

  r=length_x(2)/2
  x0=x_center(2)
  y0=y_center(2)
  rho0=d_region(1)
  T1=p_region(1)/(d_region(1)) !=kT/mu*mh
  T2=t_region(2)/0.5
  yb=0.0
  y1b=0.0
  mu=0.5
  H=T1/abs(gravity_params(2))

  DO i=1,nn
     q(i,id)=rho0*EXP(-x(i,2)/H)
     q(i,ie)=t_region(1)
     q(i,2)=0.0
     q(i,3)=0.0
     q(i,ip)=q(i,id)*T1     

!!$     IF((x(i,1)-x0)**2+(x(i,2)-y0)**2 <= r**2) THEN
!!$        q(i,ie)=t_region(2)
!!$        q(i,id)=q(i,ip)/T2
!!$     END IF

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
  DO ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  END DO

END SUBROUTINE condinit
!================================================================
!================================================================
!================================================================
!================================================================
SUBROUTINE condinit2(x,u,dx,nn,ilevel)
  USE amr_commons
  USE hydro_parameters
  USE poisson_parameters
  IMPLICIT NONE

  INTEGER::ilevel

  INTEGER ::nn                            ! Number of cells
  REAL(dp)::dx                            ! Cell size
  REAL(dp),DIMENSION(1:nvector,1:nvar)::u ! Conservative variables
  REAL(dp),DIMENSION(1:nvector,1:ndim)::x ! Cell center position.
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
  INTEGER::ivar,i,j,id,iu,iv,iw,ip
  REAL(dp)::x0,lambday,ky,lambdaz,kz
  REAL(dp)::rho0,p0,v0,v1,v2,r,y0,ie,y1,T1,T2,mu,yb,y1b

  id=1; iu=2; iv=3; iw=4; ip=ndim+2 ; ie=ndim+3
  r=length_x(2)/2
  x0=x_center(2)
  y0=y_center(2)
  DO i=1,nn
     IF((x(i,1)-x0)**2+(x(i,2)-y0)**2 <= r**2) THEN
       u(i,id)=u(i,id) + dtnew(ilevel)*mass_input
       u(i,ip)=u(i,ip) + dtnew(ilevel)*mass_input*100.0/(gamma-1.0)
       u(i,ie)=u(i,ie) + dtnew(ilevel)*mass_input*t_region(2)
     END IF
  END DO

END SUBROUTINE condinit2
!#########################################################
!#########################################################
!#########################################################
!#########################################################
SUBROUTINE impose_iso(ilevel)
  USE amr_commons
  USE hydro_commons
  IMPLICIT NONE
  INTEGER::ilevel
  
  INTEGER::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  INTEGER::info,ibound,nx_loc,idim,ivar
  REAL(dp)::dx,dx_loc,scale,d,u,v,w,e
  REAL(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  REAL(dp),DIMENSION(1:8,1:3)::xc
  REAL(dp),DIMENSION(1:3)::skip_loc

  INTEGER ,DIMENSION(1:nvector),SAVE::ind_grid,ind_cell
  REAL(dp),DIMENSION(1:nvector,1:ndim),SAVE::xx
  REAL(dp),DIMENSION(1:nvector,1:nvar),SAVE::uu
 
  IF(verbose)WRITE(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  IF(ndim>0)skip_loc(1)=DBLE(icoarse_min)
  IF(ndim>1)skip_loc(2)=DBLE(jcoarse_min)
  IF(ndim>2)skip_loc(3)=DBLE(kcoarse_min)
  scale=DBLE(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  DO ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     IF(ndim>0)xc(ind,1)=(DBLE(ix)-0.5D0)*dx
     IF(ndim>1)xc(ind,2)=(DBLE(iy)-0.5D0)*dx
     IF(ndim>2)xc(ind,3)=(DBLE(iz)-0.5D0)*dx
  END DO
  
  ! Loop over all AMR linked lists (including physical boundaries).
  DO ibound=1,nboundary+ncpu

     IF(ibound==myid)THEN
        ncache=active(ilevel)%ngrid
     ELSE IF(ibound<=ncpu)THEN
        ncache=reception(ibound,ilevel)%ngrid
     ELSE
        ncache=boundary(ibound-ncpu,ilevel)%ngrid
     END IF

     ! Loop over grids by vector sweeps
     DO igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        IF(ibound==myid)THEN
           DO i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           END DO
        ELSE IF(ibound<=ncpu)THEN
           DO i=1,ngrid
              ind_grid(i)=reception(ibound,ilevel)%igrid(igrid+i-1)
           END DO
        ELSE
           DO i=1,ngrid
              ind_grid(i)=boundary(ibound-ncpu,ilevel)%igrid(igrid+i-1)
           END DO
        END IF
     
        ! Loop over cells
        DO ind=1,twotondim

           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           DO i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           END DO
           ! Gather cell centre positions
           DO idim=1,ndim
              DO i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              END DO
           END DO
           ! Rescale position from code units to user units
           DO idim=1,ndim
              DO i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
              END DO
           END DO
      
           do ivar=1,nvar
              DO i=1,ngrid
                 uu(i,ivar)=uold(ind_cell(i),ivar)
              END DO
           end do
           CALL condinit2(xx,uu,dx_loc,ngrid,ilevel)
           ! Scatter variables
           do ivar=1,nvar
              DO i=1,ngrid
                 uold(ind_cell(i),ivar)=uu(i,ivar)
              END DO
           enddo

        END DO
        ! End loop over cells

     END DO
     ! End loop over grids

  END DO
  ! End loop over linked lists

111 FORMAT('   Entering impose_iso for level ',I2)

END SUBROUTINE impose_iso
