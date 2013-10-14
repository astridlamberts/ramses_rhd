subroutine xinner_ana
  use hydro_parameters
  use variables
  implicit none
!
  integer :: j,k
!
  ! Periodic BC
  do k=ku1,ku2
     do j=ju1,ju2      
        uin(1,iu1  ,j,k,:) = uin(1,nx-2,j,k,:)
        uin(1,iu1+1,j,k,:) = uin(1,nx-1,j,k,:)
        uin(1,iu1+2,j,k,:) = uin(1,nx  ,j,k,:)
     end do
  end do
!
  return
end subroutine xinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine xouter_ana
  use hydro_parameters
  use variables
  implicit none

  integer ::k,j

  ! Periodic BC
  do k=ku1,ku2
     do j=ju1,ju2
        uin(1,iu2-2,j,k,:) = uin(1,1 ,j,k,:)
        uin(1,iu2-1,j,k,:) = uin(1,2 ,j,k,:)
        uin(1,iu2  ,j,k,:) = uin(1,3 ,j,k,:)
     end do
  end do
!
  return
end subroutine xouter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine yinner_ana
#if NDIM>1
  use shock_params
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i,ishift

  ! Periodic BC
  do k=ku1,ku2
     do i=iu1,iu2
        ishift=i-nshift
        if (ishift>0) then
           uin(1,i,ju1  ,k,:) = uin(1,ishift,ny-2,k,:)
           uin(1,i,ju1+1,k,:) = uin(1,ishift,ny-1,k,:)
           uin(1,i,ju1+2,k,:) = uin(1,ishift,ny  ,k,:)
        else
           uin(1,i,ju1  ,k,:) = uin(1,i,ju1+3,k,:)
           uin(1,i,ju1+1,k,:) = uin(1,i,ju1+3,k,:)
           uin(1,i,ju1+2,k,:) = uin(1,i,ju1+3,k,:)
        end if
     end do
  end do
!
#endif
  return
end subroutine yinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine youter_ana
#if NDIM>1
  use shock_params
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i,ishift

  do k=ku1,ku2
     do i=iu1,iu2
        ishift=i+nshift
        if (ishift<nx+1) then
           uin(1,i,ju2-2,k,:) = uin(1,ishift,1 ,k,:)
           uin(1,i,ju2-1,k,:) = uin(1,ishift,2 ,k,:)
           uin(1,i,ju2  ,k,:) = uin(1,ishift,3 ,k,:)
        else
           uin(1,i,ju2-2,k,:) = uin(1,i,ju2-3,k,:)
           uin(1,i,ju2-1,k,:) = uin(1,i,ju2-3,k,:)
           uin(1,i,ju2  ,k,:) = uin(1,i,ju2-3,k,:)
        end if
     end do
  end do
!
#endif
  return
end subroutine youter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zinner_ana
#if NDIM > 2
  use hydro_parameters
  use variables
  implicit none

  integer ::j,i

  ! Periodic BC
  do j=ju1,ju2
     do i=iu1,iu2
        uin(1,i,j,ku1  ,:) = uin(1,i,j,nz-2,:)
        uin(1,i,j,ku1+1,:) = uin(1,i,j,nz-1,:)
        uin(1,i,j,ku1+2,:) = uin(1,i,j,nz  ,:)
     end do
  end do
!
#endif 
  return
end subroutine zinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zouter_ana
#if NDIM > 2
  use hydro_parameters
  use variables
  implicit none

  integer ::j,i

  ! Periodic BC
  do j=ju1,ju2
     do i=iu1,iu2
        uin(1,i,j,ku2-2,:) = uin(1,i,j,1 ,:)
        uin(1,i,j,ku2-1,:) = uin(1,i,j,2 ,:)
        uin(1,i,j,ku2  ,:) = uin(1,i,j,3 ,:)
     end do
  end do
!
#endif
  return
end subroutine zouter_ana
