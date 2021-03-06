module wind_parameters
use amr_parameters

implicit none

real(dp)::GM=484
real(dp)::x01=0.5d0
real(dp)::y01=0.5d0
real(dp)::z01=0.5d0
real(dp)::pe=2.0d0
real(dp)::qu=1.0d0
real(dp)::rho_0=1.7e6
logical::disk_presence=.false.
real(dp)::r_0=0.092
real(dp)::HR=0.1d0


real(dp)::H_0!=0.1d0 !0.1 r_0
real(dp)::cs_0!=0.1d0
real(dp)::p_0!=0.01d0!cs0**2*rho_0
real(dp)::P_amb!=1.0e-5!*p_0
real(dp)::rho_amb!=1.0e-7!*rho_0

real (kind=8):: Mdot2=1.0d0
real (kind=8):: vinf2=422.0d0
real (kind=8):: rstar2=0.046d0
real (kind=8):: rmask2=0.092d0
real (kind=8):: mach2=30
real (kind=8):: beta2=1.0d0
character(len=4)  :: type2='puls'
integer      :: f2=0


real (kind=8)::a=1.0d0
logical      ::rotation=.false.
real(kind=8)::x02,y02,z02,vx02,vy02
real(kind=8)::exc=0.337d0

!real(kind=8),allocatable,dimension(:)::px,py,pz



real(dp)::vinf1=422.0d0
real(dp)::Mach1=1.0d0
real(dp)::Mdot1=10.0d0
real(dp)::rstar1=0.046d0
real(dp)::rmask1=0.06d0
real(dp)::beta1=1.0d0

end module wind_parameters
