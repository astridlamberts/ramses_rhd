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


real(dp)::r_0!=1.0d0!2 rmask
real(dp)::H_0!=0.1d0 !0.1 r_0
real(dp)::cs_0!=0.1d0
real(dp)::p_0!=0.01d0!cs0**2*rho_0
real(dp)::P_amb!=1.0e-5!*p_0
real(dp)::rho_amb!=1.0e-7!*rho_0

real(dp)::vinf1=422.0d0
real(dp)::Mach1=1.0d0
real(dp)::Mdot1=10.0d0
real(dp)::rstar1=0.046d0
real(dp)::rmask1=0.06d0
real(dp)::a=1.0d0
real(dp)::beta1=1.0d0
logical::disk_presence=.false.

end module wind_parameters
