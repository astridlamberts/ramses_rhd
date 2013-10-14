module wind_parameters
!use amr_commons
!use amr_parameters


implicit none


real (kind=8):: Mdot1=1.0d0
real (kind=8):: vinf1=422.0d0
real (kind=8):: rstar1=0.046d0
real (kind=8):: r01=0.092d0
real (kind=8):: mach1=30
real (kind=8):: beta1=1.0d0
real (kind=8):: M1=15.0d0
character(len=4)    :: type1='pulsar'
integer      :: f1=0

real (kind=8):: Mdot2=1.0d0
real (kind=8):: vinf2=422.0d0
real (kind=8):: rstar2=0.046d0
real (kind=8):: r02=0.092d0
real (kind=8):: mach2=30
real (kind=8):: beta2=1.0d0
real (kind=8):: M2=15.0d0
character(len=4)  :: type2='pulsar'
integer      :: f2=0

real (kind=8)::R_rho=1000.0d0
real (kind=8)::R_P=1.0d0 
character (len=1) ::direction='d'
real (kind=8)::a=1.0d0
real (kind=8)::xcen=0.5d0
real (kind=8)::ycen=0.5d0
real (kind=8)::zcen=0.5d0
logical      ::rotation=.false.

real(kind=8),allocatable,dimension(:)::px,py,pz

real(kind=8)::x01,y01,z01,x02,y02,z02,vx01,vy01,vz01,vx02,vy02,vz02
real(kind=8)::peri=0.0d0
real(kind=8)::inc=0.0d0
real(kind=8)::node=0.0d0
real(kind=8)::theta=0.0d0
real(kind=8)::exc=0.0d0
real(kind=8)::Gstar=39.01d0

real(kind=8)::pe=2.0d0
real(kind=8)::qu=1.0d0
real(kind=8)::rho_0=1.7e6
logical::disk_presence=.false.
real(kind=8)::rin=0.092
real(kind=8)::HR=0.1d0


real(kind=8)::H_0!=0.1d0 !0.1 r_0
real(kind=8)::cs_0!=0.1d0
real(kind=8)::p_0!=0.01d0!cs0**2*rho_0
!real(kind=8)::P_amb!=1.0e-5!*p_0
!real(kind=8)::rho_amb!=1.0e-7!*rho_0











end module wind_parameters
