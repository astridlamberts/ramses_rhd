module wind_params
!use amr_commons
!use amr_parameters


implicit none


real (kind=8):: Mdot1=1.0d0
real (kind=8):: vinf1=422.0d0
real (kind=8):: rstar1=0.046d0
real (kind=8):: r01=0.092d0
real (kind=8):: mach1=30
real (kind=8):: beta1=1.0d0
character(len=4)    :: type1='star'
integer      :: f1=1

real (kind=8):: Mdot2=1.0d0
real (kind=8):: vinf2=422.0d0
real (kind=8):: rstar2=0.046d0
real (kind=8):: r02=0.092d0
real (kind=8):: mach2=30
real (kind=8):: beta2=1.0d0
character(len=4)  :: type2='pulsar'
integer      :: f2=0

real (kind=8)::a=1.0d0
real (kind=8)::period=0.0001d0
real (kind=8)::R_rho=1000.0d0
real (kind=8)::R_P=1.0d0 
real (kind=8)::x01=0.5d0
real (kind=8)::y01=0.5d0
real (kind=8)::z01=0.5d0
logical      ::rotation=.false.

real(kind=8),allocatable,dimension(:)::px,py,pz

end module wind_params
