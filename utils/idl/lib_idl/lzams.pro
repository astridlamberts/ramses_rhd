function lzams, m

ALPHA    =0.39704170
BETA     =8.52762600
GAMMA    =0.00025546
DELTA    =5.43288900
EPSILON  =5.56357900
ZETA     =0.78866060
ETA      =0.00586685 
THETA    =1.71535900
IOTA     =6.59778800
KAPPA   =10.08855000
LAMBDA   =1.01249500
MU       =0.07490166
NU       =0.01077422
XI       =3.08223400
UPSILON =17.84778000
PI       =0.00022582

lsun = 3.90d33
msun = 2.0d33
msol = m/msun

return, lsun * (ALPHA*msol^5.5 + BETA*msol^11) / $
    (GAMMA+msol^3+DELTA*msol^5+EPSILON*msol^7+ $
     ZETA*msol^8+ETA*msol^9.5)

end

