import matplotlib.pyplot as plt
import numpy as np

def f_mignone(D,M,Eprim,f):
    size=10000
    R=np.array(size)
    u2=M**2/((R+D)**2-M**2) 
    lor=(1+u2)**(0.5)
    Xsi=(R-u2/(lor+1)*D)/lor**2
    rho=D/lor
    P=(2*xsi*(xsi+2*rho))/(5*(xsi+rho)+sqrt(9*xsi**2+18*rho*xsi+25*rho**2))
    f  = R-P-Eprim


def f_mignone_prim(D,M,Eprim,f):
    u2=M**2/((R+D)**2-M**2) 
    lor=(1+u2)**(0.5)
    Xsi=(R-u2/(lor+1)*D)/lor**2
    rho=D/lor
    P=(2*xsi*(xsi+2*rho))/(5*(xsi+rho)+sqrt(9*xsi**2+18*rho*xsi+25*rho**2))

    dpdxsi=(2*xsi+2*xsi-5*P)/(5*rho+5*xsi-8*P)
    dpdrho=(2*xsi-5*P)/(5*rho+5*xsi-8*P)
    dv2dR=-2*M**2/(R+D)**3
    
    dxsidR=1/lor**2-lor/2*(D+2*lor*xsi)*dv2dR
    drhodR=D*lor/2*dv2dR
    
    dpdR=dpdxsi*dxsidR+dpdrho*drhodR
    
    f_prim=1-dpdR

