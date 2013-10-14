program euler_3d
implicit none


real::x01,y01,z01,vx01,vy01,vz01,xini1,xfin1,yini1,yfin1,zini1,zfin1,vxini1,vxfin1,vyini1,vyfin1,vzini1,vzfin1,ax1,ay1,az1,xprime
real::x02,y02,z02,vx02,vy02,vz02,xini2,xfin2,yini2,yfin2,zini2,zfin2,vxini2,vxfin2,vyini2,vyfin2,vzini2,vzfin2,ax2,ay2,az2
real::dt=0.000001d0
real::r
real::a=1.0d0
real::a1=1.0d0
integer::i
real::G=39.01d0
real::M1=38.0D0
real::M2=12.0d0
real::e=0.8d0
real::inc,noeud,peri,theta,dthetadt
real::c,costheta,p,h
real::x1,y1,r1,r2
!center of mass
real::xc=0.0d0
real::yc=0.0d0 
real::xf1,yf1,zf1,xf2,yf2,zf2,vxf1,vyf1,vzf1,vxf2,vyf2,vzf2,rf,ax1f,ay1f,az1f,ax2f,ay2f,az2f
open(20,file='test.txt')

G=G/a




!attention, mettre des radians
inc=30.d0/180.d0*3.1415d0
peri=60.d0/180.d0*3.1415d0
noeud=0.d0
theta=90.d0/180.d0*3.1415d0

xini1=M2/(M1+M2)*a*(1-e)*(cos(noeud)*cos(peri+theta)-sin(noeud)*sin(peri+theta)*cos(inc))
yini1=M2/(M1+M2)*a*(1-e)*(sin(noeud)*cos(peri+theta)+cos(noeud)*sin(peri+theta)*cos(inc))
zini1=M2/(M1+M2)*a*(1-e)*sin(peri+theta)*sin(inc)

xini2=xini1-a*(1-e)*(cos(noeud)*cos(peri+theta)-sin(noeud)*sin(peri+theta)*cos(inc))
yini2=yini1-a*(1-e)*(sin(noeud)*cos(peri+theta)+cos(noeud)*sin(peri+theta)*cos(inc))
zini2=zini1-a*(1-e)*sin(peri+theta)*sin(inc)

r=sqrt((xini1-xini2)**2+(yini1-yini2)**2+(zini1-zini2)**2)
write(*,*),r,xini1,xini2,yini1,yini2,zini1,zini2



	ax1=-G*M2/r**2*(xini1-xini2)/r
	ay1=-G*M2/r**2*(yini1-yini2)/r
        az1=-G*M2/r**2*(zini1-zini2)/r

	ax2=-G*M1/r**2*(xini2-xini1)/r
	ay2=-G*M1/r**2*(yini2-yini1)/r
        az2=-G*M1/r**2*(zini2-zini1)/r 




dthetadt=sqrt(a*(1-e**2)*G*(M1+M2)/(r**4))


r1=sqrt(xini1**2+yini1**2+zini1**2)
r2=sqrt(xini2**2+yini2**2+zini2**2)

vxini1=dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r1
vyini1=dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r1
vzini1=dthetadt*cos(peri+theta)*sin(inc)*r1

vxini2=-dthetadt*(-cos(noeud)*sin(peri+theta)-sin(noeud)*cos(peri+theta)*cos(inc))*r2
vyini2=-dthetadt*(-sin(noeud)*sin(peri+theta)+cos(noeud)*cos(peri+theta)*cos(inc))*r2
vzini2=-dthetadt*cos(peri+theta)*sin(inc)*r2





write(*,*),vxini1,vyini1,vzini1,vxini2,vyini2,vzini2


do i=1,800000
   r=sqrt((xini1-xini2)**2+(yini1-yini2)**2+(zini1-zini2)**2)


   ax1=-G*m2*(xini1-xini2)/r**3.0d0                                                                                                                                                                               
   ay1=-G*m2*(yini1-yini2)/r**3.0d0   
   az1=-G*m2*(zini1-zini2)/r**3.0d0                                                                                                                                                                               

   ax2=-G*m1*(xini2-xini1)/r**3.0d0                                                                                                                                                                               
   ay2=-G*m1*(yini2-yini1)/r**3.0d0                                                                                                                                                                               
   az2=-G*m1*(zini2-zini1)/r**3.0d0                                                                                                                                                                               


   xf1=xini1+vxini1*dt+ax1*dt**2.0d0/2.d0
   yf1=yini1+vyini1*dt+ay1*dt**2.0d0/2.d0
   zf1=zini1+vzini1*dt+az1*dt**2.0d0/2.d0

   xf2=xini2+vxini2*dt+ax2*dt**2.0d0/2.d0
   yf2=yini2+vyini2*dt+ay2*dt**2.0d0/2.d0
   zf2=zini2+vzini2*dt+az2*dt**2.0d0/2.d0

   rf=sqrt((xf1-xf2)**2+(yf1-yf2)**2+(zf1-zf2)**2)

   ax1f=-G*m2*(xf1-xf2)/rf**3.0d0  
   ay1f=-G*m2*(yf1-yf2)/rf**3.0d0  
   az1f=-G*m2*(zf1-zf2)/rf**3.0d0      
   ax2f=-G*m1*(xf2-xf1)/rf**3.0d0      
   ay2f=-G*m1*(yf2-yf1)/rf**3.0d0      
   az2f=-G*m1*(zf2-zf1)/rf**3.0d0          

   vxf1=vxini1+(ax1+ax1f)/2.d0*dt
   vyf1=vyini1+(ay1+ay1f)/2.d0*dt
   vzf1=vzini1+(az1+az1f)/2.d0*dt

   vxf2=vxini2+(ax2+ax2f)/2.d0*dt
   vyf2=vyini2+(ay2+ay2f)/2.d0*dt
   vzf2=vzini2+(az2+az2f)/2.d0*dt

   xini2=xf2
   yini2=yf2
   zini2=zf2

   vxini2=vxf2
   vyini2=vyf2
   vzini2=vzf2



   xini1=xf1
   yini1=yf1
   zini1=zf1

   vxini1=vxf1
   vyini1=vyf1
   vzini1=vzf1

   write(20,*),xini1,yini1,zini1,i,xini2,yini2,zini2!,vxini1,vyini1,vxini2,vyini2!,r,p/(1.0d0+e*xini/r),i,costheta,xprime

 enddo  
 close(20)
end program euler_3D
