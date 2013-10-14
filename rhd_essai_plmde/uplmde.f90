!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tracex(q,dq,c, qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use amr_commons, only :nstep
  implicit none

  integer::ngrid
  real(dp)::dx, dt  

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! Local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, ip
  real(dp)::dtdx,project_out
  real(dp)::cc, ccc, csq, r, u, p, a,b,v,w,eta
  real(dp)::drx, dux, dpx, dax, dvx, dwx
  real(dp)::alpham, alphap, alpha0x,alpha0y,alpha0z
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrxright, azryright, azrzright
  real(dp)::apleft,  amleft,  azrxleft, azryleft, azrzleft
  real(dp)::ni,h,lor,velsq,smallp, rp , rm,entho,vtot2,coef,vperp2
  real(dp)::D,Mx,My,E,rho,pre

  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=5
  project_out=one !zero
  ni = 1./(gamma-1.)
  smallp = 1.0e-10

!!!!
  entho = gamma/(gamma-one)
 

  
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid



              ! Cell centered values
              r   = q (l,i,j,k,1)
              u   = q (l,i,j,k,2)
              v   = q (l,i,j,k,3)
              w   = q (l,i,j,k,4)
              p   = q (l,i,j,k,5)
              h   = one +P/r*gamma/(gamma-1d0)
              csq = (gamma*p/(r*h) )
              cc  =sqrt(csq) ! relativistic sound speed


              eta =(one-u**2-csq*(v**2d0+w**2d0))**(half)
              b= h*csq
              vtot2= u**2d0+v**2d0+w**2d0
              lor=(one-vtot2)**(-half)
              a= cc*eta/(lor*r)






              ! TVD slopes in X direction
              drx = dq(l,i,j,k,1,1)
              dux = dq(l,i,j,k,2,1)
              dvx = dq(l,i,j,k,3,1)
              dwx = dq(l,i,j,k,4,1)
              dpx = dq(l,i,j,k,5,1)

!              write(*,*),drx,dux,dvx,dwx,dpx

              ! Supersonic fix for high-velocity gradients
              !!!!!!!!!!!!!!!regarder les trois directions?
              ccc = cc

              if(ABS(dux) > three*cc)ccc=zero
              if(ABS(dvx) > three*cc)ccc=zero 
              if(ABS(dwx) > three*cc)ccc=zero

!          if ((nstep .ge. 394) .and. (nstep .lt. 398)  .and. (l .eq. 161) .and. (i .eq.1) .and. (j .eq. 1) .and. (k .eq. 1))  then
 !            write(*,*),dux,dvx,three*cc,ccc,'gradients'
!endif



              ! Characteristic analysis along X direction
              alpham  =          (-half/a)*dux                       +half/b*dpx
              alpha0x = drx                                          -one/b*dpx
              alpha0y =          u*v/(one-u**two)*dux  +dvx         + v/(one-u**two)/lor**two/r/h*dpx
              alpha0z =          u*w/(one-u**two)*dux         + dwx  + w/(one-u**two)/lor**two/r/h*dpx 
              alphap  =          (half/a)*dux                        +half/b*dpx

              !eigenvalues
              coef= (one-vtot2)*(one-u**two-ccc**two*(v**two+w**two))
              spminus = (u*(one-ccc**two)-ccc*sqrt(coef))/(one-vtot2*ccc**two)*dtdx
              spplus  = (u*(one-ccc**two)+ccc*sqrt(coef))/(one-vtot2*ccc**two)*dtdx
              spzero  = (u    )*dtdx


              ! Right state


              if(((u*(one-ccc**2)+ccc*sqrt(coef))/(one-vtot2*ccc**two))>zero)spplus =-project_out
              if(((u*(one-ccc**2)-ccc*sqrt(coef))/(one-vtot2*ccc**two))>zero)spminus=-project_out              
              if( u     >zero)spzero =-project_out


              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrxright = half*(-one-spzero )*alpha0x
              azryright = half*(-one-spzero )*alpha0y
              azrzright = half*(-one-spzero )*alpha0z

!              rm=csq*spminus/(lor*(u-spminus))/r
!              rp=csq*spplus/(lor*(u-spplus))/r

              rm=u*(one-csq)-cc*eta/lor
              rp=u*(one-csq)+cc*eta/lor

              rm=cc*rm/(r*(u*cc+lor*eta))
              rp=cc*rp/(r*(u*cc-lor*eta))


              qp(l,i,j,k,1,1) = r + (apright+ azrxright+amright)
              qp(l,i,j,k,2,1) = u  + (-amright+apright)*a
              qp(l,i,j,k,5,1) = p + ( amright+apright)*b
              qp(l,i,j,k,3,1)= v + (amright*v*rm + azryright + apright*v*rp)
              qp(l,i,j,k,4,1)= w + (amright*w*rm + azrzright + apright*w*rp)

   !           if (qp(l,i,j,k,ir,1)<0.) write (*,*) 'd negative right',i,qp(l,i,j,k,ir,1) ,r
    !          if (qp(l,i,j,k,ip,1)<0.) write (*,*) 'p negative right',i,qp(l,i,j,k,ip,1),p
              qp(l,i,j,k,ir,1) = max(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = max(smallp, qp(l,i,j,k,ip,1))



              ! Left state

              spminus = (u*(one-ccc**two)-ccc*sqrt(coef))/(one-vtot2*ccc**two)*dtdx
              spplus  = (u*(one-ccc**two)+ccc*sqrt(coef))/(one-vtot2*ccc**two)*dtdx
              spzero  = (u    )*dtdx

                           
              if( (u*(one-ccc**two)+ccc*sqrt(coef))/(one-vtot2*ccc**two)   <=zero)spplus =+project_out
              if( (u*(one-ccc**two)-ccc*sqrt(coef))/(one-vtot2*ccc**two)   <=zero)spminus =+project_out
              if( u     <=zero)spzero =+project_out


              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrxleft  = half*(+one-spzero )*alpha0x
              azryleft  = half*(+one-spzero )*alpha0y
              azrzleft  = half*(+one-spzero )*alpha0z


              qm(l,i,j,k,1,1) = r +(apleft+azrxleft+amleft)
              qm(l,i,j,k,2,1) =  u + (-amleft+apleft)*a
              qm(l,i,j,k,5,1) =  p + ( amleft+apleft)*b
              qm(l,i,j,k,3,1)= v + (amleft*v*rm + azryleft + apleft*v*rp)
              qm(l,i,j,k,4,1)= w  + (w*rm*amleft + azrzleft + apleft*w*rp)

!              if (qm(l,i,j,k,ir,1)<0.) write (*,*) 'd negative left',i,qm(l,i,j,k,ir,1) ,r
!              if (qm(l,i,j,k,ip,1)<0.) write (*,*) 'p negative left',i,qm(l,i,j,k,ip,1), p
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = max(smallp, qm(l,i,j,k,ip,1))



              velsq=qp(l,i,j,k,2,1)*qp(l,i,j,k,2,1)+qp(l,i,j,k,3,1)*qp(l,i,j,k,3,1)+qp(l,i,j,k,4,1)*qp(l,i,j,k,4,1)
              if (velsq>1.d0) then
 !                write (*,*) '1. velsq>1.d0',i
                 qp(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,1)=q(l,i,j,k,ir:ip)
              endif


              velsq=qm(l,i,j,k,2,1)*qm(l,i,j,k,2,1)+qm(l,i,j,k,3,1)*qm(l,i,j,k,3,1)+qm(l,i,j,k,4,1)*qm(l,i,j,k,4,1)
              if (velsq>1.d0) then
  !               write (*,*) '1. velsq>1.d0',i
                 qp(l,i,j,k,:,1)=q(l,i,j,k,:)
                 qm(l,i,j,k,:,1)=q(l,i,j,k,:)
              endif
              


          if ((nstep .gt. 394) .and. (nstep .lt. 394)  .and. (l .eq. 161) .and. (i .eq.1) .and. (j .eq. 1) .and. (k .eq. 1))  then
             !ge 394 lt 398
                 !389 392
                 
                 r   = q (l,i,j,k,1)
                 u   = q (l,i,j,k,2)
                 v   = q (l,i,j,k,3)
                 w   = q (l,i,j,k,4)
                 p   = q (l,i,j,k,5)
                 h   = 1+P/r*gamma/(gamma-1)
                 vtot2= u**2+v**2+w**2
                 lor=(1.-vtot2)**(-1./2.)
                 D=r*lor
                 Mx=lor**2*r*h*u
                 My=lor**2*r*h*v
                 E=lor**2*r*h-p

!                 write(*,*),u,l,i,j,k
                 write(*,*),''
                 
                 write(*,*),r,u,v,p,l,sqrt(vtot2),ccc,'Qc'
                 write(*,*),qm(l,i,j,k,1,1),qm(l,i,j,k,2,1),qm(l,i,j,k,3,1),qm(l,i,j,k,5,1),'Ql'
                 write(*,*),qp(l,i,j,k,1,1),qp(l,i,j,k,2,1),qp(l,i,j,k,3,1),qp(l,i,j,k,5,1),'Qr'
                 write(*,*),D,Mx,My,E,'Uc'

                 rho=qm(l,i,j,k,1,1)
                 pre=qm(l,i,j,k,5,1)
                 h=1+gamma/(gamma-1)*pre/rho
                 lor=(1-(qm(l,i,j,k,2,1)**2+qm(l,i,j,k,3,1)**2+qm(l,i,j,k,4,1)**2))**(-1./2.)
                 D=rho*lor
                 Mx=lor**2*rho*h*qm(l,i,j,k,2,1)
                 My=lor**2*rho*h*qm(l,i,j,k,3,1)
                 E=lor**2*rho*h-pre
                 write(*,*),D,Mx,My,E,'Ul'

                 rho=qp(l,i,j,k,1,1)
                 pre=qp(l,i,j,k,5,1)
                 h=1+gamma/(gamma-1)*pre/rho
                 lor=(1-(qp(l,i,j,k,2,1)**2+qp(l,i,j,k,3,1)**2+qp(l,i,j,k,4,1)**2))**(-1./2.)
                 D=rho*lor
                 Mx=lor**2*rho*h*qp(l,i,j,k,2,1)
                 My=lor**2*rho*h*qp(l,i,j,k,3,1)
                 E=lor**2*rho*h-pre
                 write(*,*),D,Mx,My,E,'Ur'

!                 write(*,*),
                 write(*,*),''
 !                write(*,*),r*lor,lor**2*r*h*u,lor**2*r*h*v,l,i,j,k
              endif



!              if ((l >25) .and. (l<35))then
!                 write(*,*),qm(l,0,1,1,:,1),l,'qm'
!                 write(*,*),r,u,v,w,p,l
!                 write(*,*),qp(l,0,1,1,:,1),l,'qp'
!                 write(*,*)(qm(l,0,1,1,1,1)+qm(l,0,1,1,5,1)*gamma/(gamma-1))*(1-qm(l,0,1,1,2,1)**2-qm(l,0,1,1,3,1)**2)**(-1./2.)-qm(l,0,1,1,5,1),l,'Em'
!                 write(*,*)(qp(l,0,1,1,1,1)+qp(l,0,1,1,5,1)*gamma/(gamma-1))*(1-qp(l,0,1,1,2,1)**2-qp(l,0,1,1,3,1)**2)**(-1./2.)-qp(l,0,1,1,5,1),l,'Ep'
!              endif
                 !              write(*,*)(qm(l,i,j,k,1,1)+qm(l,i,j,k,5,1)*gamma/(gamma-1))*(1-qm(l,i,j,k,2,1)**2-qm(l,i,j,k,3,1)**2)**(-1./2.)-qm(l,i,j,k,5,1),l,'Em'
!              write(*,*)(qp(l,i,j,k,1,1)+qp(l,i,j,k,5,1)*gamma/(gamma-1))*(1-qp(l,i,j,k,2,1)**2-qp(l,i,j,k,3,1)**2)**(-1./2.)-qp(l,i,j,k,5,1),l,'Ep'
           end do
        end do
     end do
  end do



end subroutine tracex
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine tracexy(q,dq,c,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use amr_commons, only :nstep
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, IW, ip
  real(dp)::dtdx,dtdy,project_out
  real(dp)::cc, ccc, csq, r, u, v, p, w, h,  a
  real(dp)::alpham, alphap, alpha0x, alpha0y, alpha0z
  real(dp)::spminus, spplus, spzero
  real(dp)::dry, duy, dvy, dwy, dpy, day,say, DAX, SAX
  real(dp)::apright, amright, azrxright, azryright, azrzright, azaright
  real(dp)::apleft,  amleft,  azrxleft,  azryleft,  azrzleft,  azaleft
  real(dp)::eta,b,vtot,lor,vtot2,coef,rm,rp,smallp,velsq

  real(dp)::rho, pre,D,M,E,Mx,My


    
  dtdx = dt/dx; dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4;  ip=5
  project_out=one !zero


  

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid


              ! cell centered values
              r   = q (l,i,j,k,1)
              u   = q (l,i,j,k,2)
              v   = q (l,i,j,k,3)
              w   = q (l,i,j,k,4)
              p   = q (l,i,j,k,5)
              h   = 1+P/r*gamma/(gamma-1)
              csq = (gamma*p/(r*h) )
              cc  =sqrt(csq) ! relativistic sound speed

              eta =(1-v**2-csq*(u**2+w**2))**(1./2.)
              b= h*csq
              vtot2= u**2+v**2+w**2
              lor=(1.-vtot2)**(-1./2.)
              a= cc*eta/(lor*r)

              
              !Slopes in Y direction
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,k,k,iw,2)
              dpy = dq(l,i,j,k,ip,2)

              

              ! Supersonic fix for high-velocity gradients
              !!!!!!!!!!!!!!!regarder les trois directions?
              ccc = cc

              if(ABS(duy) > three*cc)ccc=zero
              if(ABS(dvy) > three*cc)ccc=zero 
              if(ABS(dwy) > three*cc)ccc=zero



              ! Characteristic analysis along Y direction
              alpham =             (-1./2./a)*dvy                        +1./2./b*dpy
              alpha0x =      duy +u*v/(1-v**2)*dvy         + u/(1-v**2)/lor**2/r/h*dpy
              alpha0y = dry                                                - 1./b*dpy
              alpha0z =          v*w/(1-v**2)*dvy  + dwy  + w/(1-v**2)/lor**2/r/h*dpy 
              alphap  =             (1./2./a)*dvy                        +1./2./b*dpy



              ! Top state
              !eigenvalues
              coef= (1-vtot2)*(1-v**2-ccc**2*(u**2+w**2))
              spminus = (v*(1-ccc**2)-ccc*sqrt(coef))/(1-vtot2*ccc**2)*dtdy
              spplus  = (v*(1-ccc**2)+ccc*sqrt(coef))/(1-vtot2*ccc**2)*dtdy
              spzero  = (v    )*dtdy


              if(((v*(1-ccc**2)+ccc*sqrt(coef))/(1-vtot2*ccc**2))>zero)spplus =-project_out
              if(((v*(1-ccc**2)-ccc*sqrt(coef))/(1-vtot2*ccc**2))>zero)spminus=-project_out              
              if( v     >zero)spzero =-project_out


              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham 
              azrxright = half*(-one-spzero )*alpha0x
              azryright = half*(-one-spzero )*alpha0y
              azrzright = half*(-one-spzero )*alpha0z

              rm=csq*spminus/(lor*(v-spminus))/r
              rp=csq*spplus/(lor*(v-spplus))/r

              qp(l,i,j,k,1,2) = r + (apright+ azryright+amright)
              qp(l,i,j,k,2,2) = u  + (amright*u*rm+ azrxright + apright*u*rp)
              qp(l,i,j,k,3,2) = v + (-amright + apright)*a
              qp(l,i,j,k,4,2) = w + (amright*w*rm + azrzright + apright*w*rp)
              qp(l,i,j,k,5,2) = p + ( amright+apright)*b
              
              qp(l,i,j,k,ir,2) = max(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = max(smallp, qp(l,i,j,k,ip,2))

              ! Bottom state

              !eigenvalues
              coef= (1-vtot2)*(1-v**2-ccc**2*(u**2+w**2))
              spminus = (v*(1-ccc**2)-ccc*sqrt(coef))/(1-vtot2*ccc**2)*dtdy
              spplus  = (v*(1-ccc**2)+ccc*sqrt(coef))/(1-vtot2*ccc**2)*dtdy
              spzero  = (v    )*dtdy


              if(((v*(1-ccc**2)+ccc*sqrt(coef))/(1-vtot2*ccc**2))<=zero)spplus =+project_out
              if(((v*(1-ccc**2)-ccc*sqrt(coef))/(1-vtot2*ccc**2))<=zero)spminus=+project_out              
              if( v     <zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrxleft  = half*(+one-spzero )*alpha0x
              azryleft  = half*(+one-spzero )*alpha0y
              azrzleft  = half*(+one-spzero )*alpha0z

              rm=csq*spminus/(lor*(v-spminus))/r
              rp=csq*spplus/(lor*(v-spplus))/r

              qm(l,i,j,k,1,2) = r + (apleft+ azryleft+amleft)
              qm(l,i,j,k,2,2) = u  + (amleft*u*rm+ azrxleft + apleft*u*rp)
              qm(l,i,j,k,3,2) = v + (-amleft + apleft)*a
              qm(l,i,j,k,4,2) = w + (amleft*w*rm + azrzleft + apleft*w*rp)
              qm(l,i,j,k,5,2) = p + ( amleft+apleft)*b

              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = max(smallp, qm(l,i,j,k,ip,2))


              velsq=qp(l,i,j,k,2,2)*qp(l,i,j,k,2,2)+qp(l,i,j,k,3,2)*qp(l,i,j,k,3,2)+qp(l,i,j,k,4,2)*qp(l,i,j,k,4,2)
              if (velsq>1.d0) then
 !                write (*,*) '1. velsq>1.d0',i
                 qp(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
                 qm(l,i,j,k,ir:ip,2)=q(l,i,j,k,ir:ip)
              endif


              velsq=qm(l,i,j,k,2,2)*qm(l,i,j,k,2,2)+qm(l,i,j,k,3,2)*qm(l,i,j,k,3,2)+qm(l,i,j,k,4,2)*qm(l,i,j,k,4,2)
              if (velsq>1.d0) then
  !               write (*,*) '1. velsq>1.d0',i
                 qp(l,i,j,k,:,2)=q(l,i,j,k,:)
                 qm(l,i,j,k,:,2)=q(l,i,j,k,:)
              endif

              rho=qm(l,i,j,k,1,1)
              pre=qm(l,i,j,k,5,1)
              h=1+gamma/(gamma-1)*pre/rho
              lor=(1-(qm(l,i,j,k,2,1)**2+qm(l,i,j,k,3,1)**2+qm(l,i,j,k,4,1)**2))**(-1./2.)
              D=rho*lor
              M=lor**2*rho*h*sqrt(qm(l,i,j,k,2,1)**2+qm(l,i,j,k,3,1)**2+qm(l,i,j,k,4,1)**2)
              E=lor**2*rho*h-pre
              if ((D<0).or. (E<0) .or. (M>E)) then
                 write(*,*),D,M,E,'qm1'
              endif


              rho=qm(l,i,j,k,1,2)
              pre=qm(l,i,j,k,5,2)
              h=1+gamma/(gamma-1)*pre/rho
              lor=(1-(qm(l,i,j,k,2,2)**2+qm(l,i,j,k,3,2)**2+qm(l,i,j,k,4,2)**2))**(-1./2.)
              D=rho*lor
              M=lor**2*rho*h*sqrt(qm(l,i,j,k,2,2)**2+qm(l,i,j,k,3,2)**2+qm(l,i,j,k,4,2)**2)
              E=lor**2*rho*h-pre
               if ((D<0).or. (E<0) .or. (M>E)) then
                 write(*,*),D,M,E,'qm2'
              endif


              rho=qp(l,i,j,k,1,1)
              pre=qp(l,i,j,k,5,2)
              h=1+gamma/(gamma-1)*pre/rho
              lor=(1-(qp(l,i,j,k,2,1)**2+qp(l,i,j,k,3,1)**2+qp(l,i,j,k,4,1)**2))**(-1./2.)
              D=rho*lor
              M=lor**2*rho*h*sqrt(qp(l,i,j,k,2,1)**2+qp(l,i,j,k,3,1)**2+qp(l,i,j,k,4,1)**2)
              E=lor**2*rho*h-pre
              if ((D<0).or. (E<0) .or. (M>E)) then
                 write(*,*),D,M,E,'qp1'
              endif


              rho=qp(l,i,j,k,1,2)
              pre=qp(l,i,j,k,5,2)
              h=1+gamma/(gamma-1)*pre/rho
              lor=(1-(qp(l,i,j,k,2,2)**2+qp(l,i,j,k,3,2)**2+qp(l,i,j,k,4,2)**2))**(-1./2.)
              D=rho*lor
              M=lor**2*rho*h*sqrt(qp(l,i,j,k,2,2)**2+qp(l,i,j,k,3,2)**2+qp(l,i,j,k,4,2)**2)
              E=lor**2*rho*h-pre
              if ((D<0).or. (E<0) .or. (M>E)) then
                 write(*,*),D,M,E,'qp2'
              endif

              if (u<0) then
!                 write(*,*)u
              endif
              if ((qp(l,i,j,k,2,1)*q(l,i,j,k,2)) <0) then
!                 write(*,*),'qp1 <0',qp(l,i,j,k,2,1),u,l,i,j,k
              endif

              if ((qp(l,i,j,k,2,2)*q(l,i,j,k,2)) <0) then
!                 write(*,*),'qp2 <0',qp(l,i,j,k,2,2),u,l,i,j,k
              endif

              if ((qm(l,i,j,k,2,2)*q(l,i,j,k,2)) <0) then
 !                write(*,*),'qm2 <0',qm(l,i,j,k,2,2),u,l,i,j,k
              endif

              if ((qm(l,i,j,k,2,1)*q(l,i,j,k,2)) <0) then
  !               write(*,*),'qm1 <0',qm(l,i,j,k,2,1),u,l,i,j,k
              endif

              if ((qm(l,i,j,k,1,2) .ne. r) .or. (qp(l,i,j,k,1,2) .ne. r)) then
!                 write(*,*),qm(l,i,j,k,1,2), r,qp(l,i,j,k,1,2), 'r souci'
                 endif
              if ((qm(l,i,j,k,2,2) .ne. u) .or. (qp(l,i,j,k,2,2) .ne. u)) then
!                 write(*,*),qm(l,i,j,k,2,2), u,qp(l,i,j,k,2,2), 'u pbe'
                 endif
              if ((qm(l,i,j,k,3,2) .ne. v) .or. (qp(l,i,j,k,3,2) .ne. v)) then
!                 write(*,*),qm(l,i,j,k,3,2), v,qp(l,i,j,k,3,2), 'v pbe'
                endif
              if ((qm(l,i,j,k,4,2) .ne. w) .or. (qp(l,i,j,k,4,2) .ne. w)) then
 !                write(*,*),qm(l,i,j,k,4,2), w,qp(l,i,j,k,4,2), 'w pbe'
                 endif
              if ((qm(l,i,j,k,5,2) .ne. p) .or. (qp(l,i,j,k,5,2) .ne. p)) then
 !                write(*,*),qm(l,i,j,k,5,2), p,qp(l,i,j,k,5,2), 'r pbe'
              endif

!             if ((nstep .ge. 389) .and. (nstep .lt. 392)) then

          if ((nstep .ge. 389) .and. (nstep .lt. 392)  .and. (l .eq. 161) .and. (i .eq.2) .and. (j .eq. 1) .and. (k .eq. 1))  then

                 !389 392
                 
                 r   = q (l,i,j,k,1)
                 u   = q (l,i,j,k,2)
                 v   = q (l,i,j,k,3)
                 w   = q (l,i,j,k,4)
                 p   = q (l,i,j,k,5)
                 h   = 1+P/r*gamma/(gamma-1)
                 vtot2= u**2+v**2+w**2
                 lor=(1.-vtot2)**(-1./2.)
                 D=r*lor
                 Mx=lor**2*r*h*u
                 My=lor**2*r*h*v
                 E=lor**2*r*h-p

!                 write(*,*),u,l,i,j,k
!                 write(*,*),''
                 
!                 write(*,*),r,u,v,p,l,'Qc'
!                 write(*,*),D,Mx,My,E,'Uc'
!                 write(*,*),
!                 write(*,*),''
 !                write(*,*),r*lor,lor**2*r*h*u,lor**2*r*h*v,l,i,j,k
              endif



              


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
                 a = q(l,i,j,k,n)     ! Cell centered values
                 u = q(l,i,j,k,iu)
                 v = q(l,i,j,k,iv)
                 dax = dq(l,i,j,k,n,1)    ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 sax = half*dtdy*(-v*day) ! Transverse
                 say = half*dtdx*(-u*dax) ! derivatives

                 ! Right state
                 spzero=(u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero=(u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax
                 
                 ! Top state
                 spzero=(v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero=(v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say
                 
              end do
           end do
        end do
     end do
  end do

end subroutine tracexy
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine tracexyz(q,dq,c,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c  

  ! declare local variables
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::i, j, k, l, n
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx,dtdy,dtdz,project_out
  real(dp)::cc, ccc, csq, r, u, v, w, p, a
  real(dp)::drx, dux, dvx, dwx, dpx, dax
  real(dp)::dry, duy, dvy, dwy, dpy, day
  real(dp)::drz, duz, dvz, dwz, dpz, daz
  real(dp)::alpham, alphap, alpha0r, alpha0u, alpha0v, alpha0w
  real(dp)::spminus, spplus, spzero
  real(dp)::apright, amright, azrright, azuright, azvright, azwright, azaright
  real(dp)::apleft,  amleft,  azrleft,  azuleft,  azvleft,  azwleft,  azaleft
  real(dp)::srx,sux,svx,swx,spx,sax
  real(dp)::sry,suy,svy,swy,spy,say
  real(dp)::srz,suz,svz,swz,spz,saz
    
  dtdx = dt/dx; dtdy = dt/dy; dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5
  project_out=one !zero

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
  
              ! Cell centered values
              cc  = c  (l,i,j,k)
              r   = q  (l,i,j,k,ir)
              u   = q  (l,i,j,k,iu)
              v   = q  (l,i,j,k,iv)
              w   = q  (l,i,j,k,iw)
              p   = q  (l,i,j,k,ip)
              csq = gamma*p/r

              ! TVD slopes in all 3 directions
              drx = dq(l,i,j,k,ir,1)
              dux = dq(l,i,j,k,iu,1)
              dvx = dq(l,i,j,k,iv,1)
              dwx = dq(l,i,j,k,iw,1)
              dpx = dq(l,i,j,k,ip,1)
              
              dry = dq(l,i,j,k,ir,2)
              duy = dq(l,i,j,k,iu,2)
              dvy = dq(l,i,j,k,iv,2)
              dwy = dq(l,i,j,k,iw,2)
              dpy = dq(l,i,j,k,ip,2)
              
              drz = dq(l,i,j,k,ir,3)
              duz = dq(l,i,j,k,iu,3)
              dvz = dq(l,i,j,k,iv,3)
              dwz = dq(l,i,j,k,iw,3)
              dpz = dq(l,i,j,k,ip,3)

              ! Transverse derivatives
              srx = half*dtdx*(-v*dry-w*drz - (dvy+dwz)*r      )
              spx = half*dtdx*(-v*dpy-w*dpz - (dvy+dwz)*gamma*p)
              sux = half*dtdx*(-v*duy-w*duz                    )
              svx = half*dtdx*(-v*dvy-w*dvz - (dpy)/r          )
              swx = half*dtdx*(-v*dwy-w*dwz - (dpz)/r          )

              sry = half*dtdx*(-u*drx-w*drz - (dux+dwz)*r      )
              spy = half*dtdx*(-u*dpx-w*dpz - (dux+dwz)*gamma*p)
              suy = half*dtdx*(-u*dux-w*duz - (dpx)/r          )
              svy = half*dtdx*(-u*dvx-w*dvz                    )
              swy = half*dtdx*(-u*dwx-w*dwz - (dpz)/r          )

              srz = half*dtdx*(-v*dry-u*drx - (dvy+dux)*r      )
              spz = half*dtdx*(-v*dpy-u*dpx - (dvy+dux)*gamma*p)
              suz = half*dtdx*(-v*duy-u*dux - (dpx)/r          )
              svz = half*dtdx*(-v*dvy-u*dvx - (dpy)/r          )
              swz = half*dtdx*(-v*dwy-u*dwx                    )

              ! Characteristic analysis along X direction
              alpham  = half*(dpx/csq - dux*r/cc)
              alphap  = half*(dpx/csq + dux*r/cc)
              alpha0r = drx - dpx/csq
              alpha0v = dvx
              alpha0w = dwx

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dux) > three*cc)ccc=zero

              ! Right state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)>zero)spplus =-project_out
              if((u-ccc)>zero)spminus=-project_out
              if( u     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap
              amright  = half*(-one-spminus)*alpham
              azrright = half*(-one-spzero )*alpha0r
              azvright = half*(-one-spzero )*alpha0v
              azwright = half*(-one-spzero )*alpha0w

              qp(l,i,j,k,ir,1) = r + (apright+amright+azrright)     +srx
              qp(l,i,j,k,iu,1) = u + (apright-amright         )*cc/r+sux
              qp(l,i,j,k,ip,1) = p + (apright+amright         )*csq +spx
              qp(l,i,j,k,iv,1) = v + (                azvright)     +svx
              qp(l,i,j,k,iw,1) = w + (                azwright)     +swx
              qp(l,i,j,k,ir,1) = max(smallr,qp(l,i,j,k,ir,1))

              ! Left state
              spminus = (u-ccc)*dtdx
              spplus  = (u+ccc)*dtdx
              spzero  = (u    )*dtdx
              if((u+ccc)<=zero)spplus =+project_out
              if((u-ccc)<=zero)spminus=+project_out
              if( u     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap
              amleft   = half*(+one-spminus)*alpham
              azrleft  = half*(+one-spzero )*alpha0r
              azvleft  = half*(+one-spzero )*alpha0v
              azwleft  = half*(+one-spzero )*alpha0w

              qm(l,i,j,k,ir,1) = r + (apleft+amleft+azrleft)     +srx
              qm(l,i,j,k,iu,1) = u + (apleft-amleft        )*cc/r+sux
              qm(l,i,j,k,ip,1) = p + (apleft+amleft        )*csq +spx
              qm(l,i,j,k,iv,1) = v + (              azvleft)     +svx
              qm(l,i,j,k,iw,1) = w + (              azwleft)     +swx
              qm(l,i,j,k,ir,1) = max(smallr, qm(l,i,j,k,ir,1))

              ! Characteristic analysis along Y direction
              alpham  = half*(dpy/csq - dvy*r/cc)
              alphap  = half*(dpy/csq + dvy*r/cc)
              alpha0r = dry - dpy/csq
              alpha0u = duy
              alpha0w = dwy

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dvy) > three*cc)ccc=zero

              ! Top state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)>zero)spplus =-project_out
              if((v-ccc)>zero)spminus=-project_out
              if( v     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azwright = half*(-one-spzero )*alpha0w

              qp(l,i,j,k,ir,2) = r + (apright+amright+azrright)     +sry
              qp(l,i,j,k,iv,2) = v + (apright-amright         )*cc/r+svy
              qp(l,i,j,k,ip,2) = p + (apright+amright         )*csq +spy
              qp(l,i,j,k,iu,2) = u + (                azuright)     +suy
              qp(l,i,j,k,iw,2) = w + (                azwright)     +swy
              qp(l,i,j,k,ir,2) = max(smallr,qp(l,i,j,k,ir,2))

              ! Bottom state
              spminus = (v-ccc)*dtdy
              spplus  = (v+ccc)*dtdy
              spzero  = (v    )*dtdy
              if((v+ccc)<=zero)spplus =+project_out
              if((v-ccc)<=zero)spminus=+project_out
              if( v     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azwleft  = half*(+one-spzero )*alpha0w

              qm(l,i,j,k,ir,2) = r + (apleft+amleft+azrleft)     +sry
              qm(l,i,j,k,iv,2) = v + (apleft-amleft        )*cc/r+svy
              qm(l,i,j,k,ip,2) = p + (apleft+amleft        )*csq +spy
              qm(l,i,j,k,iu,2) = u + (              azuleft)     +suy
              qm(l,i,j,k,iw,2) = w + (              azwleft)     +swy
              qm(l,i,j,k,ir,2) = max(smallr, qm(l,i,j,k,ir,2))

              ! Characteristic analysis along Z direction
              alpham  = half*(dpz/csq - dwz*r/cc)
              alphap  = half*(dpz/csq + dwz*r/cc)
              alpha0r = drz - dpz/csq
              alpha0u = duz
              alpha0v = dvz

              ! Supersonic fix for high-velocity gradients
              ccc = cc
              if(ABS(dwz) > three*cc)ccc=zero

              ! Front state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)>zero)spplus =-project_out
              if((w-ccc)>zero)spminus=-project_out
              if( w     >zero)spzero =-project_out

              apright  = half*(-one-spplus )*alphap 
              amright  = half*(-one-spminus)*alpham 
              azrright = half*(-one-spzero )*alpha0r
              azuright = half*(-one-spzero )*alpha0u
              azvright = half*(-one-spzero )*alpha0v

              qp(l,i,j,k,ir,3) = r + (apright+amright+azrright)     +srz
              qp(l,i,j,k,iw,3) = w + (apright-amright         )*cc/r+swz
              qp(l,i,j,k,ip,3) = p + (apright+amright         )*csq +spz
              qp(l,i,j,k,iu,3) = u + (                azuright)     +suz
              qp(l,i,j,k,iv,3) = v + (                azvright)     +svz
              qp(l,i,j,k,ir,3) = max(smallr,qp(l,i,j,k,ir,3))

              ! Back state
              spminus = (w-ccc)*dtdz
              spplus  = (w+ccc)*dtdz
              spzero  = (w    )*dtdz
              if((w+ccc)<=zero)spplus =+project_out
              if((w-ccc)<=zero)spminus=+project_out
              if( w     <=zero)spzero =+project_out

              apleft   = half*(+one-spplus )*alphap 
              amleft   = half*(+one-spminus)*alpham 
              azrleft  = half*(+one-spzero )*alpha0r
              azuleft  = half*(+one-spzero )*alpha0u
              azvleft  = half*(+one-spzero )*alpha0v

              qm(l,i,j,k,ir,3) = r + (apleft+amleft+azrleft)     +srz
              qm(l,i,j,k,iw,3) = w + (apleft-amleft        )*cc/r+swz
              qm(l,i,j,k,ip,3) = p + (apleft+amleft        )*csq +spz
              qm(l,i,j,k,iu,3) = u + (              azuleft)     +suz
              qm(l,i,j,k,iv,3) = v + (              azvleft)     +svz
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
                 a   =  q(l,i,j,k,n)    ! Cell centered values
                 u   =  q(l,i,j,k,iu)
                 v   =  q(l,i,j,k,iv)
                 w   =  q(l,i,j,k,iw)
                 dax = dq(l,i,j,k,n,1)  ! TVD slopes
                 day = dq(l,i,j,k,n,2)
                 daz = dq(l,i,j,k,n,3)
                 sax = half*dtdx*(-v*day-w*daz) ! Transverse
                 say = half*dtdx*(-u*dax-w*daz) ! derivatives
                 saz = half*dtdx*(-v*day-u*dax) ! 

                 
                 ! Right state
                 spzero = (u    )*dtdx
                 if(u>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*dax
                 qp(l,i,j,k,n,1) = a + azaright + sax
                 
                 ! Left state
                 spzero = (u    )*dtdx
                 if(u<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*dax
                 qm(l,i,j,k,n,1) = a + azaleft + sax

                 ! Top state
                 spzero = (v    )*dtdy
                 if(v>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*day
                 qp(l,i,j,k,n,2) = a + azaright + say
                 
                 ! Bottom state
                 spzero = (v    )*dtdy
                 if(v<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*day
                 qm(l,i,j,k,n,2) = a + azaleft + say

                 ! Front state
                 spzero = (w    )*dtdy
                 if(w>zero)spzero=-project_out
                 azaright = half*(-one-spzero )*daz
                 qp(l,i,j,k,n,3) = a + azaright + saz
                 
                 ! Back state
                 spzero = (w    )*dtdy
                 if(w<=zero)spzero=+project_out
                 azaleft = half*(+one-spzero )*daz
                 qm(l,i,j,k,n,3) = a + azaleft + saz
              end do
           end do
        end do
     end do
  end do

end subroutine tracexyz
#endif
