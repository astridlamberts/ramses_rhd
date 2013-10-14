!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  
  real(dp)::dtcell,smallp
  integer::k,idim
  
  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
  end do
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/uu(k,1)
     end do
  end do
  do idim = 1,ndim
     do k = 1, ncell
        uu(k,ndim+2) = uu(k,ndim+2)-half*uu(k,1)*uu(k,idim+1)**2
     end do
  end do

  ! Compute pressure
  do k = 1, ncell
     uu(k,ndim+2) = (uu(k,ndim+2)-uu(k,ndim+4))/uu(k,ndim+3)
!     uu(k,ndim+2) = max(uu(k,ndim+2),uu(k,1)*smallp)
  end do

  ! Debug
  if(debug)then
!     write(*,*)'OK'
     do k = 1, ncell 
        gamma=1.0/uu(k,ndim+3)+1.0
        pinf=uu(k,ndim+4)/uu(k,ndim+3)/gamma
!        if(uu(k,1).le.smallr.or.gamma*(uu(k,ndim+2)+pinf)/uu(k,1).le.0.)then
        if(uu(k,1).le.smallr)then
           write(*,*)'stop in cmpdt'
           write(*,*)'dx   =',dx
           write(*,*)'ncell=',ncell
           write(*,*)'rho  =',uu(k,1)
           write(*,*)'P    =',uu(k,ndim+2)
           write(*,*)'pinf =',pinf
           write(*,*)'gamma=',gamma
           write(*,*)'vel  =',uu(k,2:ndim+1)
           stop
        end if
     end do
  end if

  ! Compute sound speed
  do k = 1, ncell 
     gamma=1.0/uu(k,ndim+3)+1.0
     pinf=uu(k,ndim+4)/uu(k,ndim+3)/gamma
     uu(k,ndim+2)=sqrt(max(gamma*(uu(k,ndim+2)+pinf)/uu(k,1),smallc**2))
  end do

  ! Compute wave speed
  do k = 1, ncell
     uu(k,ndim+2)=dble(ndim)*uu(k,ndim+2)
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        uu(k,ndim+2)=uu(k,ndim+2)+abs(uu(k,idim+1))
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1, ncell
     uu(k,1)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        uu(k,1)=uu(k,1)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     uu(k,1)=uu(k,1)*dx/uu(k,ndim+2)**2
     uu(k,1)=MAX(uu(k,1),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell=dx/uu(k,ndim+2)*(sqrt(one+two*courant_factor*uu(k,1))-one)/uu(k,1)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,nn)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer nn
  real(dp)::ug(1:nvector,1:nvar)
  real(dp)::um(1:nvector,1:nvar)
  real(dp)::ud(1:nvector,1:nvar)
  logical ::ok(1:nvector)
  
  integer::k,idim
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error
  real(dp)::gammag,gammam,gammad,pinfg,pinfm,pinfd

  ! Convert to primitive variables
  do k = 1,nn
     ug(k,1) = max(ug(k,1),smallr)
     um(k,1) = max(um(k,1),smallr)
     ud(k,1) = max(ud(k,1),smallr)
  end do
  do idim = 1,ndim
     do k = 1,nn
        ug(k,idim+1) = ug(k,idim+1)/ug(k,1)
        um(k,idim+1) = um(k,idim+1)/um(k,1)
        ud(k,idim+1) = ud(k,idim+1)/ud(k,1)
     end do
  end do
  do k = 1,nn
     eking(k) = zero
     ekinm(k) = zero
     ekind(k) = zero
  end do
  do idim = 1,ndim
     do k = 1,nn
        eking(k) = eking(k) + half*ug(k,1)*ug(k,idim+1)**2
        ekinm(k) = ekinm(k) + half*um(k,1)*um(k,idim+1)**2
        ekind(k) = ekind(k) + half*ud(k,1)*ud(k,idim+1)**2
     end do
  end do
  do k = 1,nn
     ug(k,ndim+2) = (ug(k,ndim+2)-eking(k)-ug(k,ndim+4))/ug(k,ndim+3)
     um(k,ndim+2) = (um(k,ndim+2)-ekinm(k)-um(k,ndim+4))/um(k,ndim+3)
     ud(k,ndim+2) = (ud(k,ndim+2)-ekind(k)-ud(k,ndim+4))/ud(k,ndim+3)
  end do  

  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,nn
        dg=ug(k,1); dm=um(k,1); dd=ud(k,1)
        error=2.0d0*MAX( &
             & ABS((dd-dm)/(dd+dm+floor_d)) , &
             & ABS((dm-dg)/(dm+dg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
     do k=1,nn
     end do
  end if

  if(err_grad_p >= 0.)then
     do k=1,nn
        pg=ug(k,ndim+2); pm=um(k,ndim+2); pd=ud(k,ndim+2)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_p)), &
             & ABS((pm-pg)/(pm+pg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if

  if(err_grad_a >= 0.)then
     do k=1,nn
        pg=ug(k,ndim+3); pm=um(k,ndim+3); pd=ud(k,ndim+3)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm)), &
             & ABS((pm-pg)/(pm+pg)) )
        ok(k) = ok(k) .or. error > err_grad_a
     end do
  end if

  if(err_grad_b >= 0.)then
     do k=1,nn
        pg=ug(k,ndim+4); pm=um(k,ndim+4); pd=ud(k,ndim+4)
        error=2.0d0*MAX( &
             & ABS((pd-pm)/(pd+pm+floor_b)), &
             & ABS((pm-pg)/(pm+pg+floor_b)) )
        ok(k) = ok(k) .or. error > err_grad_b
     end do
  end if

  if(err_grad_u >= 0.)then
     do idim = 1,ndim
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
           gammag=1.0/ug(k,ndim+3)+1.0; pinfg=(gammag-1.0)/gammag*ug(k,ndim+4)
           gammam=1.0/um(k,ndim+3)+1.0; pinfm=(gammam-1.0)/gammam*um(k,ndim+4)
           gammad=1.0/ud(k,ndim+3)+1.0; pinfd=(gammad-1.0)/gammad*ud(k,ndim+4)           
           cg=sqrt(max(gammag*(ug(k,ndim+2)+pinfg)/ug(k,1),floor_u**2))
           cm=sqrt(max(gammam*(um(k,ndim+2)+pinfm)/um(k,1),floor_u**2))
           cd=sqrt(max(gammad*(ud(k,ndim+2)+pinfd)/ud(k,1),floor_u**2))
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_acoustic(qleft,qright,qgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv

  ! local variables
  integer::i,n
  real(dp)::smallp

  ! local arrays
  real(dp),dimension(1:nvector),save::rl   ,ul   ,pl   ,cl, gammal, pinfl
  real(dp),dimension(1:nvector),save::rr   ,ur   ,pr   ,cr, gammar, pinfr 
  real(dp),dimension(1:nvector),save::ro   ,uo   ,po   ,co, gammao, pinfo
  real(dp),dimension(1:nvector),save::rstar,ustar,pstar,cstar
  real(dp),dimension(1:nvector),save::wl   ,wr   ,wo   
  real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:nvector),save::frac

  if(debug)write(*,*)'enter riemann'

  ! Initial states pressure, density and velocity
  do i=1,ngrid
     rl    (i)=qleft (i,1)
     ul    (i)=qleft (i,2)
     pl    (i)=qleft (i,3)
     gammal(i)=one/qleft(i,4)+one
     pinfl (i)=qleft(i,5)/qleft(i,4)/gammal(i)
     rr    (i)=qright(i,1)
     ur    (i)=qright(i,2)
     pr    (i)=qright(i,3)
     gammar(i)=one/qright(i,4)+one
     pinfr (i)=qright(i,5)/qright(i,4)/gammar(i)
  end do

  ! Acoustic star state
  do i=1,ngrid
     cl(i) = sqrt(max(gammal(i)*(pl(i)+pinfl(i))/rl(i),smallc**2))
     cr(i) = sqrt(max(gammar(i)*(pr(i)+pinfr(i))/rr(i),smallc**2))
     wl(i) = cl(i)*rl(i)
     wr(i) = cr(i)*rr(i)
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i))) / (wl(i)+wr(i))
     ustar(i) = ((wr(i)*ur(i)+wl(i)*ul(i))+(pl(i)-pr(i))) / (wl(i)+wr(i))
!!$     pstar(i) = MAX(pstar(i),zero)
  end do

  ! Left going or right going contact wave
  do i=1,ngrid   
     sgnm(i) = sign(one,ustar(i))
  end do

  ! Left or right unperturbed state
  do i=1,ngrid
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
        co(i) = cl(i)
        gammao(i) = gammal(i)
        pinfo(i) = pinfl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
        co(i) = cr(i)
        gammao(i) = gammar(i)
        pinfo(i) = pinfr(i)
     end if
  end do

  ! Star region density and sound speed
  do i=1,ngrid
     rstar(i) = ro(i)+(pstar(i)-po(i))/co(i)**2
     rstar(i) = max(rstar(i),smallr)
     cstar(i) = sqrt(max(gammao(i)*(pstar(i)+pinfo(i))/rstar(i),smallc**2))
     cstar(i) = max(cstar(i),smallc)
  end do

  ! Head and tail speed of rarefaction
  do i=1,ngrid
     spout(i) = co   (i)-sgnm(i)*uo   (i)
     spin (i) = cstar(i)-sgnm(i)*ustar(i)
  end do

  ! Shock speed
  do i=1,ngrid
     ushock(i) = half*(spin(i)+spout(i))
!!$     ushock(i) = (rstar(i)*cstar(i) + ro(i)*co(i))/(2*rstar(i)) &
!!$          &    - sgnm(i)*ustar(i)
  end do
  do i=1,ngrid
     if(pstar(i)>=po(i))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
  end do

  ! Sample the solution at x/t=0
  do i=1,ngrid
     if(spout(i)<=zero)then      ! Initial state
        qgdnv(i,1) = ro(i)
        qgdnv(i,2) = uo(i)
        qgdnv(i,3) = po(i)
     else if(spin(i)>=zero)then  ! Star region
        qgdnv(i,1) = rstar(i)
        qgdnv(i,2) = ustar(i)
        qgdnv(i,3) = pstar(i)
     else                        ! Rarefaction
        frac(i) = spout(i)/(spout(i)-spin(i))
        qgdnv(i,1) = frac(i)*rstar(i) + (one - frac(i))*ro(i)
        qgdnv(i,2) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
        qgdnv(i,3) = frac(i)*pstar(i) + (one - frac(i))*po(i)
     end if
  end do

  ! Passive scalars
  do n = 4,nvar
     do i=1,ngrid
        if(sgnm(i)==one)then
           qgdnv(i,n) = qleft (i,n)
        else
           qgdnv(i,n) = qright(i,n)
        end if
     end do
  end do

end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_approx(qleft,qright,qgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:nvar)::qleft,qright,qgdnv

  real(dp)::rol, ul, pl, gammal, pinfl
  real(dp)::ror, ur, pr, gammar, pinfr
  real(dp)::rog, ug, pg, gammag, pinfg
  real(dp)::dir
  integer::i,iv

  ! Pressure, density and velocity
  do i=1,ngrid

     rol    = max(qleft (i,1),smallr)
     ul     = qleft (i,2)
     pl     = qleft (i,3)
     gammal = 1.0d0/qleft(i,4)+1.0d0
     pinfl  = qleft(i,5)/qleft(i,4)/gammal

     ror    = qright(i,1)
     ur     = qright(i,2)
     pr     = qright(i,3)
     gammar = 1.0d0/qright(i,4)+1.0d0
     pinfr  = qright(i,5)/qright(i,4)/gammar

     call riemsti(ror,rol,ur,ul,pr,pl,gammar,gammal,pinfr,pinfl, &
          & rog,ug,pg,gammag,pinfg,dir)

     qgdnv(i,1) = rog
     qgdnv(i,2) = ug
     qgdnv(i,3) = pg
     qgdnv(i,4) = 1.0d0/(gammag-1.0d0)
     qgdnv(i,5) = gammag*pinfg/(gammag-1.0d0)
     
     ! Passive scalar
     do iv=6,nvar
        if(dir > 0.)then
           qgdnv(i,iv)=qleft(i,iv)
        else
           qgdnv(i,iv)=qright(i,iv)
        end if
     end do

  end do

end subroutine riemann_approx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! Exact Riemann solver
subroutine riemsti(ror,rol,ur,ul,pr,pl,gammar,gammal,pinfr,pinfl,&
     &             rog,ug,pg,gammag,pinfg,dir)
  !     
  use amr_parameters,ONLY:dp
  use hydro_parameters,ONLY:smallc,smallr
  implicit real(dp)(a-h,o-z)
  real(dp),save::eps=1.d-5,v=0.d0 
  !     
  pstar=(pr+pl)*0.5d0
  ar=sqrt(gammar*(pr+pinfr)/ror)
  al=sqrt(gammal*(pl+pinfl)/rol)

  do
     rapr=(pstar+pinfr)/(pr+pinfr) 
     rapl=(pstar+pinfl)/(pl+pinfl) 
     call phi(rapr,gammar,res) 
     dmr=sqrt(max(gammar*ror*(pr+pinfr),smallc**2*ror**2))*res 
     call phi(rapl,gammal,res) 
     dml=sqrt(max(gammal*rol*(pl+pinfl),smallc**2*rol**2))*res 
     pnew=(ul-ur)/(1.d0/dml+1.d0/dmr)+(dml*pr+dmr*pl)/(dml+dmr) 
     pnew=max(pnew,eps*abs(pl+pr))
     write(*,*)pnew
     if(abs(pnew-pstar).gt.eps*abs(pl+pr))then 
        pstar=pnew 
     else
        exit
     endif
  end do
  !
  ustar=(pl-pr+dml*ul+dmr*ur)/(dml+dmr) 
  !   
  if(rapr.gt.1.d0)then 
     !     choc a dte  
     sigmar=ur+dmr/ror 
     rstarr=-dmr/(ustar-sigmar) 
  else 
     !     detente a dte 
     rstarr=ror*((pstar+pinfr)/(pr+pinfr))**(1.d0/gammar) 
  endif
  astarr=sqrt(max(gammar*(pstar+pinfr)/rstarr,smallc**2))    
  if(rapl.gt.1.d0)then 
     !     choc a gauche            
     sigmal=ul-dml/rol 
     rstarl=dml/(ustar-sigmal) 
  else 
     !     detente a gauche 
     rstarl=rol*((pstar+pinfl)/(pl+pinfl))**(1.d0/gammal) 
  endif
  astarl=sqrt(max(gammal*(pstar+pinfl)/rstarl,smallc**2))  
  !     
  !     affectation de l'etat calcule au noeud intermediaire        
  !     considere                                                                                                                         
  !     cas numero 1 . onde de choc faisant face a gauche                                                                                  
 
  if(v.lt.ustar.and.pl.le.pstar)then                       
     if(v.lt.sigmal)then                                  
        p=pl                                             
        a=al                                             
        u=ul                                             
        ro=rol                                           
     else                                               
        p=pstar                                          
        a=astarl                                         
        u=ustar                                 
        ro=rstarl                                       
     endif
     gamma=gammal    
     pinf=pinfl                                           
     dir=+1.
  endif
  !     cas numero 2 . onde de detente faisant face a gauche                                                                              
  if(v.lt.ustar.and.pl.gt.pstar)then                       
     dif=ul-al                                            
     if (v.le.dif)then                                    
        p=pl                                              
        a=al                                              
        u=ul                                              
        ro=rol                                            
     else                                                 
        dif=ustar-astarl 
        if(v.lt.dif)then                               
           u=((gammal-1.d0)*ul+2.d0*(v+al))/(gammal+1.d0)                 
           a=u-v                                           
           p1=pl*(a/al)**(2.d0*gammal/(gammal-1.d0))          
           ro1=gammal*p/a/a 
           ro=((a**2*rol**gammal)/((pinfl+pl)*gammal))**(1.d0/(gammal-1.d0)) 
           p=a**2*ro/gammal-pinfl 
        else                                              
           p=pstar                                        
           a=astarl                                       
           u=ustar                                     
           ro=rstarl                                     
        endif
     endif
     gamma=gammal 
     pinf=pinfl                                
     dir=+1.
  endif
  !     cas numero 3 . onde de choc faisant face a droite 
  if(ustar.le.v.and.pstar.ge.pr)then                       
     if(v.le.sigmar)then                                   
        p=pstar                                          
        a=astarr                                         
        u=ustar                                             
        ro=rstarr                                       
     else                                               
        p=pr                                             
        a=ar                                             
        u=ur                                             
        ro=ror                                           
     endif
     gamma=gammar 
     pinf=pinfr
     dir=-1.
  endif
  !     cas numero 4 . onde de detente faisant face a droite  
  if(ustar.le.v.and.pstar.lt.pr)then                       
     som=ustar+astarr                                     
     if(v.le.som)then                                    
        p=pstar                                           
        a=astarr                                          
        u=ustar                                         
        ro=rstarr                                        
     else                                                 
        som2=ur+ar  
        if(v.lt.som2)then                               
           u=((gammar-1.d0)*ur+2.d0*(v-ar))/(gammar+1.d0)                
           a=v-u                                          
           ro=((a**2*ror**gammar)/((pinfr+pr)*gammar))**(1.d0/(gammar-1.d0)) 
           p=a**2*ro/gammar-pinfr 
        else                                            
           p=pr                                          
           a=ar                                          
           u=ur                                          
           ro=ror                                        
        endif
     endif
     gamma=gammar 
     pinf=pinfr
     dir=-1.
  endif
  !     calcul des etats de Godunov
  gammag = gamma
  pinfg = pinf
  pg = p
  rog = ro
  ug = u
  return 
end subroutine riemsti
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine phi(r,gamma,res) 
  use amr_parameters
  implicit real(dp) (a-h,o-z)
  !
  un=1.d0-1.d-4  
  fact=(gamma-1.d0)/2.d0/gamma   
  if(r.lt.un) then 
     res=fact*(1.d0-r)/(1.d0-r**fact) 
  else 
     res=sqrt((gamma+1.d0)/2.d0/gamma*r+fact )
  endif
  return 
end subroutine phi
!-----------------------------------------------------------------------
