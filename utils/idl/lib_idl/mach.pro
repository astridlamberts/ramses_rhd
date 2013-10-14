
;program that computes the Mach number of a turbulent region with an
;isothermal soundspeed cs=1. We use the density-weighted velocity dispersion.
;x=0.0
;d = mass density 
;p = momentum density (density * speed)
moy=mean(d)
sigx=sqrt(mean((px*px)/d)/moy - (mean(px)/moy)^2)
sigy=sqrt(mean((py*py)/d)/moy - (mean(py)/moy)^2)
sigz=sqrt(mean((pz*pz)/d)/moy - (mean(pz)/moy)^2)

sigma=sqrt(sigx*sigx+sigy*sigy+sigz*sigz)

print,sigma

;M=sigma/cs
end



