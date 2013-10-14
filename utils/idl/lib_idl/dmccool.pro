function dmccool, x, t

;if x lt 1.0e-3 then x = 1.0e-3

; MacDonald & Bailey table
xmat = [-0.133, 0.105, 0.452, 0.715, 0.901 $
		      , 1.030, 1.082, 1.174, 1.257, 1.362 $
		      , 1.448, 1.523, 1.569, 1.582, 1.539 $
		      , 1.430, 1.275, 1.168, 1.092, 1.019 $
		      , 1.000, 1.004, 1.008, 0.987, 0.905 $
		      , 0.738, 0.603, 0.555, 0.552, 0.554 $
		      , 0.552, 0.535, 0.425, 0.275, 0.251 $
		      , 0.232, 0.247, 0.283, 0.322, 0.363 $
		      , 0.397]

ncool=n_elements(x*t)

; electron impact excitation luminosity (eqn 3-10)
le1 = dblarr(ncool)
if n_elements(t gt 1) then begin
	idx10=where(t gt 10)
	if (idx10[0] ne -1) then le1[idx10] = le1[idx10] + $
		2.96d-23/sqrt(T[idx10]) * exp(-92.0/T[idx10])
	idx50=where(t gt 50)
	if (idx50[0] ne -1) then le1[idx50] = le1[idx50] + $
		6.08d-23/sqrt(T[idx50]) * exp(-413.0/T[idx50]) $
		+ 3.52e-23/sqrt(T[idx50]) * (exp(-554.0/T[idx50]) + $
		1.3*exp(-961.0/T[idx50]))
	idx2d4=where(t gt 2d4)
	if (idx2d4[0] ne -1) then le1[idx2d4] = le1[idx2d4] + $
		4.14e-26*sqrt(T[idx2d4]) * exp(-22700.0/T[idx2d4]) + $
		7.13e-26*sqrt(T[idx2d4])*(1.0-2.7e-9*T[idx2d4]^2) * $
		exp(-27700.0/T[idx2d4])
endif else begin
	if (t gt 10) then le1 = le1 + $
		2.96d-23/sqrt(T) * exp(-92.0/T)
	if (t gt 50) then le1 = le1 + $
		6.08d-23/sqrt(T) * exp(-413.0/T) $
		+ 3.52e-23/sqrt(T) * (exp(-554.0/T) + $
		1.3*exp(-961.0/T))
	if (t gt 2d4) then le1 = le1 + $
		4.14e-26*sqrt(T) * exp(-22700.0/T) + $
		7.13e-26*sqrt(T)*(1.0-2.7e-9*T^2) * $
		exp(-27700.0/T)
endelse

; Hydrogen cooling
lh = dblarr(ncool)
if (idx50[0] ne -1) then lh[idx50] = $
	2.37e-27*exp(-413/T[idx50]) + 3.52e-27*(exp(-554.0/T[idx50]) $
	+ 1.4*exp(-961.0/T[idx50]))

; Neutral cooling
u = dblarr(ncool)
idxlo=where(T/157890. lt 3.16)
idxhi=where(T/157890. ge 3.16)
if (idxlo[0] ne -1) then u[idxlo] = T[idxlo]/157890.
if (idxhi[0] ne -1) then u[idxhi] = 3.16
u2 = u^2
om=.6098+1.489*u+.50755*u2-.38145*u*u2+.10196*u2*u2 $
    -.01007*u*u2*u2
dom=(1.489+2.*.50755*u-3.*.38145*u2+4.*.10196*u2*u $
       -5.*.01007*u2*u2)/157890.
p1 = dblarr(ncool)
idx1d4=where(t gt 1d4)
if (idx1d4[0] ne -1) then p1[idx1d4] = 0.5*1.41e-16*om * $
	exp(-118000./T[idx1d4]) / sqrt(T[idx1d4])

; Cooling in various regimes: < 100 K
cool = dblarr(ncool)
idx100 = where(t lt 100)
idx1d4 = where((t lt 1d4) and (t ge 100))
idxhi = where(t gt 1.27717d8)
if n_elements(x) eq 1 then begin
    if (idx100[0] ne -1) then cool[idx100] = x*le1[idx100] + lh[idx100] + $
	(1.0-x)*p1[idx100]
    if (idx1d4[0] ne -1) then cool[idx1d4] = 1.0d-23*x*2.8347e-10* $
	(T[idx1d4] - 1.0e+02)^2.3562 + x*le1[idx1d4] + lh[idx1d4] + $
	(1.0-x)*p1[idx1d4]
    if (idxhi[0] ne -1) then cool[idxhi] = x*2.3988d-4*sqrt(t[idxhi])
endif else begin
    if (idx100[0] ne -1) then cool[idx100] = x[idx100]*le1[idx100] + $
	lh[idx100] + (1.0-x[idx100])*p1[idx100]
    if (idx1d4[0] ne -1) then cool[idx1d4] = 1.0d-23*x[idx1d4]*2.8347e-10* $
	(T[idx1d4] - 1.0e+02)^2.3562 + x[idx1d4]*le1[idx1d4] + lh[idx1d4] + $
	(1.0-x[idx1d4])*p1[idx1d4]
    if (idxhi[0] ne -1) then cool[idxhi] = x[idxhi]*2.3988d-4*sqrt(t[idxhi])
endelse

idxmid=where((t ge 1d4) and (t le 1.27717d8))
if (idxmid[0] ne -1) then begin
	tlost = alog10(T[idxmid])
	ipps = floor(10.0*tlost) - 38
	idx=where(ipps gt 41)
	if (idx[0] ne -1) then ipps[idx] = 41
	idx=where(ipps lt 2)
	jaug=ipps
	if (idx[0] ne -1) then jaug[idx] = 2
	qq2 = 3.8 + 0.1*jaug
	qt2 = tlost - qq2
	qt3 = qt2 - 0.1
	qt1 = qt2 + 0.1
	qt4 = qt3 - 0.1

	xu1 = qt2*qt3*qt4/6.0e-03
	xu2 = qt1*qt3*qt4/2.0e-03
	xu3 = qt1*qt2*qt4/2.0e-03
	xu4 = qt1*qt2*qt3/6.0e-03

	tcool = -xmat[jaug-3]*xu1 + xmat[jaug-2]*xu2 - $
		xmat[jaug-1]*xu3 + xmat[jaug]*xu4
	cool[idxmid] = 1.0d-23*10.0^tcool*x[idxmid] + (1.0-x[idxmid])*p1
endif

return, cool

end
