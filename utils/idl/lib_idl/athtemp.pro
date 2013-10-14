function athtemp, d, px, py, pz, e, bx=bx, by=by, bz=bz, $
	gamma=gamma, kb=kb, mu=mu, d_n=d_n, mh=mh, cs=cs, va=va, vel=vel, p=p

if not keyword_set(gamma) then gamma=5./3.
if not keyword_set(kb) then kb=1.38d-16
if not keyword_set(mh) then mh=2.34d-24
if not keyword_set(mu) then mu=2.33*1.67d-24

vx=px/d
vy=py/d
vz=pz/d
vel=sqrt(vx^2+vy^2+vz^2)

eth=e-0.5*d*vel^2
if n_elements(bx) ne 0 then begin
	eb=0.5*(bx^2+by^2+bz^2)
	eth=eth-eb
	va=sqrt(2*eb/d)
endif
esp=eth/d
p=(gamma-1)*eth
cs=sqrt((gamma-1)*esp)

if n_elements(d_n) ne 0 then begin
	x=1-d_n/d
	mpart = (1.0-x)*mu + x*0.5*mh
	t = (gamma-1)*mpart*esp/kb
endif else begin
	t=(gamma-1)*mu*esp/kb
endelse

return, t

end

