pro solproj, ax, ay, az, px, py, pz, err=err, $
	geterr=geterr, printerr=printerr, ft=ft

sz=size(ax)
nx=sz[1]
ny=sz[2]
nz=sz[3]

; set up k vector
kx1d=sin(2*!pi*findgen(nx)/nx)
ky1d=sin(2*!pi*findgen(ny)/ny)
kz1d=sin(2*!pi*findgen(nz)/nz)
kx=fltarr(nx,ny,nz)
ky=fltarr(nx,ny,nz)
kz=fltarr(nx,ny,nz)
for j=0,ny-1 do for k=0,nz-1 do kx[*,j,k]=kx1d
for i=0,nx-1 do for k=0,nz-1 do ky[i,*,k]=ky1d
for i=0,nx-1 do for j=0,ny-1 do kz[i,j,*]=kz1d
knorm=sqrt(kx^2+ky^2+kz^2)
knorm[0,0,0]=1.0
kx=kx/knorm
ky=ky/knorm
kz=kz/knorm

; compute forward transform if needed
if keyword_set(ft) then begin
	axft=ax
	ayft=ay
	azft=az
endif else begin
	axft=fft(ax)
	ayft=fft(ay)
	azft=fft(az)
endelse

; compute dot product
aftdotk=axft*kx+ayft*ky+azft*kz

; do projection
axft=axft-aftdotk*kx
ayft=ayft-aftdotk*ky
azft=azft-aftdotk*kz

; transform back
px=fft(axft,/inverse)
py=fft(ayft,/inverse)
pz=fft(azft,/inverse)

; get error
if keyword_set(geterr) then begin
	p=fltarr(3,nx,ny,nz)
	p[0,*,*,*]=px
	p[1,*,*,*]=py
	p[2,*,*,*]=pz
	divp=div(p,/periodic)
	err=total(abs(divp))/total(sqrt(abs(ax)^2+abs(ay)^2+abs(az)^2))
	if keyword_set(printerr) then print, "Projection error = ", err
endif

end
