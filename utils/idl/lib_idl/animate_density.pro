pro animate_density, rho, idx, slice, rhomin=rhomin, rhomax=rhomax

sz=size(rho)
xinteranimate, set=[500,500,sz[1]]

if (idx eq 0) then begin
	dim1=sz[2]
	dim2=sz[3]
endif else if (idx eq 1) then begin
	dim1=sz[2]
	dim2=sz[4]
endif else begin
	dim1=sz[3]
	dim2=sz[4]
endelse
rhoimg=fltarr(dim1,dim2)
for n=0,sz[1]-1 do begin
	if (idx eq 0) then $
		for i=0,dim1-1 do for j=0,dim2-1 do $
			rhoimg[i,j]=rho[n,i,j,slice] $
	else if (idx eq 1) then $
		for i=0,dim1-1 do for j=0,dim2-1 do $
			rhoimg[i,j]=rho[n,i,slice,j] $
	else $
		for i=0,dim1-1 do for j=0,dim2-1 do $
			rhoimg[i,j]=rho[n,i,j,slice]
	rhoimg1=congrid(rhoimg,500,500)

	; get color table
	tvlct, r, g, b, /get
	ncolor=n_elements(r)
	if not keyword_set(rhomin) then rhomin=min(rhoimg1)
	if not keyword_set(rhomax) then rhomax=max(rhoimg1)
	rhoimg1=alog10(rhoimg1)
	lrhomin=alog10(rhomin)
	lrhomax=alog10(rhomax)
	rhoimg1=fix((rhoimg1-lrhomin)/(lrhomax-lrhomin)*ncolor)
	tv, rhoimg1

	xinteranimate, frame=n, window=[0,0,0,500,500]
endfor

xinteranimate

end