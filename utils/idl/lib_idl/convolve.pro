function convolve, expot, nu, b, boxsize
; Does Weinberg-Gunn convolution on expotential

table=dblarr(boxsize,boxsize,boxsize)
res=dblarr(boxsize,boxsize,boxsize)
u=dblarr(boxsize,boxsize,boxsize)

; load table
x=dindgen(boxsize)-(boxsize-1.0)/2.0
table[*,0,0]=exp(-x*x/(4*b*nu))
for j=0,boxsize-1 do begin
	for k=0,boxsize-1 do begin
		table[*,j,k]=table[*,0,0]
	endfor
endfor

; x convolution
for n=boxsize/2,3*boxsize/2-1 do begin
	if n-boxsize+1 lt 0 then begin
		min=0
	endif else begin
		min=n-boxsize+1
	endelse
	if n lt boxsize then begin
		max=n+1
	endif else begin
		max=boxsize
	endelse
	for i=min,max-1 do begin
		res[n-boxsize/2,*,*] = res[n-boxsize/2,*,*] + $
			table[i,*,*]*expot[n-i,*,*]
	endfor
endfor
u=res

; y convolution
res[*,*,*]=0.0
for n=boxsize/2,3*boxsize/2-1 do begin
	if n-boxsize+1 lt 0 then begin
		min=0
	endif else begin
		min=n-boxsize+1
	endelse
	if n lt boxsize then begin
		max=n+1
	endif else begin
		max=boxsize
	endelse
	for i=min,max-1 do begin
		res[*,n-boxsize/2,*] = res[*,n-boxsize/2,*] + $
			table[i,*,*]*u[*,n-i,*]
	endfor
endfor
u=res

; z convolution
res[*,*,*]=0.0
for n=boxsize/2,3*boxsize/2-1 do begin
	if n-boxsize+1 lt 0 then begin
		min=0
	endif else begin
		min=n-boxsize+1
	endelse
	if n lt boxsize then begin
		max=n+1
	endif else begin
		max=boxsize
	endelse
	for i=min,max-1 do begin
		res[n-boxsize/2,*,*] = res[n-boxsize/2,*,*] + $
			table[i,*,*]*u[n-i,*,*]
	endfor
endfor
u=res

return, u
end







