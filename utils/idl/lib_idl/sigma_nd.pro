function sigma_nd, arr, ctr=ctr

; function to return the second moment tensor of a 1- to 3-dimensional
; array. The keyword ctr specifies the center about which the moment
; is to be computed.

sz=size(arr)
ndim=sz[0]
if not keyword_set(ctr) then ctr=fltarr(ndim)
if ndim eq 0 then return, 0
if ndim eq 1 then return, total((dindgen(sz[1])-ctr)^2*arr)/total(arr)
if ndim eq 2 then begin
	sigma=fltarr(2,2)
	vec=dindgen(sz[2])
	for i=0L,sz[1]-1 do begin
		sigma[0,0] = sigma[0,0] + (i-ctr[0])^2*total(arr[i,vec])
		sigma[1,1] = sigma[1,1] + total((vec-ctr[1])^2*arr[i,vec])
		sigma[0,1] = sigma[0,1] + $
			(i-ctr[0])*total((vec-ctr[1])*arr[i,vec])
	endfor
	sigma = sigma/total(arr)
	sigma[1,0] = sigma[0,1]
	return, sigma
endif
if ndim eq 3 then begin
	sigma=fltarr(3,3)
	vec=dindgen(sz[3])
	for i=0L,sz[1]-1 do begin
		for j=0L,sz[2]-1 do begin
			sigma[0,0] = sigma[0,0] + $
				(i-ctr[0])^2*total(arr[i,j,vec])
			sigma[0,1] = simga[0,1] + $
				(i-ctr[0])*(j-ctr[1])*total(arr[i,j,vec])
			sigma[0,2] = simga[0,2] + $
				(i-ctr[0])*total((vec-ctr[2])*arr[i,j,vec])
			sigma[1,1] = sigma[1,1] + $
				(j-ctr[1])^2*total(arr[i,j,vec])
			sigma[1,2] = simga[1,2] + $
				(j-ctr[1])*total((vec-ctr[2])*arr[i,j,vec])
			sigma[2,2] = sigma[2,2] + $
				total((vec-ctr[2])^2*arr[i,j,vec])
		endfor
	endfor
	simga = sigma/total(arr)
	sigma[1,0] = sigma[0,1]
	sigma[2,0] = sigma[0,2]
	sigma[2,1] = sigma[1,2]
	return, sigma
endif
if ndim gt 3 then return, -1

end


