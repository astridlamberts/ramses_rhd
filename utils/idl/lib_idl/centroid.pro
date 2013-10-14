function centroid, arr

; function to find the centroid of an 1- to 3-dimensional array

sz=size(arr)
ndim=sz[0]
if ndim eq 0 then return, 0
if ndim eq 1 then return, total(dindgen(sz[1])*arr)/total(arr)
if ndim eq 2 then begin
	ctr=fltarr(2)
	vec=dindgen(sz[2])
	for i=0L,sz[1]-1 do begin
		ctr[0] = ctr[0] + i*total(arr[i,vec])
		ctr[1] = ctr[1] + total(vec*arr[i,vec])
	endfor
	ctr = ctr/total(arr)
	return, ctr
endif
if ndim eq 3 then begin
	ctr=fltarr(3)
	vec=dindgen(sz[3])
	for i=0L,sz[1]-1 do begin
		for j=0L,sz[2]-1 do begin
			ctr[0] = ctr[0] + i*total(arr[i,j,vec])
			ctr[1] = ctr[1] + j*total(arr[i,j,vec])
			ctr[2] = ctr[2] + total(vec*arr[i,j,vec])
		endfor
	endfor
	ctr = ctr/total(arr)
	return, ctr
endif
if ndim gt 3 then return, -1

end

