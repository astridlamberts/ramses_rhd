function findcm, rho

sz=size(rho)
cm=fltarr(3)

karr=lindgen(sz[3])
for i=0L,sz[1]-1 do begin
	for j=0L,sz[2]-1 do begin
		cm[0] = cm[0] + i*total(rho[i,j,karr])
		cm[1] = cm[1] + j*total(rho[i,j,karr])
		cm[2] = cm[2] + total(karr*rho[i,j,karr])
	endfor
endfor

cm = cm/total(rho)
return, cm

end