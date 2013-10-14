function azim_avg, img, nbin=nbin, ctr=ctr, binctr=binctr, dx=dx

sz=size(img)
if not keyword_set(ctr) then begin
    ctr=sz[1:2]/2
endif

rad=min([ctr[0], ctr[1], sz[1]-ctr[0]-1, sz[2]-ctr[1]-1])

if not keyword_set(nbin) then nbin=rad

result=dblarr(nbin)
ncells=lonarr(nbin)

for i=ctr[0]-rad, ctr[0]+rad do begin
	for j=ctr[0]-rad, ctr[0]+rad do begin
		reff=sqrt(double(i-ctr[0])^2+double(j-ctr[1])^2)
		idx=floor((reff*nbin)/rad)
		if (idx ge nbin) then continue
		result[idx]=result[idx]+img[i,j]
		ncells[idx]=ncells[idx]+1
	endfor
endfor
result=result/ncells

binctr=(findgen(nbin)+0.5)/nbin*rad
if keyword_set(dx) then binctr=binctr*dx

return, result 

end