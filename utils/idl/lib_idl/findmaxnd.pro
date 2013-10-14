function findmaxnd, v
; find the maximum of an N-dimensional vector

sz=size(v)
res=lonarr(sz[0])
dummy=max(v, maxidx)
for i=long(sz[0]-1),0,-1 do begin
	div=1
	for j=1,i do div=div*sz[j]
	res[i]=maxidx/div
	maxidx=maxidx mod div
endfor

return, res
end