function angdist, v1, v2

; Function to find the angular distance between two vectors or arrays
; of vectors. Arrays of vectors must be arrays of size [3,nvec]

sz=size(v1)
if sz[0] eq 1 then begin
	v1save=v1
	v1=reform(v1,3,1)
endif
sz=size(v2)
if sz[0] eq 1 then begin
	v2save=v2
	v2=reform(v2,3,1)
endif

mag1=sqrt(total(v1^2,1))
mag2=sqrt(total(v2^2,1))

v1dotv2=v1[0,*]*v2[0,*]+v1[1,*]*v2[1,*]+v1[2,*]*v2[2,*]
costh=v1dotv2/(mag1*mag2)
idx=where(costh gt 1.0)
if idx[0] ne -1 then costh[idx]=1.0
idx=where(costh lt -1.0)
if idx[0] ne -1 then costh[idx]=1.0

return, acos(th)

end
