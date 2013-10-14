function angsep, v1, v2

; Function to find the angular distance between two vectors or arrays
; of vectors. Arrays of vectors must be arrays of size [nvec,3]

sz=size(v1)
if sz[0] eq 1 then begin
	v1save=v1
	v1=reform(v1,1,3)
endif
sz=size(v2)
if sz[0] eq 1 then begin
	v2save=v2
	v2=reform(v2,1,3)
endif

mag1=sqrt(total(v1^2,2))
mag2=sqrt(total(v2^2,2))

v1dotv2=v1[*,0]*v2[*,0]+v1[*,1]*v2[*,1]+v1[*,2]*v2[*,2]
costh=v1dotv2/(mag1*mag2)
idx=where(costh gt 1.0)
if idx[0] ne -1 then costh[idx]=1.0
idx=where(costh lt -1.0)
if idx[0] ne -1 then costh[idx]=1.0

sz=size(v1save)
if sz[0] ne 0 then v1=v1save
sz=size(v2save)
if sz[0] ne 0 then v2=v2save

return, reform(acos(costh))

end
