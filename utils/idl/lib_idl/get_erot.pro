function get_erot, rho, vx, vy, vz, ncell, halfboxlen, ke=ke, cm=cm, $
	periodic=periodic

l=dblarr(3)
r=dblarr(3)
irot=0.0d
ke=0.0d
if not keyword_set(cm) then boxctr=findmaxnd(rho) $
else boxctr=findcm(rho)
for i=0,ncell-1 do begin
   if keyword_set(periodic) then begin
      if abs(i-boxctr[0]) le ncell/2 then ieff = i $
      else if i gt boxctr[0] then ieff = i - ncell $
      else ieff = i + ncell
   endif else ieff = i
   r[0] = (ieff-boxctr[0])*2.0d*halfboxlen/ncell
   for j=0,ncell-1 do begin
      if keyword_set(periodic) then begin
         if abs(j-boxctr[1]) le ncell/2 then jeff = j $
         else if j gt boxctr[1] then jeff = j - ncell $
         else jeff = j + ncell
      endif else jeff = j
      r[1] = (jeff-boxctr[1])*2.0d*halfboxlen/ncell
      for k=0,ncell-1 do begin
        if keyword_set(periodic) then begin
	   if abs(k-boxctr[2]) le ncell/2 then keff = k $
	   else if k gt boxctr[2] then keff = k - ncell $
	   else keff = k + ncell
        endif else keff = k
	r[2] = (keff-boxctr[2])*2.0d*halfboxlen/ncell

	m=(2.0d*halfboxlen/ncell)^3 * rho[i,j,k]
	l[0] = l[0] + m * (r[1]*vz[i,j,k] - r[2]*vy[i,j,k])
	l[1] = l[1] + m * (r[2]*vx[i,j,k] - r[0]*vz[i,j,k])
	l[2] = l[2] + m * (r[0]*vy[i,j,k] - r[1]*vx[i,j,k])
	irot = irot + m * (r[0]^2+r[1]^2+r[2]^2)
	ke = ke + 0.5 * m * $
		(vx[i,j,k]^2+vy[i,j,k]^2+vz[i,j,k]^2)
      endfor
   endfor
endfor

return, (l[0]^2+l[1]^2+l[2]^2) / (2.0*irot)

end
