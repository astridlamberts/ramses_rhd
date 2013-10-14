function grad, phi, periodic=periodic
; Input must be a 3 dimensional array of size [i,j,k], with all of
; i, j, and k greater than or equal to 3. Output is a 4 dimensional
; array of size [3, i-2, j-2, k-2]. Spacing between input array
; elements is taken to be 1.

sz=size(phi)
if not keyword_set(periodic) then begin
	xidx=lindgen(sz[1]-2)+1
	xidxp1=xidx+1
	xidxm1=xidx-1
	yidx=lindgen(sz[2]-2)+1
	yidxp1=yidx+1
	yidxm1=yidx-1
	zidx=lindgen(sz[3]-2)+1
	zidxp1=zidx+1
	zidxm1=zidx-1
	grad_phi=fltarr(3, sz[1]-2,sz[2]-2, sz[3]-2)
endif else begin
	xidx=lindgen(sz[1])
	xidxp1=xidx+1
	xidxp1[sz[1]-1]=0
	xidxm1=xidx-1
	xidxm1[0]=sz[1]-1
	yidx=lindgen(sz[2])
	yidxp1=yidx+1
	yidxp1[sz[2]-1]=0
	yidxm1=yidx-1
	yidxm1[0]=sz[2]-1
	zidx=lindgen(sz[3])
	zidxp1=zidx+1
	zidxp1[sz[3]-1]=0
	zidxm1=zidx-1
	zidxm1[0]=sz[3]-1
	grad_phi=fltarr(3, sz[1], sz[2], sz[3])
endelse

tmp1=phi[xidxp1,yidx,*]
tmp1=tmp1[*,*,zidx]
tmp2=phi[xidxm1,yidx,*]
tmp2=tmp2[*,*,zidx]
grad_phi[0,*,*,*]=(tmp1-tmp2) / 2.0
tmp1=phi[xidx,yidxp1,*]
tmp1=tmp1[*,*,zidx]
tmp2=phi[xidx,yidxm1,*]
tmp2=tmp2[*,*,zidx]
grad_phi[1,*,*,*]=(tmp1-tmp2) / 2.0
tmp1=phi[xidx,yidx,*]
tmp1=tmp1[*,*,zidxp1]
tmp2=phi[xidx,yidx,*]
tmp2=tmp2[*,*,zidxm1]
grad_phi[2,*,*,*]=(tmp1-tmp2) / 2.0

return, grad_phi
end

function div, psi, periodic=periodic
; The routine computes the divergence of the input vector field using
; a two-sided difference. The input must be a 4 dimensional array of
; size [3,i,j,k], with i, j, and k all greater than or equal to 3.
; Output is an array of size [i-2, j-2, k-2] if periodic is not set, or
; size [i,j,k] if it is. Spacing between input array elements is taken
; to be 1. Setting periodic causes the routine to assume periodic
; boundary conditions.
sz=size(psi)
if not keyword_set(periodic) then begin
	xidx=lindgen(sz[2]-2)+1
	xidxp1=xidx+1
	xidxm1=xidx-1
	yidx=lindgen(sz[3]-2)+1
	yidxp1=yidx+1
	yidxm1=yidx-1
	zidx=lindgen(sz[4]-2)+1
	zidxp1=zidx+1
	zidxm1=zidx-1
endif else begin
	xidx=lindgen(sz[2])
	xidxp1=xidx+1
	xidxp1[sz[2]-1]=0
	xidxm1=xidx-1
	xidxm1[0]=sz[2]-1
	yidx=lindgen(sz[3])
	yidxp1=yidx+1
	yidxp1[sz[3]-1]=0
	yidxm1=yidx-1
	yidxm1[0]=sz[3]-1
	zidx=lindgen(sz[4])
	zidxp1=zidx+1
	zidxp1[sz[4]-1]=0
	zidxm1=zidx-1
	zidxm1[0]=sz[4]-1
endelse
div_psi = reform( $
	psi[0, xidxp1, yidx, zidx] - $
	psi[0, xidxm1, yidx, zidx] + $
	psi[1, xidx, yidxp1, zidx] - $
	psi[1, xidx, yidxm1, zidx] + $
	psi[2, xidx, yidx, zidxp1] - $
	psi[2, xidx, yidx, zidxm1])
div_psi=div_psi/2.0
return, div_psi
end

function curl, psi, periodic=periodic
; Input must be a 4 dimensional array of size [3,i,j,k], with i, j,
; and k all greater than or equal to 3. Output is an array of size
; [3, i-2, j-2, k-2], or [3,i,j,k] is periodic is set. Spacing
; between input array elements is taken to be 1.

sz=size(psi)
if not keyword_set(periodic) then begin
    curl_psi=fltarr(3, sz[2]-2, sz[3]-2, sz[4]-2)
    curl_psi[0,*,*,*]=(psi[2,1:sz[2]-2,2:sz[3]-1,1:sz[4]-2] - $
	psi[2,1:sz[2]-2,0:sz[3]-3,1:sz[4]-2] - $
	psi[1,1:sz[2]-2,1:sz[3]-2,2:sz[4]-1] + $
	psi[1,1:sz[2]-2,1:sz[3]-2,0:sz[4]-3]) / 2.0
    curl_psi[1,*,*,*]= (psi[0,1:sz[2]-2,1:sz[3]-2,2:sz[4]-1] - $
	psi[0,1:sz[2]-2,1:sz[3]-2,0:sz[4]-3] - $
	psi[2,2:sz[2]-1,1:sz[3]-2,1:sz[4]-2] + $
	psi[2,0:sz[2]-3,1:sz[3]-2,1:sz[4]-2]) / 2.0
    curl_psi[2,*,*,*]= (psi[1,2:sz[2]-1,1:sz[3]-2,1:sz[4]-2] - $
	psi[1,0:sz[2]-3,1:sz[3]-2,1:sz[4]-2] - $
	psi[0,1:sz[2]-2,2:sz[3]-1,1:sz[4]-2] + $
	psi[0,1:sz[2]-2,0:sz[3]-3,1:sz[4]-2]) / 2.0
endif else begin
    curl_psi=fltarr(3, sz[2], sz[3], sz[4])

    ; points where no periodic wrap is needed
    curl_psi[0,*,1:sz[3]-2,1:sz[4]-2] = $
	(psi[2,*,2:sz[3]-1,1:sz[4]-2] - $
	psi[2,*,0:sz[3]-3,1:sz[4]-2] - $
	psi[1,*,1:sz[3]-2,2:sz[4]-1] + $
	psi[1,*,1:sz[3]-2,0:sz[4]-3]) / 2.0
    curl_psi[1,1:sz[2]-2,*,1:sz[4]-2] = $
	(psi[0,1:sz[2]-2,*,2:sz[4]-1] - $
	psi[0,1:sz[2]-2,*,0:sz[4]-3] - $
	psi[2,2:sz[2]-1,*,1:sz[4]-2] + $
	psi[2,0:sz[2]-3,*,1:sz[4]-2]) / 2.0
    curl_psi[2,1:sz[2]-2,1:sz[3]-2,*] = $
	(psi[1,2:sz[2]-1,1:sz[3]-2,*] - $
	psi[1,0:sz[2]-3,1:sz[3]-2,*] - $
	psi[0,1:sz[2]-2,2:sz[3]-1,*] + $
	psi[0,1:sz[2]-2,0:sz[3]-3,*]) / 2.0

    ; periodic wrap in x direction
    curl_psi[1,0,*,1:sz[4]-2] = $
	(psi[0,0,*,2:sz[4]-1] - $
	psi[0,0,*,0:sz[4]-3] - $
	psi[2,1,*,1:sz[4]-2] + $
	psi[2,sz[2]-1,*,1:sz[4]-2]) / 2.0
    curl_psi[2,0,1:sz[3]-2,*] = $
	(psi[1,1,1:sz[3]-2,*] - $
	psi[1,sz[2]-1,1:sz[3]-2,*] - $
	psi[0,0,2:sz[3]-1,*] + $
	psi[0,0,0:sz[3]-3,*]) / 2.0
    curl_psi[1,sz[2]-1,*,1:sz[4]-2] = $
	(psi[0,sz[2]-1,*,2:sz[4]-1] - $
	psi[0,sz[2]-1,*,0:sz[4]-3] - $
	psi[2,0,*,1:sz[4]-2] + $
	psi[2,sz[2]-2,*,1:sz[4]-2]) / 2.0
    curl_psi[2,sz[2]-1,1:sz[3]-2,*] = $
	(psi[1,0,1:sz[3]-2,*] - $
	psi[1,sz[2]-2,1:sz[3]-2,*] - $
	psi[0,sz[2]-1,2:sz[3]-1,*] + $
	psi[0,sz[2]-1,0:sz[3]-3,*]) / 2.0

    ; periodic wrap in y direction
    curl_psi[0,*,0,1:sz[4]-2] = $
	(psi[2,*,1,1:sz[4]-2] - $
	psi[2,*,sz[3]-1,1:sz[4]-2] - $
	psi[1,*,0,2:sz[4]-1] + $
	psi[1,*,0,0:sz[4]-3]) / 2.0
    curl_psi[2,1:sz[2]-2,0,*] = $
	(psi[1,2:sz[2]-1,0,*] - $
	psi[1,0:sz[2]-3,0,*] - $
	psi[0,1:sz[2]-2,1,*] + $
	psi[0,1:sz[2]-2,sz[3]-1,*]) / 2.0
    curl_psi[0,*,sz[3]-1,1:sz[4]-2] = $
	(psi[2,*,0,1:sz[4]-2] - $
	psi[2,*,sz[3]-2,1:sz[4]-2] - $
	psi[1,*,sz[3]-1,2:sz[4]-1] + $
	psi[1,*,sz[3]-1,0:sz[4]-3]) / 2.0
    curl_psi[2,1:sz[2]-2,sz[3]-1,*] = $
	(psi[1,2:sz[2]-1,sz[3]-1,*] - $
	psi[1,0:sz[2]-3,sz[3]-1,*] - $
	psi[0,1:sz[2]-2,0,*] + $
	psi[0,1:sz[2]-2,sz[3]-2,*]) / 2.0

    ; periodic wrap in z direction
    curl_psi[0,*,1:sz[3]-2,0] = $
	(psi[2,*,2:sz[3]-1,0] - $
	psi[2,*,0:sz[3]-3,0] - $
	psi[1,*,1:sz[3]-2,1] + $
	psi[1,*,1:sz[3]-2,sz[4]-1]) / 2.0
    curl_psi[1,1:sz[2]-2,*,0] = $
	(psi[0,1:sz[2]-2,*,1] - $
	psi[0,1:sz[2]-2,*,sz[4]-1] - $
	psi[2,2:sz[2]-1,*,0] + $
	psi[2,0:sz[2]-3,*,0]) / 2.0
    curl_psi[0,*,1:sz[3]-2,sz[4]-1] = $
	(psi[2,*,2:sz[3]-1,sz[4]-1] - $
	psi[2,*,0:sz[3]-3,sz[4]-1] - $
	psi[1,*,1:sz[3]-2,0] + $
	psi[1,*,1:sz[3]-2,sz[4]-2]) / 2.0
    curl_psi[1,1:sz[2]-2,*,sz[4]-1] = $
	(psi[0,1:sz[2]-2,*,0] - $
	psi[0,1:sz[2]-2,*,sz[4]-2] - $
	psi[2,2:sz[2]-1,*,sz[4]-1] + $
	psi[2,0:sz[2]-3,*,sz[4]-1]) / 2.0

    ; periodic wrap in xy direction
    curl_psi[2,0,0,*] = $
	(psi[1,1,0,*] - $
	psi[1,sz[2]-1,0,*] - $
	psi[0,0,1,*] + $
	psi[0,0,sz[3]-1,*]) / 2.0
    curl_psi[2,sz[2]-1,0,*] = $
	(psi[1,0,0,*] - $
	psi[1,sz[2]-2,0,*] - $
	psi[0,sz[2]-1,1,*] + $
	psi[0,sz[2]-1,sz[3]-1,*]) / 2.0
    curl_psi[2,0,sz[3]-1,*] = $
	(psi[1,1,sz[3]-1,*] - $
	psi[1,sz[2]-1,sz[3]-1,*] - $
	psi[0,0,0,*] + $
	psi[0,0,sz[3]-2,*]) / 2.0
    curl_psi[2,sz[2]-1,sz[3]-1,*] = $
	(psi[1,0,sz[3]-1,*] - $
	psi[1,sz[2]-2,sz[3]-1,*] - $
	psi[0,sz[2]-1,0,*] + $
	psi[0,sz[2]-1,sz[3]-2,*]) / 2.0

    ; periodic wrap in xz direction
    curl_psi[1,0,*,0] = $
	(psi[0,0,*,1] - $
	psi[0,0,*,sz[4]-1] - $
	psi[2,1,*,0] + $
	psi[2,sz[2]-1,*,0]) / 2.0
    curl_psi[1,sz[2]-1,*,0] = $
	(psi[0,sz[2]-1,*,1] - $
	psi[0,sz[2]-1,*,sz[4]-1] - $
	psi[2,0,*,0] + $
	psi[2,sz[2]-2,*,0]) / 2.0
    curl_psi[1,0,*,sz[4]-1] = $
	(psi[0,0,*,0] - $
	psi[0,0,*,sz[4]-2] - $
	psi[2,1,*,sz[4]-1] + $
	psi[2,sz[2]-1,*,sz[4]-1]) / 2.0
    curl_psi[1,sz[2]-1,*,sz[4]-1] = $
	(psi[0,sz[2]-1,*,0] - $
	psi[0,sz[2]-1,*,sz[4]-2] - $
	psi[2,0,*,sz[4]-1] + $
	psi[2,sz[2]-2,*,sz[4]-1]) / 2.0

    ; periodic wrap in yz direction
    curl_psi[0,*,0,0] = $
	(psi[2,*,1,0] - $
	psi[2,*,sz[3]-1,0] - $
	psi[1,*,0,1] + $
	psi[1,*,0,sz[4]-1]) / 2.0
    curl_psi[0,*,sz[3]-1,0] = $
	(psi[2,*,0,0] - $
	psi[2,*,sz[3]-2,0] - $
	psi[1,*,sz[3]-1,1] + $
	psi[1,*,sz[3]-1,sz[4]-1]) / 2.0
    curl_psi[0,*,0,sz[4]-1] = $
	(psi[2,*,1,sz[4]-1] - $
	psi[2,*,sz[3]-1,sz[4]-1] - $
	psi[1,*,0,0] + $
	psi[1,*,0,sz[4]-2]) / 2.0
    curl_psi[0,*,sz[3]-1,sz[4]-1] = $
	(psi[2,*,0,sz[4]-1] - $
	psi[2,*,sz[3]-2,sz[4]-1] - $
	psi[1,*,sz[3]-1,0] + $
	psi[1,*,sz[3]-1,sz[4]-2]) / 2.0

endelse

return, curl_psi
end
