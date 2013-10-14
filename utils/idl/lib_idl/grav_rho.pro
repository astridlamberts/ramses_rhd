function grav_rho, rho, pot=pot

; routine to compute the gravitational binding energy of a density
; grid, rho, using a PM method

twopi=3.1415927*2.0
pe = 0.0d
sz=long(n_elements(rho[0,0,*]))

; start by taking the Fourier transform of the density field
rhotemp=double(rho-mean(rho))
rhotemp=complex(rhotemp)
rhoft=fft(rhotemp)*sz^3

; set up a table of trigonometric terms we'll need
trigtable=cos(twopi*findgen(sz)/sz)
ghatscale=twopi/sz^3

; now do the convolution
for k=0,sz-1 do begin
   for j=0,sz-1 do begin
	for i=0,sz/2+1 do begin
		ghat = 3. - trigtable[i] - trigtable[j] - trigtable[k]
		if ghat ne 0 then ghat = -ghatscale/ghat
		pe = pe + total(abs(rhoft[i,j,k])^2*ghat)
		rhoft[i,j,k] = rhoft[i,j,k]*ghat
		if (sz-i-1 gt 0) then rhoft[sz-i-1,j,k] = rhoft[sz-i-1,j,k]*ghat
	endfor
   endfor
endfor

if keyword_set(pot) then pot = float(fft(rhoft, 1, /double))

return, pe
end