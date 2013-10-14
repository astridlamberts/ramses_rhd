function beamsmear, img, dx, dist, beamsize, sigmapix=sigmapix

; This function takes an image and applies a beam to it. The image is
; a 2d array passed in as img, with pixels of size dx in physical space
; (in cm). The distance to the object (in kpc) is specified by dist.
; The beamsize (in arcsec) is the FWHM of the beam. The function returns
; a smeared image of the same size as the original.

; First compute the size of the beam in pixels
if not keyword_set(sigmapix) then begin
	kpc = 3.09d21
	fwhmpix = (2*!pi*dist*kpc * beamsize/(360.*60.*60.)) / dx
	sigmapix = fwhmpix/(2.*sqrt(2.*alog(2.)))
endif

; Now create the smearing kernel. Make the size the smallest power
; of 2 (in pixels) that goes out to at leat 4 sigma.
kersizemin = ceil(4*sigmapix)
kersize = 2^(ceil(alog(kersizemin)/alog(2)))
ker = dblarr(kersize, kersize)
for i=0, kersize-1 do begin
	ctr = 0.5*(kersize-1)
	isqr = (i - ctr)^2
	jsqr = (findgen(kersize) - ctr)^2
	ker[i,*] = exp(-(isqr+jsqr)/(2*sigmapix^2))
endfor
ker = ker/total(ker)

; Convolve with the kernel
return, convol(img, ker)

end
