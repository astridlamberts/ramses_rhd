function gaussmear, img, sigma

; This function takes an image and smears it with a Gaussian.

; Create the smearing kernel. Make the size the smallest power
; of 2 (in pixels) that goes out to at leat 4 sigma.
kersizemin = ceil(4*sigma)
kersize = 2^(ceil(alog(kersizemin)/alog(2)))
ker = dblarr(kersize, kersize)
for i=0, kersize-1 do begin
	ctr = 0.5*(kersize-1)
	isqr = (i - ctr)^2
	jsqr = (findgen(kersize) - ctr)^2
	ker[i,*] = exp(-(isqr+jsqr)/(2*sigma^2))
endfor
ker = ker/total(ker)

; Convolve with the kernel
return, convol(img, ker)

end
