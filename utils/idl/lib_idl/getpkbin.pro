function getPkBin, fx, binlimits, nbin=nbin, nlogbin=nlogbin, $
	binctr=binctr, log=log, power=power
; The function computes the power spectral density (or power, if power
; is set) of the input 3d array fx. The PSD is computed in bins whose
; limits are set by the array binlimits. If nbin is set, the data are
; automatically binned into nbin linear-spaced bins, and binlimits is
; returned with the limits on these bins. nlogbin does the same thing,
; but uses logarithmically spaced bins. The bin centers are returned
; in the keyword binctr. If neither nbin nor nlogbin is set, and log
; is not set, the centers are computed linearly. If neither nbin nor
; nlogbin is set and log is set, they are computed logarithmically.

; First set up the bins
sz=n_elements(fx[0,0,*])
if keyword_set(nbin) and keyword_set(nlogbin) then begin
	print, 'Error: only one of nbin and nlogbin may be specified.'
	return, 0
endif
if keyword_set(nbin) then begin
	kmax=sz*sqrt(3.0)/2.0
	binlimits=findgen(nbin+1)/nbin * kmax
	binctr=0.5*(binlimits[0:nbin-1]+binlimits[1:nbin])
endif
if keyword_set(nlogbin) then begin
	kmax=sz*sqrt(3.0)/2.0
	binlimits=exp(findgen(nlogbin+1)*alog(kmax)/nlogbin)
	binctr=sqrt(binlimits[0:nlogbin-1]*binlimits[1:nlogbin])
endif
numbin=n_elements(binlimits)-1
if not keyword_set(nlogbin) and not keyword_set(nbin) then begin
	if keyword_set(log) then binctr = $
		0.5*(binlimits[0:numbin-1]+binlimits[1:numbin]) $
	else binctr = $
		sqrt(binlimits[0:numbin-1]*binlimits[1:numbin])
endif

; Compute |FT(fx)|^2
fxtemp=complex(fx)
fxft=fft(fxtemp,-1)
fxft=float(fxft*conj(fxft))*sz^3

; Now go through the grid, putting each element in the proper bin
pk=fltarr(numbin)
binlimitssqr=binlimits^2
for i=0L,sz/2 do begin
   for j=0L,sz/2 do begin
      for k=0L,sz/2 do begin
	bin=max(where(binlimitssqr le (i^2+j^2+k^2)))
	if (bin ge 0) and (bin lt numbin) then begin
		pk[bin] = pk[bin] + fxft[i,j,k]
		if k ne 0 then pk[bin] = pk[bin] + fxft[i,j,sz-k]
		if j ne 0 then pk[bin] = pk[bin] + fxft[i,sz-j,k]
		if i ne 0 then pk[bin] = pk[bin] + fxft[sz-i,j,k]
		if (k ne 0) and (j ne 0) then $
			pk[bin] = pk[bin] + fxft[i,sz-j,sz-k]
		if (k ne 0) and (i ne 0) then $
			pk[bin] = pk[bin] + fxft[sz-i,j,sz-k]
		if (j ne 0) and (i ne 0) then $
			pk[bin] = pk[bin] + fxft[sz-i,sz-j,k]
		if (k ne 0) and (j ne 0) and (i ne 0) then $
			pk[bin] = pk[bin] + fxft[sz-i,sz-j,sz-k]
	endif
      endfor
   endfor
endfor

; if power is unset, divide by the bin volume to give PSD rather than power
if not keyword_set(power) then pk = pk / ((4.0/3.0)*!pi* $
	(binlimits[1:numbin]^3-binlimits[0:numbin-1]^3))

; return
return, pk
end


pro getPk, fx, kvec, pk, totpower=totpower
; This function computes the power spectral density pk of a function
; fx. kvec is a list of |k| values for which pk is tabulated. If the
; keyword totpower is set, rather than computing the power spectral
; density the routine computes the total power at frequency |k|. The
; difference is that the power spectral density is the total power
; normalized by the k-space volume element at frequency |k|.

; First take fft of fx
fxtemp=complex(fx)
fxft=fft(fxtemp,-1)
sz=n_elements(fx[0,0,*])

; Square f~(\vec{k}) and discard complex parts
fxft=float(fxft*conj(fxft))

; Figure out all the values of k and store to kvec
distarr=lonarr(sz^3, /nozero)
fxft1d=fltarr(sz^3, /nozero)
for i=0, sz-1 do begin
	for j=0, sz-1 do begin
		for k=0, sz-1 do begin
			distarr[i*sz^2+j*sz+k] = (i < (sz-i))^2 + $
				(j < (sz-j))^2 + (k < (sz-k))^2
			fxft1d[i*sz^2+j*sz+k] = fxft[i,j,k]
		endfor
	endfor
endfor
sortlist=sort(distarr)
kvecsqr=distarr[sortlist]

; get pk, averaging values of fxft at the same k
pk=fltarr(n_elements(distarr), /nozero)
pkptr=0L
pk[0]=fxft1d[sortlist[0]]
for n=1L, n_elements(distarr)-1 do begin
	if kvecsqr[n] eq kvecsqr[n-1] then begin
		pk[pkptr]=pk[pkptr]+fxft1d[sortlist[n]]
						; add to pk
	endif else begin
		pkptr=pkptr+1			; move pointer
		pk[pkptr]=fxft1d[sortlist[n]]	; store new pk
		kvecsqr[pkptr]=kvecsqr[n]	; store new k^2
	endelse
endfor

; extract subarrays that contain useful data and set the correct k
pk=pk[0:pkptr]
kvec=sqrt(kvecsqr[0:pkptr])*sqrt(3.0)/(2.0*sqrt(kvecsqr[pkptr]))

; we now have the total power. If we want the power spectral density,
; we now have to divide by 4*pi*k^2.
if not keyword_set(totpower) then $
	pk[1:pkptr]=pk[1:pkptr]/(4*3.1415927*kvec[1:pkptr]^2)

end


function spectralReshape, fx, idx, totpower=totpower
; given a function fx defined as an [n,n,n] array, return a new
; function on an [n,n,n] array that has the same total power and
; Fourier-domain phases as fx, but has a power spectrum that fits
; k^-idx. There is a subtelty in doing this: because the fft exists
; on a discretely sampled grid, there may be a varying number of
; points on the grid corresponding to a given physical k. For example,
; if a [n,n,n] element array is passed, there is only one point
; corresponding to the k=n*sqrt(3), while there are 3 points
; corresponding to k=1. If one doesn't correct for this effect, then
; the function returned will have a total power that goes as k^index,
; rather than a power-spectral-density that goes as k^index. If the
; keyword totpower is set, the routine produces a total power
; k^index. If not, it produces a power-spectral-density k^index.

sz=long(n_elements(fx[0,0,*]))

; Compute total power in fx
pow=total(fx^2)

; Take fft of fx
fxtemp=complex(fx)
fxnewft=fft(fxtemp,-1)

; compute weights
weight=fltarr(sz,sz,sz)		; array to hold weights
m=lindgen(sz)
m = m < (sz-m)
for l=0L, sz/2 do begin
	for n=0L, sz/2 do begin

		ksqr = m^2 + l^2 + n^2		; k value

		; compute row of weights
		ksqr1=ksqr+(ksqr eq 0)
			; avoid divide by zero
		weightrow = ksqr1^(idx/4.0)
		weightrow = weightrow*(ksqr ne 0)
			; eliminate 0 frequency

		; do normalization if totpower isn't set
		if not keyword_set(totpower) then begin
			weightrow=weightrow*sqrt(ksqr)
		endif

		; set weight in all 8 octants, 2 at a time
		weight[l,n,*] = weightrow
		if n ne 0 then weight[l,sz-n,*] = weightrow
		if l ne 0 then weight[sz-l,n,*] = weightrow
		if (n ne 0) and (l ne 0) then $
			weight[sz-l,sz-n,*] = weightrow
	endfor
endfor

; compute multiplicity factors
; There are likely to be several distinct gridpoints at any given
; value of k. There is no easy way to predict how many gridpoints fall
; on a given k-sphere. Therefore we have to generate all the k's that
; occur in our grid and count the number of repetitions of each to get
; a multiplicity for that k. 

mult=intarr(sz,sz,sz)		; array to hold multiplicities
dist1d=intarr(sz^3)
n=lindgen(sz)
nsqr = (n < (sz-n))^2
for l=0L, sz-1 do begin		; generate all distances
	for m=0L, sz-1 do begin
		dist1d[(l*sz^2+m*sz) : (l*sz^2+m*sz+sz-1)] = $
			(l < (sz-l))^2 + (m < (sz-m))^2 + nsqr
	endfor
endfor
dist1d=dist1d[sort(dist1d)]	; sort distances
uniqlist=uniq(dist1d)		; find unique distances
dist1d=dist1d[uniqlist]		; keep only unique distances
mult1d=[1, uniqlist[1:n_elements(uniqlist)-1] - $
	uniqlist[0:n_elements(uniqlist)-2]]
				; compute multiplicities
for l=0L, sz-1 do begin		; put multiplicities into matrix
	for m=0L, sz-1 do begin
		for n=0L, sz-1 do begin
			dist = (l < (sz-l))^2 + (m < (sz-m))^2 + $
				(n < (sz-n))^2
			idx = where(dist eq dist1d)
				; find index in distance list
			mult[l,m,n] = mult1d[idx]
				; store multiplicity
		endfor
	endfor
endfor

; apply multiplicity correction
weight=weight/sqrt(mult)

; Assign weights to new function
fxnewft=fxnewft/abs(fxnewft)*weight

; Go back to space domain
fnew=fft(fxnewft,1)
fnew=float(fnew)

; Renormalize fnew to keep power constant
pownew=total(fnew^2)
fnew=fnew*sqrt(pow/pownew)

return, fnew
end
