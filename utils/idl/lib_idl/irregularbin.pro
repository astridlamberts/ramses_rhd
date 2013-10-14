function irregularBin, x, fx, binlimits, noaverage=noaverage
; given a function fx evaluated at a set of points x, compute the
; average of the function in a series of bins specified by the array
; binlimits, which gives the x values for the bin edges. binlimits
; must be sorted into ascending order.

nx=n_elements(x)
nbin=n_elements(binlimits)-1
binval=fltarr(nbin)
numbin=fltarr(nbin)
for i=0, nx-1 do begin
	binidx=max(where(binlimits le x[i]))
	if (binidx ge 0) and (binidx lt nbin) then begin
		binval[binidx]=binval[binidx]+fx[i]
		numbin[binidx]=numbin[binidx]+1
	endif
endfor

if not keyword_set(noaverage) then begin
	numbin=numbin+1*(numbin eq 0)	; avoid divide by zero
	binval=binval/numbin
endif

return, binval
end
