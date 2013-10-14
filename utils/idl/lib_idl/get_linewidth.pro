function get_linewidth, sig, v, thresh=thresh

; computes the 1 sigma linewidth of a given signal. v is the velocity bins of
; each data point. If it is not passed, the linewidth is returned in
; units of velocity bins.

; Set up v
if n_elements(v) eq 0 then v = findgen(n_elements(sigsub))

; Remove regions below the signal threshhold
if keyword_set(thresh) then begin
	if total(where(sig ge thresh)) eq -1 then return, 0
	sigsub=sig[where(sig ge thresh)]
	vsub=v[where(sig ge thresh)]
	if n_elements(sigsub) lt 10 then return, 0
endif else begin
	sigsub=sig
	vsub=v
endelse

; Don't try to fit flat signals; that causes gaussfit to barf
if mean(sigsub) eq max(sigsub) then return, 0

; Now fit with a Gaussian
res = gaussfit(vsub, sigsub, gfit, nterms=3)

; return linewidth
return, gfit[2]

end

