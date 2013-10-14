function autocorr, fx, binctr=binctr, nbin=nbin, nlogbin=nlogbin, $
	binlo=binlo, binhi=binhi, xmin=xmin, xmax=xmax
; The function computes the autocorrelation function of the input 3d
; array fx in bins in displacement x whose limits are set by the array
; binlimits. On return, the center of each bin is returned in the
; keyword binctr. The user may specifying the binning in the
; following ways:
;
; 1. Set nbin or nlogbin to specify the number of bins. Setting nbin
; makes the bins linear, and setting nlogbin makes them logarithmic in
; spacing. The default locations of the bin edges are at [0,1] cells
; for [linear, logarithmic] bins for the low edge, and ngrid/2 cells for
; the high edge, where ngrid minimum linear dimension of the input
; grid fx. However, these defaults may be overridden by setting xmin and
; xmax.
;
; 2. Manually specify the bin edges in the arrays binlo and binhi. In
; this case the arrays give lower and upper limits on each bin. Note
; that it must be the case that binlo[i+1] <= binhi[i] < binhi[i+1]
; for the routine to function properly, i.e. the bins must be
; increasing and non-overlapping. If binlo and binhi are not set, then
; they will be returned containing the low and high bin
; limits. Conversely, binlo and binhi will only be used if neither
; nbin nor nlogbin is set. Otherwise binlo and binhi will be ignored
; and overwritten.


; see how binning is being specified, and check for sanity
if keyword_set(nbin) and keyword_set(nlogbin) then begin
    message, 'Error: only one of nbin and nlogbin may be specified.'
    return, 0
endif
if not (keyword_set(nbin) or keyword_set(nlogbin)) then begin
    if not (keyword_set(binlo) and keyword_set(binhi)) then begin
        message, 'Error: must set nbin, nlogbin, or both binlo and binhi'
        return, 0
    endif
    if n_elements(binlo) ne n_elements(binhi) then begin
        message, 'Error: binlo and binhi must have same number of elements'
        return, 0
    endif
    numbin=n_elements(binlo)
    if total(binlo[1:nbin-1] le binhi[0:nbin-2]) ne 0 then begin
        message, 'Error: bins must be non-overlapping'
        return, 0
    endif
    if (total(binlo eq sort(binlo)) ne 0) or $
      (total(binhi eq sort(binhi)) ne 0) then begin
        message, 'Error: bins limits must be sorted in increasing order'
        return, 0
    endif
    if (n_elements(binlo) ne n_elements(uniq(binlo))) or $
      (n_elements(binhi) ne n_elements(uniq(binhi))) then begin
        message, 'Error: bin limits must be increasing'
        return, 0
    endif
endif else begin
    ngrid=min((size(fx))[1:3])
    if not keyword_set(xmax) then xmax=ngrid/2.0
    if not keyword_set(xmin) then begin
        if keyword_set(nbin) then xmin=0.0 else xmin=1.0
    endif
    if keyword_set(nbin) then begin
        numbin=nbin
        dxbin=(xmax-xmin)/numbin
        binlo=xmin+dxbin*findgen(numbin)
        binhi=xmin+dxbin*(findgen(numbin)+1)
        binctr=0.5*(binlo+binhi)
    endif else begin
        numbin=nlogbin
        dlogxbin=(xmax/xmin)^(1.0/numbin)
        binlo=xmin*dlogxbin^findgen(numbin)
        binhi=xmin*dlogxbin^(findgen(numbin)+1)
        binctr=sqrt(binlo*binhi)
    endelse
endelse

; Compute |FT(fx)|^2
wksp=fft(fx,-1)
wksp=wksp*conj(wksp)

; Now transform back to real space to get the vector autocorrelation
wksp=fft(wksp,1,/overwrite)
wksp=float(wksp)

; Now go through the grid, adding each element in the proper bin
acbin=fltarr(numbin)
cnt=lonarr(numbin)
sz=(size(wksp))[1:3]
wksp1=fltarr(sz[0], sz[1], sz[2])
for i=0,sz[0]-1 do for j=0,sz[1]-1 do $
  wksp1[i,j,*]=sqrt(float(i)^2+float(j)^2+findgen(sz[2])^2)
wksp1[0,0,0]=1.0e-6
if keyword_set(nbin) then begin
    wksp1=floor((temporary(wksp1)-xmin)/dxbin)
endif else if keyword_set(nlogbin) then begin
    wksp1=floor(alog(temporary(wksp1)/xmin)/alog(dlogxbin))
endif else begin
    wksp1=intarr(sz[0], sz[1], sz[2])
    for i=0,sz[0]-1 do begin
        for j=0,sz[1]-1 do begin
            for k=0,sz[2]-1 do begin
                wksp1[i,j,k] = where((wksp1[i,j,k] ge binlo) and $
                                     (wksp1[i,j,k] lt binhi))
            endfor
        endfor
    endfor
endelse
for n=0,numbin-1 do begin
    idx = where(wksp1 eq n)
    acbin[n] = total(wksp[idx])
    cnt[n] = n_elements(idx)
endfor

acbin=acbin/cnt

return, acbin

end
