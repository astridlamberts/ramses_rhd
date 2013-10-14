function pdf, arr, nbin, binctr=binctr, logctr=logctr, max=max, min=min, $
	linear=linear, binlim=binlim

; set keywords to defaults
if n_elements(max) eq 0 then max=(1.0+1.0d-10)*max(arr)
if n_elements(min) eq 0 then min=min(arr)

; set up bins
if not keyword_set(binlim) then begin
	if not keyword_set(linear) then begin
		binlim=min*(max/min)^(findgen(nbin+1)/nbin)
		pdf_return=dblarr(nbin)
		binctr=min*(max/min)^((findgen(nbin)+0.5)/nbin)
		meanval=mean(arr)
		logctr=alog(binctr/meanval)
	endif else begin
		binlim=min+(max-min)*(findgen(nbin+1)/nbin)
		pdf_return=dblarr(nbin)
		binctr=min+(max-min)*((findgen(nbin)+0.5)/nbin)
	endelse
endif else begin
	nbin=n_elements(binlim)-1
	if not keyword_set(linear) then $
		binctr=sqrt(binlim[0:nbin-1]*binlim[1:nbin]) $
	else $
		binctr=0.5*(binlim[0:nbin-1]+binlim[1:nbin])
	pdf_return=dblarr(nbin)
endelse

; compute total number of cells
ncells=n_elements(arr)

; go through the bins
for i=0,nbin-1 do begin
    masklo=arr ge binlim[i]
    maskhi=arr lt binlim[i+1]
    mask=masklo and maskhi
    pdf_return[i]=total(mask)/double(ncells)
    if (keyword_set(linear)) then $
	pdf_return[i] = pdf_return[i]/(binlim[i+1]-binlim[i]) $
    else $
	pdf_return[i] = pdf_return[i]/alog10(binlim[i+1]/binlim[i])
endfor

; return
return, pdf_return

end
