function readdat, npart, nfile, fname

x=fltarr(3,npart)
x1=fltarr(nfile,3,npart)
for n=0, nfile-1 do begin
	if n lt 10 then begin
		openr, fp, fname+'.00'+strtrim(string(n),1), /get_lun
	endif else begin
		openr, fp, fname+'.0'+strtrim(string(n),1), /get_lun
	endelse
	readf, fp, dummy
	readf, fp, x
	free_lun, fp
	x1[n,*,*]=x
endfor

return, x1

end


