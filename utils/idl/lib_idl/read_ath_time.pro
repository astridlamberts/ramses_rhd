function read_ath_time, fname, stepnum=stepnum

; Extract time and step number from an athena rst file

; Read to end of header block
openr, fp, fname, /get_lun
dummy=''
while dummy ne 'N_STEP' do begin
	readf, fp, dummy
	if eof(fp) then begin
		message, 'Error: EOF detected'
		return, -1
	endif
endwhile

; Read step number
stepnum=0L
readu, fp, stepnum

; Read time
readf, fp, dummy
readf, fp, dummy
time=0.0d0
readu, fp, time

free_lun, fp

return, time

end

