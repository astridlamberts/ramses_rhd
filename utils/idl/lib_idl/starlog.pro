; get wc of file
spawn, 'wc SinkLog', res
res=strtrim(res, 1)
nline=reform(long(res))

; allocate storage arrays
mstar=dblarr(nline)
rstar=dblarr(nline)
time=dblarr(nline)
burnstate=intarr(nline)
mdotstar=dblarr(nline)
mdeutstar=dblarr(nline)

; open file
openr, fp, 'SinkLog', /get_lun

; read data
strline=''
ddummy=0.0d0
idummy=0
for n=0L, nline[0]-1 do begin

    ; note: fortran IO is sufficiently bad that it is easier just to
    ; read the data as a string and do my own parsing
    readf, fp, strline

    ; break up the line using spaces as delimiters
    fields=strsplit(strline, /extract)

    ; store the data
    reads, fields[0], ddummy
    time[n] = ddummy
    reads, fields[2], ddummy
    mstar[n] = ddummy
    reads, fields[13], ddummy
    rstar[n] = ddummy
    reads, fields[14], ddummy
    mdeutstar[n] = ddummy
    reads, fields[16], ddummy
    mdotstar[n] = ddummy
    reads, fields[17], idummy
    burnstate[n] = idummy

endfor

; close file
free_lun, fp

end
