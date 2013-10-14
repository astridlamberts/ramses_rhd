openr, fp, 'SinkLog.bin', /get_lun

npart=0L
nline=0L

readu, fp, nline
readu, fp, npart

data=dblarr(npart+1, nline)

readu, fp, data

free_lun, fp

time=reform(data[0,*])
msink=data[1:npart,*]

end
