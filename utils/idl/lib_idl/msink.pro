; get wc of files
spawn, 'wc msink.txt', res
res=strtrim(res, 1)
nline=fix(res)

msink=dblarr(nline)
time=dblarr(nline)

openr, fp, 'msink.txt', /get_lun
readf, fp, msink
free_lun, fp

openr, fp, 'time.txt', /get_lun
readf, fp, time
free_lun, fp

mdot=(msink[1:nline-1]-msink[0:nline-2])/(time[1:nline-1]-time[0:nline-2])
thalf=(time[1:nline-1]+time[0:nline-2])/2.0

plot, time, msink-msink[0], psym=-2

print, linfit(time, msink)

end