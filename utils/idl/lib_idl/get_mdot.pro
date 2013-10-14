function get_mdot, time=time, msink=msink

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

res=linfit(time, msink)
return, res[1]

end