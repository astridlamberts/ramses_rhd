function readtext, fname

; get wc of file
spawn, 'wc '+fname, res
res=strtrim(res, 1)
nline=fix(res)

if nline eq 0 then return, -1

out=dblarr(nline)

openr, fp, fname, /get_lun
readf, fp, out
free_lun, fp

return, out

end