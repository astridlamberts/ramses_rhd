pro rst_transformation,fname,dtout=dtout,cfl=cfl,parallel=parallel

; some dummy variables
dummy=''
dummy1=''
dummy2=''
dummy3=''
dummy4=''
dummy5=''
dummy6=''

; handle parallelism
if n_elements(parallel) eq 3 then $
  nproc=parallel[0]*parallel[1]*parallel[2] $
else nproc=1
for n=0,nproc-1 do begin
;print,'parallelok'
                            ; construct file name
    if n ne 0 then begin
        fnamesplit=strsplit(fname, '.', /extract)
print,'a'
        nsplit=n_elements(fnamesplit)
print,'b'
        fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
print,'c'
        newname=''
print,'d'      
  for i=0,nsplit-2 do newname+=fnamesplit[i]+'.'
print,'e'      
  newname+=fnamesplit[nsplit-1]
print,'f' 
    ;   fnamesplit=strsplit(outname, '.', /extract)
;print,'g'       
; nsplit=n_elements(fnamesplit)
print,'h'       
; fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
     ;   newoutname=''
     ;   for i=0,nsplit-2 do newoutname+=fnamesplit[i]+'.'
     ;   newoutname+=fnamesplit[nsplit-1]
    endif else begin
        newname=fname
      ;  newoutname=outname
    endelse

 ; open input file
    openu, fp, newname, /get_lun  ;opens the file for update, get_lun lets IDL assign a free logical number between 100 and 128 to each file
;    openw, fpout, newoutname, /get_lun
    print, 'Converting '+newname+' 

;changing the file
 while dummy ne '<par_end>' do begin

  readf, fp, dummy ;reads a variable in the file labelled 'fp' and will asign it to 'dummy'
            if strmid(dummy,0,7) eq 'cour_no' then begin 
            pos=strpos(dummy,'=')+1 ;strpos finds the first occurence of a substring in a string
            outstr=strmid(dummy,0,pos+1)+string(CFL)

            endif else if strmid(dummy,0,3) eq 'dt ' then begin
            pos=strpos(dummy,'=')+1
            outstr=strmid(dummy,0,pos+1)+string(dtout)
         
            endif else if strmid(dummy,0,4) eq 'time' then begin
            pos=strpos(dummy,'=')+1
            outstr=strmid(dummy,0,pos+1)+string(dtout)
            endif
 endwhile
; close files
    free_lun, fp

endfor
end

