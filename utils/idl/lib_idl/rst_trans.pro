pro rst_trans, fname, outname, dtout, nprad=nprad, raddir=raddir,$ 
radint=radint, parallel=parallel, seed=seed, swap=swap, ib_x1=ib_x1,$
         ib_x2=ib_x2,ib_x3=ib_x3,ob_x1=ob_x1, ob_x2=ob_x2,ob_x3=ob_x3,CFL=CFL

; converts a restart file to a new one with just one or 2 different feautures. ?;o rescaling is done, just changing some values in the restart file

; get new velocity scaling
;print,temp

;kb=1.38d-16
;mu=3.9d-24
;gamma=5./3.
;cs=double(sqrt(kb*temp/mu))
;print,cs
; random number generator stuff
;IM1=2147483563LL
;IA1=40014LL
;IQ1=53668LL
;IR1=12211LL
;NTAB_RAN2=32
;if not keyword_set(seed) then idum=2674LL else idum=long64(seed)
;idum2=idum
;iy=0LL
;iv=lon64arr(NTAB_RAN2)
;for j=NTAB_RAN2+7,0,-1 do begin
;    k=idum/IQ1
;    idum=IA1*(idum-k*IQ1)-k*IR1
;    if idum lt 0 then idum += IM1
;    if j lt NTAB_RAN2 then iv[j]=idum
;endfor
;iy=iv[0]

; get srad
;if not keyword_set(srad) and keyword_set(rs) then begin
 ;   alpha4=2.59d-13
;    mu=2.34d-24
;    srad=4.d/3.d*!pi*double(rs)^3*alpha4*(double(dbar)/mu)^2
;endif

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

                                ; construct file name
    if n ne 0 then begin
        fnamesplit=strsplit(fname, '.', /extract)
        nsplit=n_elements(fnamesplit)
        fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
        newname=''
        for i=0,nsplit-2 do newname+=fnamesplit[i]+'.'
        newname+=fnamesplit[nsplit-1]
        fnamesplit=strsplit(outname, '.', /extract)
        nsplit=n_elements(fnamesplit)
        fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
        newoutname=''
        for i=0,nsplit-2 do newoutname+=fnamesplit[i]+'.'
        newoutname+=fnamesplit[nsplit-1]
    endif else begin
        newname=fname
        newoutname=outname
    endelse

                                ; open input and output files
    openr, fp, newname, /get_lun
    openw, fpout, newoutname, /get_lun
    print, 'Converting '+newname+' to '+newoutname+'...'

                                ; parse header
    outputblock=0
    while dummy ne '<par_end>' do begin
        readf, fp, dummy
        if strmid(dummy,0,2) eq 'Nx' then begin
            dim=strmid(dummy,2,1)
            pos=strpos(dummy,'=')+1
            if dim eq 1 then nx=fix(strmid(dummy,pos)) $
            else if dim eq 2 then ny=fix(strmid(dummy,pos)) $
                 else nz=fix(strmid(dummy,pos))
            outstr=dummy
       ; endif else if strmid(dummy,0,5) eq 'x1min' then begin
       ;     pos=strpos(dummy,'=')
       ;     outstr=strmid(dummy,0,pos+1)+string(-lbox/2)
       ; endif else if strmid(dummy,0,5) eq 'x2min' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(-lbox/2)
       ; endif else if strmid(dummy,0,5) eq 'x3min' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(-lbox/2)
       ; endif else if strmid(dummy,0,5) eq 'x1max' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(lbox/2)
       ; endif else if strmid(dummy,0,5) eq 'x2max' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(lbox/2)
       ; endif else if strmid(dummy,0,5) eq 'x3max' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(lbox/2)
       ; endif else if strmid(dummy,0,4) eq 'tlim' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+'1.0e15'
        endif else if strmid(dummy,0,7) eq 'cour_no' then begin 
            pos=strpos(dummy,'=')+1
            outstr=strmid(dummy,0,pos+1)+string(CFL)

       ; endif else if strmid(dummy,0,6) eq 'ibc_x1' then begin
       ;     pos=strpos(dummy,'=')+1
       ;     outstr=strmid(dummy,0,pos+1)+string(ib_x1)
   ; endif else if strmid(dummy,0,6) eq 'ibc_x2' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+string(ib_x2)
   ; endif else if strmid(dummy,0,6) eq 'ibc_x3' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+string(ib_x3)
   ; endif else if strmid(dummy,0,6) eq 'obc_x1' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+string(ob_x1)
   ; endif else if strmid(dummy,0,6) eq 'obc_x2' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+string(ob_x2)
   ;       endif else if strmid(dummy,0,6) eq 'obc_x3' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+string(ob_x3)

  

     ;   endif else if strmid(dummy,0,1) eq '<' then begin
     ;       if strmid(dummy,0,7) eq '<output' then outputblock=1 $
     ;       else outputblock=0
            if strmid(dummy,0,9) eq '<par_end>' then begin
                printf, fpout, '<ionradiation>'
                printf, fpout, 'ray_number        = 2        # ray number'
                printf, fpout, 'min_tree_level    = 2        # minimum level to use in ray tree'
                printf, fpout, 'sigma_ph          = 6.3e-18  # ionization cross section at threshhold'
                printf, fpout, 'm_H               = 2.34e-24 # gas mass per H atom'
                printf, fpout, 'mu                = 3.90e-24 # mean mass per particle in neutral gas'
                printf, fpout, 'e_gamma           = 3.84e-12 # energy deposited per photoionization (= 2.4 eV)'
                printf, fpout, 'alpha_C           = 3.0e-3   # carbon abundance'
                printf, fpout, "k_B               = 1.38e-16 # Boltzmann's constant"
                printf, fpout, 'time_unit         = 1.0      # number of seconds per unit of code time'
                printf, fpout, 'max_de_iter       = 0.5      # maximum change in total energy due to'
                printf, fpout, 'max_de_therm_iter = 0.5      # maximum change in thermal energy due to'
                printf, fpout, 'max_dx_iter       = 0.5      # maximum change in ionization fraction'
                printf, fpout, 'max_de_step       = 10       # maximum change in total energy due to'
                printf, fpout, 'max_de_therm_step = 10       # maximum change in thermal energy due to'
                printf, fpout, 'max_dx_step       = -1       # maximum change in ionization fraction'
                printf, fpout, 'tfloor            = 10.0     # temperature floor for radiation'
                printf, fpout, 'tceil             = 1.0e5    # temperature for radiation'
                printf, fpout, 'maxiter           = 5000     # maximum number of radiation iterations'
                printf, fpout, 'rebuild_interval  = 3        # frequency of ray tree rebuild'
                printf, fpout, ''
            endif
            outstr=dummy
    ;    endif else if strmid(dummy,0,3) eq 'num' then begin
    ;        pos=strpos(dummy,'=')+1
    ;        outstr=strmid(dummy,0,pos+1)+'1'
    ;    endif else if strmid(dummy,0,5) eq 'nstep' then begin
    ;        pos=strpos(dummy,'=')+1
    ;        outstr=strmid(dummy,0,pos+1)+'0'
        endif else if strmid(dummy,0,3) eq 'dt ' then begin
            pos=strpos(dummy,'=')+1
          ;  if outputblock eq 1 then begin
                outstr=strmid(dummy,0,pos+1)+string(dtout)
          ;  endif else begin
          ;      dt=double(strmid(dummy,pos+1))
          ;      outstr=strmid(dummy,0,pos+1)+string(dt*cellsz/cs)
          ;  endelse
        endif else if strmid(dummy,0,4) eq 'time' then begin
            pos=strpos(dummy,'=')+1
            if outputblock eq 1 then begin
                outstr=strmid(dummy,0,pos+1)+string(dtout)
            endif else begin
                outstr=strmid(dummy,0,pos+1)+'0.00000'
            endelse
      ;  endif else if strmid(dummy,0,10) eq 'problem_id' then begin
      ;      pos=strpos(dummy,'=')+1
      ;      fnamesplit=strsplit(outname, '.', /extract)
      ;      nsplit=n_elements(fnamesplit)
      ;      outstr=strmid(dummy,0,pos+1)+fnamesplit[nsplit-3]
      ;  endif else if strmid(dummy,0,7) eq 'problem' then begin
      ;      pos=strpos(dummy,'=')+1
      ;      fnamesplit=strsplit(outname, '.', /extract)
      ;      nsplit=n_elements(fnamesplit)
      ;      outstr=strmid(dummy,0,pos+1)+fnamesplit[nsplit-3]
      ;  endif else if strmid(dummy,0,12) eq 'restart_file' then begin
      ;      pos=strpos(dummy,'=')+1
      ;      fnamesplit=strsplit(outname, '.', /extract)
      ;      nsplit=n_elements(fnamesplit)
      ;      outstr=strmid(dummy,0,pos+1)+fnamesplit[nsplit-3]
   ;     endif else if strmid(dummy,0,18) eq 'ionizing radiation' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+'yes'
   ;     endif else if strmid(dummy,0,8) eq 'eq_state' then begin
   ;         pos=strpos(dummy,'=')+1
   ;         outstr=strmid(dummy,0,pos+1)+'adiabatic'
   ;     endif else if strmid(dummy,0,10) eq 'iso_csound' then begin
   ;         outstr='gamma   = 1.666667     # gamma = c_p / c_v'
   ;      endif else begin
   ;         outstr=dummy
        endif  ;endelse
   ;     printf, fpout, outstr
    endwhile
   ; cellsz=double(lbox/max([nx,ny,nz]))

                                ; set time step number
   ; readf, fp, dummy
   ; printf, fpout, dummy
   ; nstep=0L
   ; readu, fp, nstep
   ; writeu, fpout, 0L
   ; readf, fp, dummy
   ; printf, fpout, dummy

                                ; set time to 0 and scale dt based on
                                ; new size and velocity scaling
   ; readf, fp, dummy
   ; printf, fpout, dummy
   ; time=0.0d0
   ; readu, fp, time
   ; writeu, fpout, 0.0d0
   ; readf, fp, dummy
   ; printf, fpout, dummy
   ; readf, fp, dummy
   ; printf, fpout, dummy
   ; dt=0.0d0
   ; readu, fp, dt
   ; writeu, fpout, dt*cellsz/cs

                                ; set grid size for parallel runs
   ; if n_elements(parallel) eq 3 then begin
   ;     nx = nx/parallel[0]
   ;     ny = ny/parallel[1]
   ;     nz = nz/parallel[2]
   ; endif

                                ; allocate memory and read data
    ;d =dblarr(nx,ny,nz)
    ;px=dblarr(nx,ny,nz)
    ;py=dblarr(nx,ny,nz)
    ;pz=dblarr(nx,ny,nz)
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readu, fp, d
    ;d=double(d*dbar)
    ;writeu, fpout, d
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readu, fp, px
    ;px=double(px*dbar*cs)
    ;writeu, fpout, px
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readu, fp, py
    ;py=double(py*dbar*cs)
    ;writeu, fpout, py
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readf, fp, dummy
    ;printf, fpout, dummy
    ;readu, fp, pz
    ;pz=double(pz*dbar*cs)
    ;writeu, fpout, pz
    ;if keyword_set(mhd) then begin
     ;   bx=dblarr(nx+1,ny,nz)
     ;   by=dblarr(nx,ny+1,nz)
     ;   bz=dblarr(nx,ny,nz+1)
     ;   readf, fp, dummy1
     ;   readf, fp, dummy2
     ;   readu, fp, bx
     ;   bx=double(bx*sqrt(dbar)*cs)
     ;   readf, fp, dummy3
     ;   readf, fp, dummy4
     ;   readu, fp, by
     ;   by=double(by*sqrt(dbar)*cs)
     ;   readf, fp, dummy5
     ;   readf, fp, dummy6
     ;   readu, fp, bz
     ;   bz=double(bz*sqrt(dbar)*cs)
   ; endif
   ; printf, fpout, ''
   ; printf, fpout, 'ENERGY'
   ; e=double(0.5*(px^2+py^2+pz^2)/d + 0.5*(bx^2+by^2+bz^2) + $
   ;          d*kb*temp/(mu*(gamma-1)))
   ; writeu, fpout, e
   ; printf, fpout, dummy1
   ; printf, fpout, dummy2
   ; writeu, fpout, bx
   ; printf, fpout, dummy3
   ; printf, fpout, dummy4
   ; writeu, fpout, by
   ; printf, fpout, dummy5
   ; printf, fpout, dummy6
   ; writeu, fpout, bz
   ; printf, fpout, ''
   ; printf, fpout, 'RADIATOR PLANE LIST'
   ; writeu, fpout, long(nprad)  
   ; for m=0,nprad-1 do begin      
   ;         writeu, fpout, long(raddir[m])
   ;          writeu, fpout, double(radint[m])
   ; endfor
   ; printf, fpout, ''
   ; printf, fpout, 'SCALAR 0'
   ; writeu, fpout, d
;  ;  printf, fpout, ''
;    printf, fpout, 'USER_DATA'

; transfer user data
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    onebyte = bytarr(1)
    while (not eof(fp)) do begin
      readu, fp, onebyte
      writeu, fpout, onebyte
    endwhile
    
                                ; close files
    free_lun, fp
    free_lun, fpout

endfor

end
