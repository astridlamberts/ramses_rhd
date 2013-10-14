pro gridshift_ath, fname, outname, shift=shift, parallel=parallel, dmax=dmax, $
                   isothermal=isothermal, mhd=mhd, ion=ion, nscalar=nscalar

; some dummy variables
dummy=''
dummy1=''
dummy2=''
dummy3=''
dummy4=''
dummy5=''
dummy6=''


; read data

; open file
openr, fp, fname, /get_lun

; read parameters block to get problem size and limits
dummy=''
xmin=dblarr(3)
xmax=dblarr(3)
while dummy ne '<par_end>' do begin
    readf, fp, dummy
    if strmid(dummy,0,2) eq 'Nx' then begin
        dim=strmid(dummy,2,1)
        pos=strpos(dummy,'=')+1
        if dim eq 1 then nx=fix(strmid(dummy,pos)) $
        else if dim eq 2 then ny=fix(strmid(dummy,pos)) $
        else nz=fix(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x1min' then begin
        pos=strpos(dummy,'=')+1
        xmin[0]=double(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x2min' then begin
        pos=strpos(dummy,'=')+1
        xmin[1]=double(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x3min' then begin
        pos=strpos(dummy,'=')+1
        xmin[2]=double(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x1max' then begin
        pos=strpos(dummy,'=')+1
        xmax[0]=double(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x2max' then begin
        pos=strpos(dummy,'=')+1
        xmax[1]=double(strmid(dummy,pos))
    endif else if strmid(dummy,0,5) eq 'x3max' then begin
        pos=strpos(dummy,'=')+1
        xmax[2]=double(strmid(dummy,pos))
    endif
endwhile
xvec=(findgen(nx)+0.5)*(xmax[0]-xmin[0]-1)/nx-xmin[0]
yvec=(findgen(ny)+0.5)*(xmax[1]-xmin[1]-1)/ny-xmin[1]
zvec=(findgen(nz)+0.5)*(xmax[2]-xmin[2]-1)/nz-xmin[2]

; read time step number
readf, fp, dummy
nstep=0L
readu, fp, nstep
readf, fp, dummy

; read time and dt
readf, fp, dummy
time=0.0d0
readu, fp, time
readf, fp, dummy
readf, fp, dummy
dt=0.0d0
readu, fp, dt

; get number of quantities
if keyword_set(nscalar) then nadvect=nscalar else nadvect=0
if keyword_set(ion) then nadvect=nadvect+1
nvar = 4
if not keyword_set(isothermal) then nvar = nvar + 1
if keyword_set(mhd) then nvar = nvar + 3
nvar = nvar+nadvect

; check if this is a parallel run
if n_elements(parallel) eq 3 then begin
    nxgrid = nx/parallel[0]
    nygrid = ny/parallel[1]
    nzgrid = nz/parallel[2]
    nxb=nxgrid
    nyb=nygrid
    nzb=nzgrid
    block = dblarr(nxgrid,nygrid,nzgrid)
    nproc = parallel[0]*parallel[1]*parallel[2]
endif else begin
    nxb = nx
    nyb = ny
    nzb = nz
endelse

; allocate memory and read data
d =dblarr(nx,ny,nz)
px=dblarr(nx,ny,nz)
py=dblarr(nx,ny,nz)
pz=dblarr(nx,ny,nz)
if n_elements(parallel) ne 3 then begin
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, d
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, px
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, py
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, pz
    if not keyword_set(isothermal) then begin
        readf, fp, dummy
        readf, fp, dummy
        e=dblarr(nx,ny,nz)
        readu, fp, e
    endif
    if keyword_set(mhd) then begin
        bx=dblarr(nx+1,ny,nz)
        by=dblarr(nx,ny+1,nz)
        bz=dblarr(nx,ny,nz+1)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, bx
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, by
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, bz
    endif
    if keyword_set(ion) then begin
        readf, fp, dummy
        readf, fp, dummy
        d_n=dblarr(nx,ny,nz)
        readu, fp, d_n
    endif
    if keyword_set(nscalar) then begin
        s=dblarr(nscalar,nx,ny,nz)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, s
    endif
endif else begin
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, block
    d[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, block
    px[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, block
    py[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, block
    pz[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    if not keyword_set(isothermal) then begin
        e=dblarr(nx,ny,nz)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        e[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    endif
    if keyword_set(mhd) then begin
	xfblock = dblarr(nxgrid+1,nygrid,nzgrid)
	yfblock = dblarr(nxgrid,nygrid+1,nzgrid)
	zfblock = dblarr(nxgrid,nygrid,nzgrid+1)
        bx=dblarr(nx+1,ny,nz)
        by=dblarr(nx,ny+1,nz)
        bz=dblarr(nx,ny,nz+1)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, xfblock
        bx[0:nxgrid,0:nygrid-1,0:nzgrid-1]=xfblock
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, yfblock
        by[0:nxgrid-1,0:nygrid,0:nzgrid-1]=yfblock
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, zfblock
        bz[0:nxgrid-1,0:nygrid-1,0:nzgrid]=zfblock
    endif
    if keyword_set(ion) then begin
        d_n=dblarr(nx,ny,nz)
        dummyd=0.0d0
        dummyf=0.0
        dummyl=0L
        nrad=0L
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, nrad
        xrad=dblarr(nrad, 3)
        srad=dblarr(nrad)
        for m=1,nrad do begin
            for l=1,3 do begin
                readu, fp, dummyd
                xrad[m-1,l-1]=dummyd
            endfor
            readu, fp, dummyd
            srad[m-1]=dummyd
            readu, fp, dummyl
            for l=1,9 do readu, fp, dummyf
        endfor
        readf, fp, dummy
        readf, fp, dummy
        for m=1,70 do readu, fp, dummyl
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        d_n[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
    endif
    if keyword_set(nscalar) then begin
        s=dblarr(nscalar,nx,ny,nz)
        sblock=dblarr(nscalar,nxgrid,nygrid,nzgrid)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, sblock
        s[*,0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=sblock
    endif
endelse

; close file
free_lun, fp

; if parallel, read other files
if n_elements(parallel) eq 3 then begin
    xgrid=0
    ygrid=0
    zgrid=0
    for n=1,nproc-1 do begin

        ; construct file name
        fnamesplit=strsplit(fname, '.', /extract)
        nsplit=n_elements(fnamesplit)
        fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
        newname=''
        for i=0,nsplit-2 do newname+=fnamesplit[i]+'.'
        newname+=fnamesplit[nsplit-1]

        ; figure out where this block goes
        xgrid++
        if xgrid eq parallel[0] then begin
            xgrid=0
            ygrid++
            if ygrid eq parallel[1] then begin
                ygrid=0
                zgrid++
            endif
        endif
        xlo=nxgrid*xgrid
        ylo=nygrid*ygrid
        zlo=nzgrid*zgrid

        ; open and read header
        openr, fp, newname, /get_lun
        while dummy ne '<par_end>' do readf, fp, dummy

        ; read data
        readf, fp, dummy
        readu, fp, nstep
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, time
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, dt
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        d[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        px[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        py[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, block
        pz[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        if not keyword_set(isothermal) then begin
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, block
            e[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        endif
        if keyword_set(mhd) then begin
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, xfblock
            bx[xlo:xlo+nxgrid,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=xfblock
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, yfblock
            by[xlo:xlo+nxgrid-1,ylo:ylo+nygrid,zlo:zlo+nzgrid-1]=yfblock
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, zfblock
            bz[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid]=zfblock
        endif
        if keyword_set(ion) then begin
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, nrad
            for m=1,nrad do begin
                for l=1,4 do readu, fp, dummyd
                readu, fp, dummyl
                for l=1,9 do readu, fp, dummyf
            endfor
            readf, fp, dummy
            readf, fp, dummy
            for m=1,70 do readu, fp, dummyl
            readf, fp, dummy
            readf, fp, dummy
            readu, fp, block
            d_n[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
        endif
        if keyword_set(nscalar) then begin
            readf, fp, dummy
            readf, fp, dummy
            sblock=dblarr(nscalar,nxgrid,nygrid,nzgrid)
            readu, fp, sblock
            s[*xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=sblock
        endif

        ; close file
        free_lun, fp

    endfor
endif


; get shift if maxd is set, and shift velocities too
if keyword_set(dmax) then begin
    pos=findmaxnd(d)
    shift=[nx,ny,nz]/2-pos
    v=[px[pos[0],pos[1],pos[2]], py[pos[0],pos[1],pos[2]], $
       pz[pos[0],pos[1],pos[2]]] / d[pos[0],pos[1],pos[2]]
    px=d*(px/d-v[0])
    py=d*(py/d-v[1])
    pz=d*(pz/d-v[2])
endif

; now we have the data stored globally, so shift it
d=shift(d, shift)
px=shift(px, shift)
py=shift(py, shift)
pz=shift(pz, shift)
if not keyword_set(isothermal) then e=shift(e, shift)
if keyword_set(mhd) then begin
    bx=shift(bx, shift)
    by=shift(by, shift)
    bz=shift(bz, shift)
endif
if keyword_set(ion) then d_n=shift(d_n, shift)
if keyword_set(nscalar) then s=shift(s, [0, shift[0], shift[1], shift[2]])

if n_elements(parallel) ne 3 then parallel=[1,1,1]

; now write data to output files
nproc=parallel[0]*parallel[1]*parallel[2]
xb=-1
yb=0
zb=0
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

                                ; figure out where this block goes
    xb++
    if xb eq parallel[0] then begin
        xb=0
        yb++
        if yb eq parallel[1] then begin
            yb=0
            zb++
        endif
    endif
    xlo=nxb*xb
    ylo=nyb*yb
    zlo=nzb*zb

                                ; open input and output files
    openr, fp, newname, /get_lun
    openw, fpout, newoutname, /get_lun
    print, 'Converting '+newname+' to '+newoutname+'...'

                                ; copy header
    while dummy ne '<par_end>' do begin
        readf, fp, dummy
        outstr=dummy
        printf, fpout, outstr
    endwhile
    readf, fp, dummy
    printf, fpout, dummy
    nstep=0L
    readu, fp, nstep
    writeu, fpout, nstep
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    time=0.0d0
    readu, fp, time
    writeu, fpout, time
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    dt=0.0d0
    readu, fp, dt
    writeu, fpout, dt

                                ; set grid size for parallel runs
    block = dblarr(nxb, nyb, nzb)

    ; write data
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    readu, fp, block
    writeu, fpout, d[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
        readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    readu, fp, block
    writeu, fpout, px[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    readu, fp, block
    writeu, fpout, py[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy
    readu, fp, block
    writeu, fpout, pz[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    if not keyword_set(isothermal) then begin
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, block
        writeu, fpout, e[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    endif
    if keyword_set(mhd) then begin
	xfblock = dblarr(nxb+1,nyb,nzb)
	yfblock = dblarr(nxb,nyb+1,nzb)
	zfblock = dblarr(nxb,nyb,nzb+1)
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, xfblock
        writeu, fpout, bx[xlo:xlo+nxb,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, yfblock
        writeu, fpout, by[xlo:xlo+nxb-1,ylo:ylo+nyb,zlo:zlo+nzb-1]
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, zfblock
        writeu, fpout, bz[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb]
    endif
    if keyword_set(ion) then begin
        dummyd=0.0d0
        dummyf=0.0
        dummyl=0L
        nrad=0L
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, nrad
        writeu, fpout, nrad
        xrad=dblarr(nrad, 3)
        srad=dblarr(nrad)
        for m=1,nrad do begin
            for l=1,3 do begin
                readu, fp, dummyd
                writeu, fpout, dummyd
                xrad[m-1,l-1]=dummyd
            endfor
            readu, fp, dummyd
            writeu, fpout, dummyd
            srad[m-1]=dummyd
            readu, fp, dummyl
            writeu, fpout, dummyl
            for l=1,9 do begin
                readu, fp, dummyf
                writeu, fpout, dummyf
            endfor
        endfor
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        for m=1,70 do begin
            readu, fp, dummyl
            writeu, fpout, dummyl
        endfor
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        readu, fp, block
        writeu, fpout, d_n[xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    endif
    if keyword_set(nscalar) then begin
        readf, fp, dummy
        printf, fpout, dummy
        readf, fp, dummy
        printf, fpout, dummy
        sblock=dblarr(nscalar,nxb,nyb,nzb)
        readu, fp, sblock
        writeu, fpout, s[*,xlo:xlo+nxb-1,ylo:ylo+nyb-1,zlo:zlo+nzb-1]
    endif
    readf, fp, dummy
    printf, fpout, dummy
    readf, fp, dummy
    printf, fpout, dummy

    free_lun, fp
    free_lun, fpout
    
endfor

end
