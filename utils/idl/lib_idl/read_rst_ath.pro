pro read_rst_ath, fname, d, px, py, pz, e=e, isothermal=isothermal, $
	mhd=mhd, s=s, d_n=d_n, ion=ion, nscalar=nscalar, $
	bx=bx, by=by, bz=bz, time=time, dt=dt, parallel=parallel, $
	xvec=xvec, yvec=yvec, zvec=zvec, nofloor=nofloor, nofloat=nofloat, $
                  nrad=nrad, xrad=xrad, srad=srad

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
	block = dblarr(nxgrid,nygrid,nzgrid)
	nproc = parallel[0]*parallel[1]*parallel[2]
endif

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

; as a convenience, floor negative d_n's unless requested otherwise
if keyword_set(ion) and not keyword_set(nofloor) then begin
    badidx=where(d_n le 0.0)
    if (badidx[0] ne -1) then begin
        goodidx=where(d_n gt 0.0)
        d_n[badidx]=d[badidx]*min(d_n[goodidx]/d[goodidx])
    endif
endif

; average b field to cell centers
if keyword_set(mhd) then begin
    bx=0.5*(bx[0:nx-1,*,*]+bx[1:nx,*,*])
    by=0.5*(by[*,0:ny-1,*]+by[*,1:ny,*])
    bz=0.5*(bz[*,*,0:nz-1]+bz[*,*,1:nz])
endif

; convert to float to save memory
if not keyword_set(nofloat) then begin
    d=float(d)
    px=float(px)
    py=float(py)
    pz=float(pz)
    if not keyword_set(isothermal) then e=float(e)
    if keyword_set(ion) then d_n=float(d_n)
    if keyword_set(mhd) then begin
        bx=float(bx)
        by=float(by)
        bz=float(bz)
    endif
endif

end
