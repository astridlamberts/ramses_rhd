pro read_vtk_ath, fname, d, px, py, pz, e=e, isothermal=isothermal, $
	mhd=mhd, s=s, d_n=d_n, ion=ion, nscalar=nscalar, $
	bx=bx, by=by, bz=bz, time=time, dt=dt, parallel=parallel, $
	xvec=xvec, yvec=yvec, zvec=zvec, nofloor=nofloor

; open file
openr, fp, fname, /get_lun

; read time
dummy=''
readf, fp, dummy
readf, fp, dummy
pos=strpos(dummy,'=')+1
time=float(strmid(dummy, pos))

; make sure data is in binary format
readf, fp, dummy
if dummy ne 'BINARY' then message, 'Error: can only read binary vtk files'

; read size of data cube
readf, fp, dummy
readf, fp, dummy
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
nx=fix(dummy)-1
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
ny=fix(dummy)-1
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
nz=fix(dummy)-1

; read computational grid
readf, fp, dummy
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
xorig=fix(dummy)
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
yorig=fix(dummy)
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
zorig=fix(dummy)

; read cell spacing
readf, fp, dummy
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
dx=fix(dummy)
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
dy=fix(dummy)
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
dz=fix(dummy)

; read number of cells
readf, fp, dummy
pos=strpos(dummy, ' ')+1
dummy=strmid(dummy, pos)
ncell=long(dummy)

; check if this is a parallel run
if n_elements(parallel) eq 3 then begin
	block = fltarr(nx,ny,nz)
	xblock = fltarr(nx)
	yblock = fltarr(ny)
	zblock = fltarr(nz)
	nxgrid = nx
	nygrid = ny
	nzgrid = nz
	nx = nx*parallel[0]
	ny = ny*parallel[1]
	nz = nz*parallel[2]
	nproc = parallel[0]*parallel[1]*parallel[2]
endif

; allocate memory and read data
d =fltarr(nx,ny,nz)
px=fltarr(nx,ny,nz)
py=fltarr(nx,ny,nz)
pz=fltarr(nx,ny,nz)
vel=fltarr(3,nx,ny,nz)
xvec=fltarr(nx)
yvec=fltarr(ny)
zvec=fltarr(nz)
if n_elements(parallel) ne 3 then begin
    xvec=xorig+dx*(findgen(nx)+0.5)
    yvec=yorig+dy*(findgen(ny)+0.5)
    zvec=zorig+dz*(findgen(nz)+0.5)
    readf, fp, dummy
    readf, fp, dummy
    readu, fp, d
    swap_endian_inplace, d, /swap_if_little_endian
    readf, fp, dummy
    readu, fp, vel
    swap_endian_inplace, vel, /swap_if_little_endian
    px=reform(vel[0,*,*,*])*d
    py=reform(vel[1,*,*,*])*d
    pz=reform(vel[2,*,*,*])*d
    if not keyword_set(isothermal) then begin
        readf, fp, dummy
        readf, fp, dummy
	e=fltarr(nx,ny,nz)
	readu, fp, e
        swap_endian_inplace, e, /swap_if_little_endian
    endif
    if keyword_set(mhd) then begin
        bx=fltarr(nx,ny,nz)
        by=fltarr(nx,ny,nz)
        bz=fltarr(nx,ny,nz)
        readf, fp, dummy
        readu, fp, bx
        readu, fp, by
        readu, fp, bz
        swap_endian_inplace, bx, /swap_if_little_endian
        swap_endian_inplace, by, /swap_if_little_endian
        swap_endian_inplace, bz, /swap_if_little_endian
    endif
    if keyword_set(nscalar) then begin
        s=fltarr(nscalar,nx,ny,nz)
        readf, fp, dummy
        readf, fp, dummy
        readu, fp, s
        swap_endian_inplace, s, /swap_if_little_endian
    endif
endif else begin
	readu, fp, xblock
	xvec[0:nxgrid-1]=xblock
	readu, fp, yblock
	yvec[0:nygrid-1]=yblock
	readu, fp, zblock
	zvec[0:nzgrid-1]=zblock
	readu, fp, block
	d[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	readu, fp, block
	px[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	readu, fp, block
	py[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	readu, fp, block
	pz[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	if not keyword_set(isothermal) then begin
		e=fltarr(nx,ny,nz)
		readu, fp, block
		e[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	endif
	if keyword_set(mhd) then begin
		bx=fltarr(nx,ny,nz)
		by=fltarr(nx,ny,nz)
		bz=fltarr(nx,ny,nz)
		readu, fp, block
		bx[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
		readu, fp, block
		by[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
		readu, fp, block
		bz[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	endif
	if keyword_set(ion) then begin
		d_n=fltarr(nx,ny,nz)
		readu, fp, block
		d_n[0:nxgrid-1,0:nygrid-1,0:nzgrid-1]=block
	endif
	if keyword_set(nscalar) then begin
		s=fltarr(nscalar,nx,ny,nz)
		sblock=fltarr(nscalar,nxgrid,nygrid,nzgrid)
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
		readu, fp, ndata
		readu, fp, eos
		readu, fp, tdt

		; read data
		readu, fp, xblock
		xvec[xlo:xlo+nxgrid-1]=xblock
		readu, fp, yblock
		yvec[ylo:ylo+nygrid-1]=yblock
		readu, fp, zblock
		zvec[zlo:zlo+nzgrid-1]=zblock
		readu, fp, block
		d[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		readu, fp, block
		px[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		readu, fp, block
		py[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		readu, fp, block
		pz[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		if not keyword_set(isothermal) then begin
			readu, fp, block
			e[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		endif
		if keyword_set(mhd) then begin
			readu, fp, block
			bx[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
			readu, fp, block
			by[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
			readu, fp, block
			bz[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		endif
		if keyword_set(ion) then begin
			readu, fp, block
			d_n[xlo:xlo+nxgrid-1,ylo:ylo+nygrid-1,zlo:zlo+nzgrid-1]=block
		endif
		if keyword_set(nscalar) then begin
			sblock=fltarr(nscalar,nxgrid,nygrid,nzgrid)
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

end
