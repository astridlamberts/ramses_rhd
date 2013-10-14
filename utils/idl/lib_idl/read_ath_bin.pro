pro read_ath_bin, fname, nx, ny, nz, d, px, py, pz, e, isothermal=isothermal, mhd=mhd, s=s, d_n=d_n, ray_radiation=ray_radiation, nscalar=nscalar, parallel=parallel, bx=bx, by=by, bz=bz

; open file
openr, fp, fname, /get_lun

; read size of data cube
ndata=LONARR(4)
readu, fp, ndata
nvar=ndata[3]
nx=ndata[0]
ny=ndata[1]
nz=ndata[2]

; if parallel, increase data cube size appropriately
if n_elements(parallel) eq 3 then begin
	block = fltarr(nx,ny,nz)
	nxblock = nx
	nyblock = ny
	nzblock = nz
	nx = nx*parallel[0]
	ny = ny*parallel[1]
	nz = nz*parallel[2]
	nproc = parallel[0]*parallel[1]*parallel[2]
endif

; check we have the expected number of variables
nvar = 4
if not keyword_set(isothermal) then nvar = nvar + 1
if keyword_set(nscalar) then nvar = nvar + nscalar
if keyword_set(ray_radiation) then nvar = nvar + 1
if keyword_set(mhd) then nvar = nvar + 3
if (nvar ne ndata[3]) then begin
	print, 'Error: expected '+strtrim(string(nvar),2)+' components, '+ $
		'found '+strtrim(string(ndata[3]),2)+'!'
	return
endif

; read eos
eos=fltarr(2)
readu, fp, eos
gamm1=eos[0]
isocs=eos[1]

; allocate memory and read
d =fltarr(nx,ny,nz)
px=fltarr(nx,ny,nz)
py=fltarr(nx,ny,nz)
pz=fltarr(nx,ny,nz)
if n_elements(parallel) ne 3 then begin
	readu, fp, d
	readu, fp, px
	readu, fp, py
	readu, fp, pz
	if not keyword_set(isothermal) then begin
		e=fltarr(nx,ny,nz)
		readu, fp, e
	endif
	if keyword_set(mhd) then begin
	    bx=fltarr(nx,ny,nz)
	    by=fltarr(nx,ny,nz)
	    bz=fltarr(nx,ny,nz)
	    readu, fp, bx
	    readu, fp, by
	    readu, fp, bz
	endif
	if keyword_set(ray_radiation) then begin
		d_n=fltarr(nx,ny,nz)
		readu, fp, d_n
	endif
	if keyword_set(nscalar) then begin
		s=fltarr(nscalar,nx,ny,nz)
		readu, fp, s
	endif
endif else begin
	readu, fp, block
	d[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	readu, fp, block
	px[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	readu, fp, block
	py[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	readu, fp, block
	pz[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	if not keyword_set(isothermal) then begin
		e=fltarr(nx,ny,nz)
		readu, fp, block
		e[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	endif
	if keyword_set(mhd) then begin
		bx=fltarr(nx,ny,nz)
		by=fltarr(nx,ny,nz)
		bz=fltarr(nx,ny,nz)
		readu, fp, block
		bx[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
		readu, fp, block
		by[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
		readu, fp, block
		bz[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	endif
	if keyword_set(ray_radiation) then begin
		d_n=fltarr(nx,ny,nz)
		readu, fp, block
		d_n[0:nxblock-1,0:nyblock-1,0:nzblock-1]=block
	endif
	if keyword_set(nscalar) then begin
		s=fltarr(nscalar,nx,ny,nz)
		sblock=fltarr(nscalar,nxblock,nyblock,nzblock)
		readu, fp, sblock
		s[*,0:nxblock-1,0:nyblock-1,0:nzblock-1]=sblock
	endif
endelse

; close file
free_lun, fp

; if parallel, read other files
if n_elements(parallel) eq 3 then begin
	xblock=0
	yblock=0
	zblock=0
	for n=1,nproc-1 do begin

		; construct file name
		fnamesplit=strsplit(fname, '.', /extract)
		nsplit=n_elements(fnamesplit)
		fnamesplit[nsplit-3]=fnamesplit[nsplit-3]+'-id'+strtrim(string(n),2)
		newname=''
		for i=0,nsplit-2 do newname+=fnamesplit[i]+'.'
		newname+=fnamesplit[nsplit-1]

		; figure out where this block goes
		xblock++
		if xblock eq parallel[0] then begin
			xblock=0
			yblock++
			if yblock eq parallel[1] then begin
				yblock=0
				zblock++
			endif
		endif
		xlo=nxblock*xblock
		ylo=nyblock*yblock
		zlo=nzblock*zblock

		; open and read header
		openr, fp, newname, /get_lun
		readu, fp, ndata
		readu, fp, eos

		; read data
		readu, fp, block
		d[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		readu, fp, block
		px[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		readu, fp, block
		py[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		readu, fp, block
		pz[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		if not keyword_set(isothermal) then begin
			readu, fp, block
			e[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		endif
		if keyword_set(mhd) then begin
			readu, fp, block
			bx[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
			readu, fp, block
			by[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
			readu, fp, block
			bz[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		endif
		if keyword_set(ray_radiation) then begin
			readu, fp, block
			d_n[xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=block
		endif
		if keyword_set(nscalar) then begin
			sblock=fltarr(nscalar,nxblock,nyblock,nzblock)
			readu, fp, sblock
			s[*xlo:xlo+nxblock-1,ylo:ylo+nyblock-1,zlo:zlo+nzblock-1]=sblock
		endif

		; close file
		free_lun, fp

	endfor

endif

end
