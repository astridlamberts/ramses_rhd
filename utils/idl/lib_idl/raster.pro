pro raster, arr, xlim, ylim, zrange=zrange, $
	log=log, anisotropic=anisotropic, $
	xrange=xrange, yrange=yrange, $
	background=background, imgout=imgout, $
	noplot=noplot, _extra=extra

; read keywords
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0

; check we have a 2d array, or one that reduces to one, then get limits
sz=size(arr)
if sz[0] ne 2 then begin
    reformflag=1
    szorig=sz
    arr=reform(arr)
    sz=size(arr)
    if sz[0] ne 2 then begin
        arr=reform(arr, szorig[1:szorig[0]])
        print, 'Error: input array must be two-dimensional!'
        return
    endif
endif else reformflag=0
if n_elements(xlim) eq 0 then xlim1=[0,sz[1]-1] $
else if n_elements(xlim) eq 1 then xlim1=[0.,xlim] $
else xlim1=xlim

if n_elements(ylim) eq 0 then ylim1=[0,sz[2]-1] $
else if n_elements(ylim) eq 1 then ylim1=[0.,ylim] $
else ylim1=ylim

; get range
if not keyword_set(xrange) then xrange=xlim1
if not keyword_set(yrange) then yrange=ylim1

; get some information about the graphics environment
bangmulti = !p.multi
tempmulti = !p.multi
tempmulti[1:3] = tempmulti[1:3] + (tempmulti[1:3] eq 0)

; Convert data to normal coordinates. If we're on a device that uses
; windows, make sure a valid window is open before attempting this.
devname = !d.name
devname = strmid(devname, 0, 3)
devname = strupcase(devname)
case devname of
	'MAC': if !d.window eq -1 then window, 0
	'WIN': if !d.window eq -1 then window, 0
	'X': if !d.window eq -1 then window, 0
	ELSE:
endcase

; see how many colors are available
ncolors = !d.table_size

; draw axes
plot, [0,0], [0,0], /nodata, xrange=xrange, yrange=yrange, $
	isotropic=isotropic, _extra = extra

; If necessary, clip array to retain only part that is in 
; requested range
if total(xrange ne xlim1) + total(yrange ne ylim1) ne 0 then begin

	; figure out pixel sizes
	sz=size(arr)
	inputpx=float([xlim1[1]-xlim1[0], ylim1[1]-ylim1[0]])/sz[1:2]

	; clip image if necessary
	if xlim1[0] lt xrange[0] then begin
		xlimcl = fix((xrange[0]-xlim1[0])/inputpx[0])
		arrclip = arr[xlimcl:sz[1]-1,*]
		xlim1[0]=xrange[0]
	endif else arrclip=arr
	if ylim1[0] lt yrange[0] then begin
		ylimcl = fix((yrange[0]-ylim1[0])/inputpx[1])
		arrclip = arrclip[*,ylimcl:sz[2]-1]
		ylim1[0]=yrange[0]
	endif
	if xlim1[1] gt xrange[1] then begin
		xlimcl = (xrange[1]-xlim1[0])/inputpx[0]
		arrclip = arrclip[0:xlimcl,*]
		xlim1[1]=xrange[1]
	endif
	if ylim1[1] gt yrange[1] then begin
		ylimcl = (yrange[1]-ylim1[0])/inputpx[1]
		arrclip = arrclip[*,0:ylimcl]
		ylim1[1]=yrange[1]
	endif

	; get pixel sizes of clipped array
	sz=size(arrclip)
	inputpx=[xlim1[1]-xlim1[0], ylim1[1]-ylim1[0]]/sz[1:2]

	; record if limits match after clipping
	limmatch = (total(xrange ne xlim1) + total(yrange ne ylim1)) eq 0

endif else begin
	arrclip=arr
	limmatch=1
endelse

; What we do here depends on whether we have a device with scalable
; pixels or not. If we have scalable pixels, we can just construct
; an image box where each cell of the grid
; represents a single pixel. We then rely on the scaling of pixels
; to make that image correctly overlay the plot window. If we don't
; have scalable pixels, the number of cells per pixel will in be a
; complicated function of the size of the graphics device and the
; resolution of the amr grid. For reference, X windows does not have
; scalable pixels, postscript does.

scalablePixels = !d.flags mod 2

if scalablePixels then begin

	img=arrclip

endif else begin

	; we do not have scalable pixels

	; get max and min for image box in pixels
	imglo = convert_coord(xrange[0], yrange[0], /data, /to_device)
	imghi = convert_coord(xrange[1], yrange[1], /data, /to_device)
	imglo=round(imglo[0:1])
	imghi=round(imghi[0:1])
	if (imghi[0] - imglo[0]) lt 2 then begin
		print, 'Error: requested x axis range is too small.'
		return
	endif
	if (imghi[1] - imglo[1]) lt 2 then begin
		print, 'Error: requested y axis range is too small.'
		return
	endif

	; Do image box and array box match? If so, this case is simple.
	if limmatch eq 1 then begin

		; map input array to image box
		img=congrid(arrclip, imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)

	endif else begin

		; ranges don't match, so set up an image box and then
		; put the input image into part of it

		; prepare image array
		img=dblarr(imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)
		if keyword_set(background) then img = img + background
		npxout=[imghi[0]-imglo[0], imghi[1]-imglo[1]] + 1
		outputpx=float([xrange[1]-xrange[0], yrange[1]-yrange[0]]) / $
	    		npxout

		; rescale input image
		arrscale = congrid(arrclip, inputpx[0]/outputpx[0]*sz[1], $
			inputpx[1]/outputpx[1]*sz[2])
		szscale=size(arrscale)

		; figure out offset to place this image into output
		; box
		offset = [xlim1[0]-xrange[0], ylim1[0]-yrange[0]] / $
			outputpx[0]

		; handle off by one errors in rounding pixel numbers
		if offset[0]+szscale[1] gt npxout[0] then $
			offset[0] = offset[0] - 1
		if offset[1]+szscale[2] gt npxout[1] then $
			offset[1] = offset[1] - 1

		; put result into image box
		img[ offset[0] : offset[0]+szscale[1]-1, $
		     offset[1] : offset[1]+szscale[2]-1 ] = arrscale

	endelse
endelse

; if requested, take the log of the image
if keyword_set(log) then begin
	img=alog10(img)
	; do something sensible if the log produces NaN's
	idx1=where(1-finite(img))
	idx2=where(finite(img))
	if (idx1[0] ne -1) and (idx2[0] ne -1) then img[idx1]=min(img[idx2])-0.5*mean(img[idx2])
endif

; set the range to be used to assign colors
if n_elements(zrange) eq 0 then begin
	maxrange=max(img)
	minrange=min(img)
endif else begin
	maxrange=zrange[1]
	minrange=zrange[0]
endelse

; put pixels into specified range
if n_elements(zrange) ne 0 then begin
	img = (img lt minrange) * minrange + (img ge minrange) * img
	img = (img gt maxrange) * maxrange + (img le maxrange) * img
endif

; scale the image to the color table
imgout = fix((img - minrange) / (maxrange - minrange) * (ncolors-1))

; Turn off device decomposition for scalable pixels
if ((!d.name eq 'X') or (!d.name eq 'WIN')) then begin
	if !Version.Release ge 5.2 then device, get_decomposed=thisdecomposed
	device, decomposed=0
endif

; Now display the image. How we do this again depends on whether we
; have scalable pixels or not.

if scalablePixels then begin

	; we have scalable pixels

	; get box limits in normal coordinates
	imglonorm = convert_coord(xlim1[0], ylim1[0], /data, /to_norm)
	imghinorm = convert_coord(xlim1[1], ylim1[1], /data, /to_norm)

	; show image
	tv, imgout, imglonorm[0], imglonorm[1], $
		xsize = imghinorm[0] - imglonorm[0], $
		ysize = imghinorm[1] - imglonorm[1], $
		/norm

endif else begin

	; we do not have scalable pixels
	tv, imgout, imglo[0], imglo[1]
endelse

; turn decomposition back to previous setting
if ((!d.name eq 'X') or (!d.name eq 'WIN')) then begin
	if !Version.Release ge 5.2 then device, decomposed=thisdecomposed
endif

if keyword_set(verbose) then begin
	print, 'imglo0 = ', imglo[0] / !D.X_SIZE
	print, 'imglo1 = ', imglo[1] / !D.Y_SIZE
	print, 'imghi0 = ', imghi[0] / !D.X_SIZE
	print, 'imglo1 = ', imghi[1] / !D.Y_SIZE
endif

; overlay the plot axes again
!p.multi=bangmulti
plot, [0,0], [0,0], /noerase, /nodata, xrange=xrange, yrange=yrange, $
	isotropic=isotropic, _extra = extra

; increment !p.multi, since the plot call with /noerase set doesn't do
; that automatically
if !p.multi[0] eq 0 then !p.multi[0]=tempmulti[1]*tempmulti[2]*tempmulti[3]
!p.multi[0] = !p.multi[0] - 1

; return original array if we reformed it
if reformflag eq 1 then arr=reform(arr, szorig[1:szorig[0]])

return

end


