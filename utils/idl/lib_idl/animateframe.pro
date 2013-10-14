nframes=500

; get list of image files
spawn, 'ls *.frm', frmlst
nfrm=n_elements(frmlst)

; get list of image times
times=dblarr(nfrm)
tbreak=dblarr(nfrm+1)
time=0.0d0
for n=0, nfrm-1 do begin
    pltname=strmid(frmlst[n], 0, strlen(frmlst[n])-4)
    tname=pltname+'.time'
    openr, fp, tname, /get_lun
    readf, fp, time
    times[n]=time
    free_lun, fp
endfor

; get time breaks between images
tbreak[0]=times[0]
for n=1,nfrm-1 do tbreak[n]=0.5*(times[n-1]+times[n])
tbreak[nfrm]=times[nfrm-1]

; get time step
dt=tbreak[nfrm]/(nframes-1)
frameptr=0

; read first image into memory
print, 'Reading image file '+frmlst[frameptr]+', file 1 of '+strtrim(string(nframes),2)+'...'
openr, fp, frmlst[frameptr], /get_lun
wsize=intarr(2)
readu, fp, wsize
if wsize[0] lt 0 then eswap=1 else eswap=0
if eswap then wsize=swap_endian(wsize)
frm=bytarr(wsize[0], wsize[1])
readu, fp, frm
if eswap then frm=swap_endian(frm)
free_lun, fp

; flip frame: for some reason this seems to be necessary to avoid
; upside-down movies
frm=reverse(frm, 2)

for n=0, nframes-1 do begin

    ; get time and read next image file if necessary
    time=n*dt
    if (time gt tbreak[frameptr+1]) then begin
        frameptr=max(where(time gt tbreak))
        print, 'Reading image file '+frmlst[frameptr]+', file '+strtrim(string(n+1),2)+' of '+strtrim(string(nframes),2)+'...'
        openr, fp, frmlst[frameptr], /get_lun
        wsize=intarr(2)
        readu, fp, wsize
        if wsize[0] lt 0 then eswap=1 else eswap=0
        if eswap then wsize=swap_endian(wsize)
        frm=bytarr(wsize[0], wsize[1])
        readu, fp, frm
        if eswap then frm=swap_endian(frm)
        free_lun, fp

	; flip frame: for some reason this seems to be necessary to avoid
	; upside-down movies
	frm=reverse(frm, 2)

    endif
       
    ; initialize mpeg if necessary
    if n eq 0 then begin
	mpegid=mpeg_open(wsize, filename='movie.mpg', quality=80)
	loadct, 39
	tvlct, r, g, b, /get
    endif

    ; convert to 24 bit
    img24=bytarr(3,wsize[0],wsize[1])
    img24[0,*,*]=r(frm)
    img24[1,*,*]=g(frm)
    img24[2,*,*]=b(frm)

    ; add frame
    mpeg_put, mpegid, image=img24, frame=n

endfor

; save mpeg movie
mpeg_save, mpegid
mpeg_close, mpegid

end
