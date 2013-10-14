pro fieldline2d, bx, by, xlim, ylim, seedpts=seedpts, nseed=nseed, $
	stepsize=stepsize, maxpts=maxpts, arrowpt=arrowpt, $
	periodic=periodic, linelen=linelen, color=color

sz=size(bx)
nx=sz[1]
ny=sz[2]
if n_elements(xlim) eq 0 then xlim=[0,nx]
if n_elements(ylim) eq 0 then ylim=[0,ny]
dx=float(xlim[1]-xlim[0])/nx
dy=float(ylim[1]-ylim[0])/ny
delta=[dx,dy]
bxmin=[xlim[0],ylim[0]]
if n_elements(stepsize) eq 0 then stepsize=0.01
step=min(delta)*stepsize
if not keyword_set(maxpts) then maxpts=max([nx,ny])/stepsize
if n_elements(arrowpt) eq 0 then arrowpt=0.5
if n_elements(color) eq 0 then color=!d.table_size-11
if n_elements(periodic) eq 0 then periodic=[0,0]

; if necessary, generate seed points, placing them uniformly along the x and y
; borders of the computational domain
if not keyword_set(seedpts) then begin
	if not keyword_set(nseed) then nseed=20
	seedpts=fltarr(nseed,2)
	seedpts[0:nseed/2-1,0] = (findgen(nseed/2)+0.5)*(xlim[1]-xlim[0]) / $
		(nseed/2) + xlim[0]
	seedpts[0:nseed/2-1,1] = ylim[0]+dy/100.
	seedpts[nseed/2:nseed-1,0] = xlim[0]+dx/100.
	seedpts[nseed/2:nseed-1,1] = (findgen(nseed-nseed/2)+0.5) * $
		(ylim[1]-ylim[0]) / $
		(nseed-nseed/2) + ylim[0]
endif else nseed=n_elements(seedpts)/2

for n=0,nseed-1 do begin

	; preparation work
	line=fltarr(maxpts,2)
	lineneg=fltarr(maxpts,2)

	; walk field line in both directions until reaching edge of domain
	; or the length limits
	pos=reform(seedpts[n,*])
	line[0,*]=pos
	npts=1
	while ((pos[0] gt xlim[0]) and (pos[0] lt xlim[1]) and $
	       (pos[1] gt ylim[0]) and (pos[1] lt ylim[1]) and $
	       (npts lt maxpts)) do begin

		; flag for periodic looping
		loop = 0

		; find cell center below this point
		idx0 = fix(floor((pos - bxmin)/delta - 0.5))
		idx1 = idx0+[1,0]
		idx2 = idx0+[0,1]
		idx3 = idx0+[1,1]

		; get weight for interpolation
		wgt = (pos - (bxmin + delta*(idx0+0.5)))/delta

		; handle special case when point is near grid edge
		lopt=periodic*([nx,ny]-1)
		hipt=(1-periodic)*([nx,ny]-1)
		idx0 = idx0*(idx0 ge 0) + lopt*(idx0 lt 0)
		idx1 = idx1*(idx1 ge 0) + lopt*(idx1 lt 0)
		idx1 = idx1*(idx1 lt [nx,ny]) + hipt*(idx1 ge [nx,ny])
		idx2 = idx2*(idx2 ge 0) + lopt*(idx2 lt 0)
		idx2 = idx2*(idx2 lt [nx,ny]) + hipt*(idx2 ge [nx,ny])
		idx3 = idx3*(idx3 lt [nx,ny]) + hipt*(idx3 ge [nx,ny])

		; get current direction
		xvec = (1-wgt[0])*(1-wgt[1])*bx[idx0[0],idx0[1]] + $
		       wgt[0]*(1-wgt[1])*bx[idx1[0],idx1[1]] + $
		       (1-wgt[0])*wgt[1]*bx[idx2[0],idx2[1]] + $
		       wgt[0]*wgt[1]*bx[idx3[0],idx3[1]]
		yvec = (1-wgt[0])*(1-wgt[1])*by[idx0[0],idx0[1]] + $
		       wgt[0]*(1-wgt[1])*by[idx1[0],idx1[1]] + $
		       (1-wgt[0])*wgt[1]*by[idx2[0],idx2[1]] + $
		       wgt[0]*wgt[1]*by[idx3[0],idx3[1]]
		mag=sqrt(xvec^2+yvec^2)

		; increment position
		pos = pos + [xvec,yvec]/mag*stepsize

		; loop if periodic
		if total(periodic) ne 0 then begin
		    if (periodic[0] eq 1) then begin
			if (pos[0] le xlim[0]) then begin
				pos[0]=pos[0]+(xlim[1]-xlim[0])
				loop = 1
			endif
			if (pos[0] ge xlim[1]) then begin
				pos[0]=pos[0]-(xlim[1]-xlim[0])
				loop = 1
			endif
		    endif
		    if (periodic[1] eq 1) then begin
			if (pos[1] le ylim[0]) then begin
				pos[1]=pos[1]+(ylim[1]-ylim[0])
				loop = 1
			endif
			if (pos[1] ge ylim[1]) then begin
				pos[1]=pos[1]-(ylim[1]-ylim[0])
				loop = 1
			endif
		    endif
		endif

		; record position
		line[npts,*] = pos
		npts = npts+1

		if ((npts gt 1) and (loop eq 0)) then $
		    plots, [line[npts-2:npts-1,0]], $
			   [line[npts-2:npts-1,1]], $
		           noclip=0, color=color

	endwhile

	pos=reform(seedpts[n,*])
	nptsneg=0
	while ((pos[0] gt xlim[0]) and (pos[0] lt xlim[1]) and $
	       (pos[1] gt ylim[0]) and (pos[1] lt ylim[1]) and $
	       (nptsneg lt maxpts)) do begin

		; flag for periodic looping
		loop = 0

		; find cell center below this point
		idx0 = fix((pos - bxmin)/delta - 0.5)
		idx1 = idx0+[1,0]
		idx2 = idx0+[0,1]
		idx3 = idx0+[1,1]

		; get weight for interpolation
		wgt = (pos - (bxmin + delta*(idx0+0.5)))/delta

		; handle special case when point is near grid edge
		lopt=periodic*([nx,ny]-1)
		hipt=(1-periodic)*([nx,ny]-1)
		idx0 = idx0*(idx0 ge 0) + lopt*(idx0 lt 0)
		idx1 = idx1*(idx1 ge 0) + lopt*(idx1 lt 0)
		idx1 = idx1*(idx1 lt [nx,ny]) + hipt*(idx1 ge [nx,ny])
		idx2 = idx2*(idx2 ge 0) + lopt*(idx2 lt 0)
		idx2 = idx2*(idx2 lt [nx,ny]) + hipt*(idx2 ge [nx,ny])
		idx3 = idx3*(idx3 lt [nx,ny]) + hipt*(idx3 ge [nx,ny])

		; get current direction
		xvec = (1-wgt[0])*(1-wgt[1])*bx[idx0[0],idx0[1]] + $
		       wgt[0]*(1-wgt[1])*bx[idx1[0],idx1[1]] + $
		       (1-wgt[0])*wgt[1]*bx[idx2[0],idx2[1]] + $
		       wgt[0]*wgt[1]*bx[idx3[0],idx3[1]]
		yvec = (1-wgt[0])*(1-wgt[1])*by[idx0[0],idx0[1]] + $
		       wgt[0]*(1-wgt[1])*by[idx1[0],idx1[1]] + $
		       (1-wgt[0])*wgt[1]*by[idx2[0],idx2[1]] + $
		       wgt[0]*wgt[1]*by[idx3[0],idx3[1]]
		mag=sqrt(xvec^2+yvec^2)

		; increment position
		pos = pos - [xvec,yvec]/mag*stepsize

		; loop if periodic
		if total(periodic) ne 0 then begin
		    if (periodic[0] eq 1) then begin
			if (pos[0] le xlim[0]) then begin
				pos[0]=pos[0]+(xlim[1]-xlim[0])
				loop = 1
			endif
			if (pos[0] ge xlim[1]) then begin
				pos[0]=pos[0]-(xlim[1]-xlim[0])
				loop = 1
			endif
		    endif
		    if (periodic[1] eq 1) then begin
			if (pos[1] le ylim[0]) then begin
				pos[1]=pos[1]+(ylim[1]-ylim[0])
				loop = 1
			endif
			if (pos[1] ge ylim[1]) then begin
				pos[1]=pos[1]-(ylim[1]-ylim[0])
				loop = 1
			endif
		    endif
		endif

		; record position
		lineneg[nptsneg,*] = pos
		nptsneg = nptsneg+1

		if ((nptsneg gt 1) and (loop eq 0)) then $
		    plots, [lineneg[nptsneg-2:nptsneg-1,0]], $
			   [lineneg[nptsneg-2:nptsneg-1,1]], $
		           noclip=0, color=color

	endwhile

	; join lines
	xline = [reverse(lineneg[0:nptsneg-1,0]), line[0:npts-1,0]]
	yline = [reverse(lineneg[0:nptsneg-1,1]), line[0:npts-1,1]]

	; draw the lines
	;plots, xline, yline, noclip=0, color=color

	; draw arrows
	npts=npts+nptsneg
	if arrowpt[0] ge 0.0 then begin
		for i=0,n_elements(arrowpt)-1 do $
		    arrow, xline[npts*arrowpt[i]], yline[npts*arrowpt[i]], $
			xline[npts*arrowpt[i]+1], yline[npts*arrowpt[i]+1], $
			/data, color=color
	endif

endfor

end

