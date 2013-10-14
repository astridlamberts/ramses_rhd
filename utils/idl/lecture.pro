*; routine qui sert a la lecture des sorties, elle lit la grille
; (rd_amr), les variables hydro (rd_hydro), met les cellules
; dans l'ordre (amr2cell) et affiche ce qu'on veut qu'elle
; affiche (tv2d avec type)

pro lecture, grid, hydro,type=type, log=log, xr=xr,yr=yr, save=data,clt=clt,nout=nout
;data=0.0
data=fltarr(256,256)
rd_amr,grid,nout=nout
rd_hydro,hydro,nout=nout
amr2cell, grid, hydro, cell

tv2d, grid, hydro, type=type, log=log, xr=xr,yr=yr, save=data,clt=clt

azrrr


end
