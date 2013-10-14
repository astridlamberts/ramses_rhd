;how to plot a raster and add a decent colorbar

;b_x=bx[*,*,128]
;b_y=by[*,*,128]

loadct,13 ;gives the rainbow colors
psopen,'run8.2dneutral.ps',/color

;t=athtemp(d,d*vx,d*vy,d*vz,e,bx=bx,by=by,bz=bz,d_n=dn)
;print, mean(t)
 raster,dn[*,*,0],/log,[-5,5],[-5,5],xsty=1,ysty=1 ,zra=[-23,-21], title='neutral density in the XY plane in the middle of the box', xtitle='X axis (pc)',ytitle='Y axis (pc)', charsize=2.2
 
;vel,b_x,b_y, length=0.2,title='B in the XY plane in the middle of the box';,xtitle='B along X',xtitle='B along Y';, charsize=2.8

;dont add too much divisions because there is not eneough precision on
;the scale
colorbar,position=[0.18,0.88,0.95,0.95], maxrange=-20, minrange=-23, divisions=2

psclose








end
