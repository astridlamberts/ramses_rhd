;raster,d[*,*,100],[-5,5],[-5,5],xsty=1,ysty=1,zra=[0,5]


set_plot,'ps'
device,filename='dn_700000yr.ps',/color,bits_per_pixel=8
loadct,13
!p.thick=2
!x.thick=2
!y.thick=2
Xsize=15
raster,dn[*,*,100],[-5,5],[-5,5],xsty=1,ysty=1,zra=[-24,-21],/log,title='Neutral density after 700 000 yr in slice 100', xtitle='distance (pc)',ytitle='distance(pc)'
margin=0.12
colorbar,maxrange=-21,minrange=-24,/vertical

psclose

end
