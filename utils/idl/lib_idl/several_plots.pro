;size_domain=256
dt=2.0e12
pos=findgen(256)*0.078
read_bin_ath, 'ifront_1e9_100_long_8.0030.bin',d,px,py,pz,e=e,d_n=d_n,/ion
t=athtemp(d,px,py,pz,e,d_n=d_n)

psopen,'variables.ps'
!p.multi=[0,2,2]
!X.MARGIN=[8,1]
!Y.MARGIN=[3,3]
plot,pos,d[*,0,0],/ylog,title='density',xtitle='distance(pc)',ytitle='d (g/cm^3)',charsize=2.5
plot,pos,d_n[*,0,0],/ylog,title='neutral density',xtitle='distance(pc)',ytitle='d (g /cm^3)',charsize=2.5
plot,pos,t[*,0,0],/ylog,title='temperature',xtitle='distance(pc)',ytitle='T (K)',charsize=2.5


plot,pos,e[*,0,0],/ylog,title='energy',xtitle='distance(pc)',ytitle='energy (erg=10^-7J)',charsize=2.5











psclose

end
