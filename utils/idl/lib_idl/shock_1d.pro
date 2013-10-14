pro shock

rd_1d,data
out1=*data[1]
pos=findgen(out1.nc)/out1.nc

psopen,'shock_1D.ps'
!p.multi=[0,2,2]
!X.MARGIN=[8,1]
!Y.MARGIN=[3,3]


plot,pos,out1.d,psym=2,xtitle='position x',ytitle='densite',title='rho ',charsize=2.5

plot,pos,out1.u,psym=2,xtitle='position x',ytitle='vitesse',title='u',charsize=2.5

plot,pos,out1.P,psym=2,xtitle='position x',ytitle='pression',title='P',charsize=2.5

plot,pos,out1.l,xtitle='position x',ytitle='L',title='niveau de raffinement',charsize=2.5

psclose

end
