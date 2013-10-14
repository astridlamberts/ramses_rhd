
pos=findgen(256)*0.078
psopen,'density_gap.ps'
read_bin_ath, 'ifront_256_3e10_100.0035.bin',d,px,py,pz,e=e,d_n=d_n,/ion
plot,pos,d[*,0,0],/ylog,title='density for different timesteps',xtitle='distance (pc)',ytitle='density (g/cm^3)',charsize=3


read_bin_ath, 'ifront_256_3e10_100.0020.bin',d,px,py,pz,e=e,d_n=d_n,/ion
oplot,pos,d[*,0,0],linestyle=2
read_bin_ath, 'ifront_256_3e10_100.0005.bin',d,px,py,pz,e=e,d_n=d_n,/ion
oplot,pos,d[*,0,0],linestyle=1

psclose
end
