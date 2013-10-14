nbre_iterations=30      ;number of outputs
nbre_zones=256           ; number of cells
dt=5.0e9                ;timestep
Flux=1.0e9              ;flux
c_sound=7.80e5           ;sound speed in ionized region
n_H=double(100)          ;initial number density
alpha_b=3.92e-13         ;recombination coefficient
size_domain=3.09e17     ;size of the physical domain in cm
shock=findgen(nbre_iterations+1)
strom_length=Flux/(n_H*n_H*alpha_b)       ;stromgrem_length (cm)
strom_length_pc=strom_length/(3.09d18) ;stromgrem_length (pc)
print,'the stromgrem_length is :  ' ,strom_length ,'   cm'
time_s=findgen(nbre_iterations+1)*dt
time_ky=findgen(nbre_iterations+1)*dt/(3.16e10)
anal_front2=strom_length_pc*(1-exp(-time_s*alpha_B*n_H))
anal_front=(Flux*time_s)/(n_H*3.09e18)
anal_pos=strom_length_pc*(1+((1.25*time_s*c_sound)/(strom_length)))^(0.8)



for i=1,9 do begin   ;first 9 outputs    
    ;window, i
   read_bin_ath, 'no_HD_1e9.000'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    ; print,i0
    ;plot, d[*,0,0]   
    g=d_n(*,0,0)                                    
    h=d(*,0,0) ;dont call it E or it overwrites the energy!!!!
    y=(g/h)*1.0d0
    shock[i]=max(where(y lt 0.9 ))
;   print,i-1, ((shock[i-1]-anal_pos[i-1]))
 ;  print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
;print, total(d)
endfor


if nbre_iterations lt 100 then begin
last = nbre_iterations
endif  else begin
last = 99
endelse

for i=10,last do begin   ;next 99 outputs
    ;window, 
    read_bin_ath, 'no_HD_100.00'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    h=d(*,0,0)
    g=d_n(*,0,0)
    y=(g/h)*1.0d0
    shock[i]=max(where(y lt 0.9))
; print,i-1, ((shock[i-1]-anal_pos[i-1]))
; print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
  ;  print ,i,shock[i-1]
;print, total(d) ; check if  mass is conserved
endfor
;stop
shock=(shock*size_domain/(3.09e18))/(nbre_zones)



;plot,time_ky,shock, psym=-2
;oplot,time_ky, anal_front2

;window,1
;psopen,'R_front_1e9.ps'
set_plot,'ps'
device,filename='R_front_1e9.ps',/color,bits_per_pixel=8
loadct,13
!p.thick=2
!x.thick=2
!y.thick=2
Xsize=15
margin=0.12


plot,time_ky, shock,psym=-2,title='position of a R-front ',xtitle='time (kyr)', ytitle='front position (pc)',charsize=2.5,symsize=3
oplot,time_ky, anal_front,linestyle=1,color=20
oplot,[0,5],[0.0823,0.0823],linestyle=2,color=170
oplot,time_ky,anal_front2, color=255

oplot,[1.0,1.3],[0.016,0.016], psym=-2,symsize=3
xyouts, 1.5,0.016, 'simaultion',charsize=2.8
oplot, [1.0,1.3],[0.02,0.02],color=255
xyouts, 1.5,0.02, 'front with recombinations',charsize=2.8
oplot, [1.0,1.3],[0.024,0.024],linestyle=1,color=20
xyouts, 1.5,0.024, 'front without recombinations',charsize=2.8
oplot,[1.0,1.3],[0.028,0.028],linestyle=2,color=170
xyouts, 1.5,0.028, 'Stromgren length',charsize=2.8


;device,/close
;set_plot,'win'
psclose
end

