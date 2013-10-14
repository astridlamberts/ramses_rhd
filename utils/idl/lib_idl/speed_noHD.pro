;Astrid Lamberts
;program that plots the position if a shock in function of time and
;its analytic solution

nbre_iterations=20       ;number of outputs
nbre_zones=256           ; number of cells
dt=5e8                ;timestep
Flux=1.0e8              ;flux
Flux2=1.0e9
c_sound=7.80e5           ;sound speed in ionized region
n_H=double(100)          ;initial number density
alpha_b=3.9e-13
;3.89e-13         ;recombination coefficient
size_domain=2*3.09e18     ;size of the physical domain in cm
size_domain2=2*3.09e17   
shock=findgen(nbre_iterations+1)   ; array with shock position at each timestep  
shock2=findgen(nbre_iterations+1)   ; array with shock position at each timestep  

strom_length=Flux/(n_H*n_H*alpha_b)       ;stromgrem_length (cm)
strom_length_pc=strom_length/(3.09d18) ;stromgrem_length (pc)
print,'the stromgrem_length is :  ' ,strom_length ,'   cm'
time_s=findgen(nbre_iterations+1)*dt
time_My=findgen(nbre_iterations+1)*dt/(3.16e10)
;anal_pos=strom_length*(1+((1.25*time_s*c_sound)/(strom_length)))^(0.8)
; determination of shock position for all timesteps
; in the first timestep there is no shock
shock[0]=0.0
shock2[0]=0.0
;anal_front=(Flux*time_s)/n_H

for i=1,9 do begin   ;first 9 outputs    
    ;window, i
   read_bin_ath, 'no_HD_norec_1e10.000'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    ; print,i0
    ;plot, d[*,0,0]   
    g=d_n(*,0,0)                                    
    h=d(*,0,0) ;dont call it E or it overwrites the energy!!!!
    y=(g/h)*1.0d0
    shock[i]=max(where(y lt 0.99 ))
;   print,i-1, ((shock[i-1]-anal_pos[i-1]))
 ;  print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
;print, total(d)
endfor

;print, 'a'
if nbre_iterations lt 100 then begin
last = nbre_iterations
endif  else begin
last = 99
endelse

for i=10,last do begin   ;next 99 outputs
    ;window, 
    read_bin_ath, 'no_HD_norec_1e10.00'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    h=d(*,0,0)
    g=d_n(*,0,0)
    shock[i]=max(where((g/h) lt 0.99))
; print,i-1, ((shock[i-1]-anal_pos[i-1]))
; print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
  ;  print ,i,shock[i-1]
;print, total(d) ; check if  mass is conserved
endfor
;print,' b'

if (nbre_iterations ge 100) then begin
for i=100,(nbre_iterations) do begin   ;other outputs
    ;window, i
   read_bin_ath, 'no_HD2_256_1e10_100.0'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
  ;  plot, d[*,0,0]                                           
  ; print,x
    h=d(*,0,0)
g=d_n(*,0,0)
 ;   x=double(max(f))
    shock[i]=max(where((g/h) lt 0.95))
 ;print,i-1,( (shock[i-1]-anal_pos[i-1]))
  ;   print ,i,shock[i-1]
; print,i,((shock[i-1]*2*3.09d19)/(nbre_zones*3.09d18))
;print, total(d)
endfor
endif

for i=1,9 do begin   ;first 9 outputs    
    ;window, i
   read_bin_ath, 'no_HD_norec_1e9.000'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    ; print,i0
    ;plot, d[*,0,0]   
    g=d_n(*,0,0)                                    
    h=d(*,0,0) ;dont call it E or it overwrites the energy!!!!
    y=(g/h)*1.0d0
    shock2[i]=max(where(y lt 0.99 ))
;   print,i-1, ((shock[i-1]-anal_pos[i-1]))
 ;  print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
;print, total(d)
endfor

;print, 'a'
if nbre_iterations lt 100 then begin
last = nbre_iterations
endif  else begin
last = 99
endelse

for i=10,last do begin   ;next 99 outputs
    ;window, 
    read_bin_ath, 'no_HD_norec_1e9.00'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    h=d(*,0,0)
    g=d_n(*,0,0)
    shock2[i]=max(where((g/h) lt 0.99))
; print,i-1, ((shock[i-1]-anal_pos[i-1]))
; print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
  ;  print ,i,shock[i-1]
;print, total(d) ; check if  mass is conserved
endfor
;print,' b'

if (nbre_iterations ge 100) then begin
for i=100,(nbre_iterations) do begin   ;other outputs
    ;window, i
   read_bin_ath, 'no_HD2_256_1e10_100.0'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
  ;  plot, d[*,0,0]                                           
  ; print,x
    h=d(*,0,0)
g=d_n(*,0,0)
 ;   x=double(max(f))
    shock[i]=max(where((g/h) lt 0.95))
 ;print,i-1,( (shock[i-1]-anal_pos[i-1]))
  ;   print ,i,shock[i-1]
; print,i,((shock[i-1]*2*3.09d19)/(nbre_zones*3.09d18))
;print, total(d)
endfor
endif

; conversion into pc 
shock2=(shock2*size_domain2/(3.09e18))/(nbre_zones)
shock=(shock*size_domain/(3.09e18))/(nbre_zones)
;fit=linfit(time_ky,shock)
;y=fit[0]+time_ky*fit[1]
;for i=1, nbre_iterations-1 do begin
;print, i, shock[i]-anal_pos[i]
;endfor

;for the recombination case 
;analytic solution for the front position as a function of time
;t_rec=1.0/(alpha_B*n_H)
;anal_front2=strom_length_pc*(1-exp(-time_s*alpha_B*n_H))
anal_front=(Flux*time_s)/(n_H*3.09e18)
anal_front2=(Flux2*time_s)/(n_H*3.09e18)
;anal_pos=strom_length*(1+((1.25*time_s*c_sound)/(strom_length)))^(0.8)
;analytic curve

psopen, 'test.ps'


;set_plot,'ps'
!p.thick=2
!x.thick=2
!y.thick=2
;!p.multi=[0,2,0]
;Xsize=15
;!x.margin=[8,4]
;!y.margin=[6,6]

;plot,time_ky,anal_front,title='Position of R-front without recombinations ',xtitle='time (kyr)', ytitle='front position (pc)',charsize=2.3
;plot,time_ky, anal_front2,title='position of a R-front ',xtitle='time (kyr)', ytitle='front position (pc)',charsize=2.3
;oplot,time_ky, shock, psym=-4,linestyle=1
;oplot,linfit(time_ky,shock),linestyle=1
;oplot,[0.02,0.05],[0.25,0.25]
;xyouts, 0.07,0.25, 'analytic solution',charsize=2.8
;oplot, [0.02,0.05],[0.265,0.265],psym=-2, linestyle=1
;xyouts, 0.07,0.265, 'simulation',charsize=2.8

;plot,time_ky,shock2,psym=-2,title='Position of R-front F=10^9 ',xtitle='time (kyr)', ytitle='front position (pc)',charsize=2.3, linestyle=1

;oplot,time_ky, anal_front2

;oplot,[0.02,0.05],[0.025,0.025]
;xyouts, 0.07,0.025, 'analytic solution',charsize=2.8
;oplot, [0.02,0.05],[0.0265,0.0265],psym=-2, linestyle=1
;xyouts, 0.07,0.0265, 'simulation',charsize=2.8




psclose

;window,3

;test=abs(y-(shock-0.00248))/y
;test1=anal_front[1:20]
;time_test=time_ky[0:20]
;window,1
;oplot,time_s, anal_front
;psclose
;psopen,'first_stage_1e10_100_dif.ps'
;!p.multi=[0,1,2];
;!X.MARGIN=[12,1]
;!Y.MARGIN=[3,3];
;plot,[0,0.3],[3.9e-3,3.9e-3],linestyle=3,title='absolute difference ',xtitle='time(kyr)',ytitle='pc',charsize=3
;oplot,time_ky,abs(shock-anal_front),psym=-2;,title='absolute difference ',xtitle='time(kyr)',ytitle='pc',charsize=3
;oplot,[0,0.3],[7.8e-4,7.8e-4],linestyle=3
;plot,time_test,abs(test1-shock[1:20])/test1,psym=-2,/ylog,title='relative difference',xtitle='time(kyr)',ytitle='relative difference',charsize=3

;psclose
;stop
end

;.compile athtemp
;t=athtemp (d,px,py,pz,e,d_n=d_n)

;krumholz/lib/mpich-1.2.7/bin/mpirun -np 2 athena -i athinput    
;read_bin_ath, 'cool1e-2_noHD.0001.bin',d,px,py,pz,e=e,d_n=d_n,/ion,par =[2,1,1]
;device, true_color=24,retain=2

;Astrid Lamberts
;program that plots the position if a shock in function of time and
;its analytic solution


