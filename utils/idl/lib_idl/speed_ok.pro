;Astrid Lamberts
;program that plots the position if a shock in function of time and
;its analytic solution

nbre_iterations=40
 ;number of outputs
nbre_zones=1024   ; number of cells
dt=5.0e12         
   ;timestep
Flux=3.0e8                ;flux
c_sound=7.758d5            ;sound speed in ionized region
n_H=double(500)          ;initial number density
alpha_b=4.146d-13          ;recombination coefficient
size_domain=4*3.09d19     ;size of the physical domain in cm
shock=findgen(nbre_iterations+1)   ; array with shock position at eacxh timestep  

time_s=findgen(nbre_iterations+1)*dt
anal_front=(Flux*time_s)/n_H
; determination of shock position for all timesteps
; in the first timestep there is no shock
shock[0]=0.0

for i=1,9 do begin   ;first 9 outputs    
    ;window, i
    read_bin_ath, 'test_check_cel7.000'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    ; print,i0
    ;plot, d[*,0,0]   
    g=d_n(*,0,0)                                    
    h=d(*,0,0) ;dont call it E or it overwrites the energy!!!!
    y=(g/h)*1.0d0
    shock[i]=max(where(y lt 0.99 ))
  
   print,i, (shock[i]-anal_front[i])
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
    read_bin_ath, 'test_check_cel7.00'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
    h=d(*,0,0)
    g=d_n(*,0,0)
    shock[i]=max(where((g/h) lt 0.99))
; print,i,((shock[i-1]*2*3.09d19)/nbre_zones*3.09d18)
  ;  print ,i,shock[i-1]
;print, total(d) ; check if  mass is conserved
 print,i, (shock[i]-anal_front[i])
endfor
;print,' b'

if (nbre_iterations ge 100) then begin
for i=100,(nbre_iterations) do begin   ;other outputs
    ;window, i
   read_bin_ath, 'test_check_cel7.0'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
  ;  plot, d[*,0,0]                                           
  ; print,x
    h=d(*,0,0)
g=d_n(*,0,0)
 ;   x=double(max(f))
    shock[i]=max(where((g/h) lt 0.95))
 print,i, (shock[i]-anal_front[i])
  ;   print ,i,shock[i-1]
; print,i,((shock[i-1]*2*3.09d19)/(nbre_zones*3.09d18))
;print, total(d)
endfor
endif


; conversion into cm 
shock=(shock*size_domain)/(nbre_zones)

;linear fit to find mean velocity
;tim=findgen(nbre_iterations)*dt
velocity=linfit(time_s,shock)
;plot,shock,xtitle='position (cells)',ytitle ='time(timesteps=2.0e11)',title='position of the shock versus time'
;psclose
print,'the mean velocity is' ,velocity[1] ,'  cm per s'


strom_length=Flux/(n_H*n_H*alpha_b)       ;stromgrem_length (cm)
strom_length_pc=strom_length/(3.09d18) ;stromgrem_length (pc)
print,'the stromgrem_length is :  ' ,strom_length ,'   cm'

;analytic solution for the shock position as a function of time
;anal_front=(Flux*time_s)/n_H
anal_pos=strom_length*(1+((1.25*time_s*c_sound)/(strom_length)))^(0.8)
;analytic curve
;window,2
;psopen, 'front_pos_3e9_1000_512_80.ps'
window,4
plot,time_s, shock, psym=-2;,title='front position for the  analytic and simulated case (*) with F= 3e9 n_h=1000 and 512 cells', xtitle='time (s)', ytitle='front position (cm)' ;', /ylog
;oplot, time_s,anal_pos
oplot, time_s,anal_front;,title='shock position for analytic and simulated case (*)', xtitle='time in timesteps', ytitle='position in parsecs' 

;psclose
;window,3
;plot, time_s,shock, psym=-2
;oplot, time_s,anal_pos

;psopen, 'shock_pos_3-4_dt.ps'
;plot, time_s,anal_pos,title='1e-2'
;plot, time_s,shock, psym=-2,title='3e-04'
;plot, time_s,anal_pos
;psclose

;shock speed computation
;b=(strom_length_pc)^(-0.5)
;position=findgen(300)*0.1+0.00001
;vs=c_sound/sqrt((sqrt(position)*b -0.25))
;plot,vs

;for j=1,15 do begin
;i=(j)*20  ;finds which timestep it is
;v_shock[j]=(2*3.09d19*(shock[i+10]-shock[i-10]))/(20*dt*512)
;;print,j,shock[j],v_shock[j]
;po[j]=shock[i]
;pos[j]=po[j]*((3.09d19*2)/(512*3.09d17))
;endfor
;v_shock[0]=(2*3.09d19*(shock[10]))/(10*dt*512)
;pos[0]=0

stop
end

;.compile athtemp
;t=athtemp (d,px,py,pz,e,d_n=d_n)

;krumholz/lib/mpich-1.2.7/bin/mpirun -np 2 athena -i athinput    
;read_bin_ath, 'cool1e-4___cs_noT.0001.bin',d,px,py,pz,e=e,d_n=d_n,/ion,par =[2,1,1]
;device, true_color=24,retain=2

;Astrid Lamberts
;program that plots the position if a shock in function of time and
;its analytic solution


