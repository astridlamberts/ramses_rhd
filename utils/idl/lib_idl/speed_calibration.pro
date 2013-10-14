;Astrid Lamberts
;program that plots the position if a shock in function of time and
;its analytic solution

nbre_iterations=35               ;number of outputs
nbre_zones=256                   ; number of cells
dt=5.0e12                        ;timestep between outputs (s)
Flux=1.0e8                       ;initial flux (photons per s per cm²)
c_sound=7.75e5                   ;sound speed in ionized region(cm/s)
n_H=double(100)                  ;initial number density (atoms per cm³)
alpha_b=2.59e-13                 ;recombination coefficient (cm³/s)
size_domain=3.09e19             ;size of the physical domain (cm)
shock=findgen(nbre_iterations)   ; array with shock position at each timestep  
time_s=findgen(nbre_iterations)*dt ;array with the physical time in s at each timestep
time_My=findgen(nbre_iterations)*dt/(3.16e13);array with the physical time in million years at each timestep


;extracting the shock position for each output

for i=1,9 do begin     ;first 9 outputs    
    read_bin_ath, 'no_turb_B_perp.000'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,bx=bx,by=by,bz=bz,/ion,/mhd,par =[4,1,1] 
   ; g=d_n(*,0,0)                                    
   ; h=d(*,0,0) 
   ; y=(g/h)*1.0d0
   ; shock[i-1]=max(where(  y    lt 0.99))
shock[i-1]=where(d[*,0,0] eq max (d[*,0,0]))
 ;   print, total(d) ; (to check if the mass is conserved)
endfor

;checking how many iterations there are
if nbre_iterations lt 100 then begin
   last = nbre_iterations
endif  else begin
   last = 99
endelse

for i=10,last do begin   ;next 99 outputs
    read_bin_ath, 'no_turb_B_perp.00'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,bx=bx,by=by,bz=bz,/ion,/mhd,par =[4,1,1]
   ; h=d(*,0,0)
  ;  g=d_n(*,0,0)
 ;   shock[i-1]=max(where(g/h lt 0.99 ))
shock[i-1]=where(d[*,0,0] eq max (d[*,0,0]))
  ;  print, total(d) ; check if  mass is conserved
endfor


if (nbre_iterations ge 100) then begin
   for i=100,(nbre_iterations) do begin   ;last outputs
      read_bin_ath, 'testcool5.0'+strtrim(i,2)+'.bin',d,px,py,pz,e=e,d_n=d_n,/ion;,par =[2,1,1]
      h=d(*,0,0)
      g=d_n(*,0,0)
      shock[i-1]=max(where((g/h) lt 0.95))
      ;print, total(d)
   endfor
end
; conversion into pc 
shock=(shock*size_domain)/(nbre_zones*3.09e18)

;computing the analytic solution
strom_length=Flux/(n_H*n_H*alpha_b)       ;stromgrem_length (cm)
strom_length_pc=strom_length/(3.09d18)    ;stromgrem_length (pc)
print,'the stromgrem_length is :  ' ,strom_length ,'   cm'

;anal_pos=strom_length_pc*(1+((1.25*time_s*c_sound)/(strom_length)))^(0.8)

plot,time_My,shock,psym=-2

;making the plots

;psopen, 'front_pos_1e9_100.ps'
;window,2
;!X.MARGIN=[10,3]
;!Y.MARGIN=[4,2]
;plot,time_My, shock, psym=-2,title='position of a D-front as a function of time ', xtitle='time (Myr)', ytitle='front position (pc)',charsize=2.8 , symsize=3;', /ylog
;oplot, time_My,anal_pos

;oplot, [0.2,0.4],[5,5]
;xyouts, 0.5,5, 'analytic solution',charsize=3.0
;oplot, [0.2,0.4],[5.5,5.5],psym=-2, symsize=3
;xyouts, 0.5,5.5, 'simulation',charsize=3.0

;window,1
;psclose
;psopen,'front_pos_1e8_100_dif.ps'
;!p.multi=[0,1,2]
;!X.MARGIN=[12,1]
;!Y.MARGIN=[3,3]

;plot,time_My,abs(anal_pos-shock),psym=-2,title='absolute difference ',xtitle='time(Myr)',ytitle='pc',charsize=3
;oplot,[0,40],[0.117,0.117],linestyle=3
;plot,time_My,abs(anal_pos-shock)/anal_pos,psym=-2,/ylog,title='relative difference',xtitle='time(Myr)',ytitle='relative difference',charsize=3

;psclose

end




;  read_bin_ath, 'no_turb_B.0001.bin',d,px,py,pz,e=e,d_n=d_n,bx=bx,by=by,bz=bz,/ion,/mhd,par =[4,1,1] 
