;program which computes the recombination rate at a given timestep. It
;takes into account the ;differences in temperatuure in each cell

size_domain=double(3.09e19)
nbre_cells=128
mu_H=double(2.34e-24) ;it sound like we don't need the mean gas mass per
;particle as we want the number of recombinations/s/cm²
limit= 0.5
i=0
recombination_rate=double(0.0)
size_cell=size_domain/nbre_cells

;print,'a'
while (limit lt 0.9999) do begin

    n_ion=double((d[i,0,0]-d_n[i,0,0])/mu_H)
;print,'b'
alpha_b=double(2.59e-13)
;    alpha_b=2.59e-13*(t[i,0,0]/1.0e4)^(-0.7)
;print,'c'
    recomb_rate=double(n_ion*n_ion*alpha_b)
;print,recomb_rate
    recombination_rate=double(recombination_rate+recomb_rate)
;print,'e'
    limit=double(d_n[i,0,0]/d[i,0,0])
    i=i+1
print,recombination_rate
print,i
endwhile
;stop
;to get the rate per cm²
recombination_rate=recombination_rate*size_cell

print,'the total recombination rate is',recombination_rate

end
