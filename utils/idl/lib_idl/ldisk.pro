function ldisk, m, mdot, r

fk = 0.5
g=6.67d-8
return, (1.0-fk)*g*m*mdot/r

end
