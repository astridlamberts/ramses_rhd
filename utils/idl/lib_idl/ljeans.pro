function ljeans, rho, t

kb=1.38e-16
mp=1.67e-24
mu=14./6.
G=6.67e-8
pi=3.1415927
cs=sqrt(kb*t/(mu*mp))
return, cs*sqrt(pi/(G*rho))

end
